import numpy as np
from scipy.misc import derivative as deriv
from simengine3D.utility import tilde, column

class ConstraintFunc():
    """All of the constraints contain an f(t) term. This class stores the 
    function f(t) and its derivatives. If no functions are given as input
    arguments, the functions automatically default to returning zero"""
    
    def __init__(self, f = None, f_prime = None, f_pprime = None):
        
        #if no function is provided set f and its derivatives to return zero
        if f == None: 
            zero_func = lambda t: 0
            
            self.f = zero_func
            self.f_prime = zero_func
            self.f_pprime = zero_func
            
        else:
            self.f = f
            
            if f_prime == None:
                #numerically differentiate f if a function isn't given
                self.f_prime = lambda t: deriv(f, t, dx = 1e-6, n = 1)
            else:
                self.f_prime = f_prime
                
            if f_pprime == None:
                #numerically take second derivative of f if function isn't given
                self.f_pprime = lambda t: deriv(f, t, dx = 1e-6, n = 2)
            else:
                self.f_pprime = f_pprime
        
    def val(self, t):
        return self.f(t) #return value of f(t)
    
    def deriv1(self, t):
        return self.f_prime(t) #return first derivative of f(t) at t
    
    def deriv2(self, t):
        return self.f_pprime(t) #return second derivative of f(t) at t
    

class Gcon():
    """Parent class for all of the constraint primitives"""
    
    def __init__(self, body_i, body_j, constraint_func = None):
        """body_i and body_j are both Body objects and constraint_func is
        is a ConstraintFunc object"""
        
        self.body_i = body_i
        self.body_j = body_j
        
        self.i_idx = self.body_i.idx
        self.j_idx = self.body_j.idx
        
        if constraint_func == None:
            self.func = ConstraintFunc()
        else:
            self.func = constraint_func
            
        self.DOF_constrained = 1 #all primitive constraints constrain 1 DOF
        
        self.lagrange = None
    
    def B(self, p, a_bar):
        """Returns 3x4 matrix that is created by the B operator.
                p = the time varying vector
                a_bar = the constant local vector"""  
                
        e0 = p[0]
        e = p[1::]        
        e_tilde = tilde(e)
        a_bar_tilde = tilde(a_bar)
        
        t0 = e0*np.eye(3) + e_tilde
        
        v1 = t0 @ a_bar
        m1 = e @ a_bar.T - t0 @ a_bar_tilde
        
        M = 2*np.hstack((v1, m1))
        
        return M
    
    def dij_calc(self, sp_i_bar, sq_j_bar):
        
        rp = self.body_i.r + self.body_i.A @ sp_i_bar
        rq = self.body_j.r + self.body_j.A @ sq_j_bar
        
        dij = rq - rp
        
        return dij
    
class GconDP1(Gcon):
    
    def __init__(self, body_i, ai_bar, body_j, aj_bar, constraint_func = None):
        
        super().__init__(body_i, body_j, constraint_func) #pass arguments to 
                                                          #parent class
                                                          #constructor
        #store both of the local coordinate vectors ensuring they are both
        #numpy arrays
        self.ai_bar = column(ai_bar)
        self.aj_bar = column(aj_bar)
        
    def val(self, t):
        
        Ai = self.body_i.A
        Aj = self.body_j.A
        
        return (self.ai_bar.T @ Ai.T @ Aj @ self.aj_bar - self.func.val(t)).item()
    
    def vel_rhs(self, t):
        
        return self.func.deriv1(t)
    
    def accel_rhs(self, t):
        
        #extract vectors into local variables for code readability
        pi = self.body_i.p
        pj = self.body_j.p
        pi_dot = self.body_i.p_dot
        pj_dot = self.body_j.p_dot
        
        ai = self.body_i.A @ self.ai_bar
        aj = self.body_j.A @ self.aj_bar
        
        ai_dot = self.B(pi, self.ai_bar) @ pi_dot
        aj_dot = self.B(pj, self.aj_bar) @ pj_dot
        
        t1 = -ai.T @ self.B(pj_dot, self.aj_bar) @ pj_dot
        t2 = -aj.T @ self.B(pi_dot, self.ai_bar) @ pi_dot
        t3 = -2*ai_dot.T @ aj_dot
        
        gamma = t1 + t2 + t3 + self.func.deriv2(t)
        
        return gamma.item()
    
    def partial_p(self):
        
        pi = self.body_i.p
        pj = self.body_j.p
        Ai = self.body_i.A
        Aj = self.body_j.A
        
        ai = Ai @ self.ai_bar
        aj = Aj @ self.aj_bar
        
        dphi_dpi = aj.T @ self.B(pi, self.ai_bar)
        dphi_dpj = ai.T @ self.B(pj, self.aj_bar)
        
        return dphi_dpi.flatten(), dphi_dpj.flatten()
    
    def partial_r(self):
        
        dphi_dri = np.zeros(3)
        dphi_drj = np.zeros(3)
        
        return dphi_dri, dphi_drj
    
    def pi(self):
        
        dphi_dpi, dphi_dpj = self.partial_p()
        
        pi_i = 0.5*dphi_dpi @ self.body_i.E.T
        pi_j = 0.5*dphi_dpj @ self.body_j.E.T
        
        return pi_i, pi_j
    
    def reaction_force(self, body = 'i'):
        
        #note that reaction forces act at center of gravity, not joint
        phi_ri, phi_rj = self.partial_r()
        
        if body.lower() == 'j':
            partial = np.atleast_2d(phi_rj)
        else:
            partial = np.atleast_2d(phi_ri)
            
        F = -partial.T @ self.lagrange
        
        return F
    
    def reaction_torque(self, body = 'i'):
        
        #note that reaction torques act at center of gravity, not joint
        Pi_i, Pi_j = self.pi()
        
        if body.lower() == 'j':
            Pi = np.atleast_2d(Pi_j)
        else:
            Pi = np.atleast_2d(Pi_i)
            
        T = -Pi.T @ self.lagrange
        
        return T
        
    
class GconDP2(Gcon):
    
    def __init__(self, body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None):
        
        super().__init__(body_i, body_j, constraint_func) #pass arguments to 
                                                          #parent class
                                                          #constructor
        #store both of the local coordinate vectors ensuring they are both
        #numpy arrays
        self.ai_bar =   column(ai_bar)
        self.sp_i_bar = column(sp_i_bar)
        self.sq_j_bar = column(sq_j_bar)
        
    def val(self, t):
        
        Ai = self.body_i.A
        dij = self.dij_calc(self.sp_i_bar, self.sq_j_bar)
        
        return (self.ai_bar.T @ Ai.T @ dij - self.func.val(t)).item()
    
    def vel_rhs(self, t):
        
        return self.func.deriv1(t)
    
    def accel_rhs(self, t):
        
        #extract vectors into local variables for code readability
        pi = self.body_i.p
        pj = self.body_j.p
        pi_dot = self.body_i.p_dot
        pj_dot = self.body_j.p_dot
        ri_dot = self.body_i.r_dot
        rj_dot = self.body_j.r_dot
        
        Ai = self.body_i.A
        ai = Ai @ self.ai_bar
        ai_dot = self.B(pi, self.ai_bar) @ pi_dot
        
        dij = self.dij_calc(self.sp_i_bar, self.sq_j_bar)
        
        term1 =  rj_dot + self.B(pj, self.sq_j_bar) @ pj_dot
        term2 = -ri_dot - self.B(pi, self.sp_i_bar) @ pi_dot
        dij_dot = term1 + term2
        
        t1 = -ai.T @ self.B(pj_dot, self.sq_j_bar) @ pj_dot
        t2 =  ai.T @ self.B(pi_dot, self.sp_i_bar) @ pi_dot
        t3 = -dij.T @ self.B(pi_dot, self.ai_bar) @ pi_dot
        t4 = -2*ai_dot.T @ dij_dot
        
        gamma = t1 + t2 + t3 + t4 + self.func.deriv2(t)
        
        return gamma.item()
    
    def partial_p(self):
        #comes from slide 22 of lecture 9
        pi = self.body_i.p
        pj = self.body_j.p
        
        Ai = self.body_i.A
        ai = Ai @ self.ai_bar
        
        dij = self.dij_calc(self.sp_i_bar, self.sq_j_bar)
        
        dphi_dpi = dij.T @ self.B(pi, self.ai_bar) - ai.T @ self.B(pi, self.sp_i_bar)
        dphi_dpj = ai.T @ self.B(pj, self.sq_j_bar)
        
        return dphi_dpi.flatten(), dphi_dpj.flatten()
        
    
    def partial_r(self):
        #comes from slide 22 of lecture 9
        
        Ai = self.body_i.A
        ai = Ai @ self.ai_bar
        
        dphi_dri = -ai.T
        dphi_drj =  ai.T
        
        return dphi_dri.flatten(), dphi_drj.flatten()
    
    def pi(self):
        
        dphi_dpi, dphi_dpj = self.partial_p()
        
        pi_i = 0.5*dphi_dpi @ self.body_i.E.T
        pi_j = 0.5*dphi_dpj @ self.body_j.E.T
        
        return pi_i, pi_j
    
    def reaction_force(self, body = 'i'):
        
        #note that reaction forces act at center of gravity, not joint
        phi_ri, phi_rj = self.partial_r()
        
        if body.lower() == 'j':
            partial = np.atleast_2d(phi_rj)
        else:
            partial = np.atleast_2d(phi_ri)
            
        F = -partial.T @ self.lagrange
        
        return F
    
    def reaction_torque(self, body = 'i'):
        
        #note that reaction torques act at center of gravity, not joint
        Pi_i, Pi_j = self.pi()
        
        if body.lower() == 'j':
            Pi = np.atleast_2d(Pi_j)
        else:
            Pi = np.atleast_2d(Pi_i)
            
        T = -Pi.T @ self.lagrange
        
        return T
        
        
class GconD(Gcon):
    
    def __init__(self, body_i, sp_i_bar, body_j, sq_j_bar, constraint_func = None):
        
        super().__init__(body_i, body_j, constraint_func) #pass arguments to 
                                                          #parent class
                                                          #constructor
        #store the local vectors ensuring they are numpy column vectors
        self.sp_i_bar = column(sp_i_bar)
        self.sq_j_bar = column(sq_j_bar)
        
    def val(self, t):
        
        dij = self.dij_calc(self.sp_i_bar, self.sq_j_bar)
        
        return (dij.T @ dij - self.func.val(t)).item()
    
    def vel_rhs(self, t):
        
        return self.func.deriv1(t)
    
    def accel_rhs(self, t):
        
        #extract vectors into local variables for code readability
        pi = self.body_i.p
        pj = self.body_j.p
        pi_dot = self.body_i.p_dot
        pj_dot = self.body_j.p_dot
        ri_dot = self.body_i.r_dot
        rj_dot = self.body_j.r_dot
        
        dij = self.dij_calc(self.sp_i_bar, self.sq_j_bar)
        
        term1 =  rj_dot + self.B(pj, self.sq_j_bar) @ pj_dot
        term2 = -ri_dot - self.B(pi, self.sp_i_bar) @ pi_dot
        dij_dot = term1 + term2
        
        t1 =  -2*dij.T @ self.B(pj_dot, self.sq_j_bar) @ pj_dot
        t2 =   2*dij.T @ self.B(pi_dot, self.sp_i_bar) @ pi_dot
        t3 = -2*dij_dot.T @ dij_dot
        
        gamma = t1 + t2 + t3 + self.func.deriv2(t)
        
        return gamma.item()
    
    def partial_p(self):
        
        #comes from slide 23, lecture 9
        pi = self.body_i.p
        pj = self.body_j.p
        
        dij = self.dij_calc(self.sp_i_bar, self.sq_j_bar)
        
        dphi_dpi = -2*dij.T @ self.B(pi, self.sp_i_bar)
        dphi_dpj =  2*dij.T @ self.B(pj, self.sq_j_bar)
        
        return dphi_dpi.flatten(), dphi_dpj.flatten()
        
    def partial_r(self):

        #comes from slide 23, lecture 9
        
        dij = self.dij_calc(self.sp_i_bar, self.sq_j_bar)
        
        dphi_dri = -2*dij.T
        dphi_drj =  2*dij.T
        
        return dphi_dri.flatten(), dphi_drj.flatten()
    
    def pi(self):
        
        dphi_dpi, dphi_dpj = self.partial_p()
        
        pi_i = 0.5*dphi_dpi @ self.body_i.E.T
        pi_j = 0.5*dphi_dpj @ self.body_j.E.T
        
        return pi_i, pi_j
    
    def reaction_force(self, body = 'i'):
        
        #note that reaction forces act at center of gravity, not joint
        phi_ri, phi_rj = self.partial_r()
        
        if body.lower() == 'j':
            partial = np.atleast_2d(phi_rj)
        else:
            partial = np.atleast_2d(phi_ri)
            
        F = -partial.T @ self.lagrange
        
        return F
    
    def reaction_torque(self, body = 'i'):
        
        #note that reaction torques act at center of gravity, not joint
        Pi_i, Pi_j = self.pi()
        
        if body.lower() == 'j':
            Pi = np.atleast_2d(Pi_j)
        else:
            Pi = np.atleast_2d(Pi_i)
            
        T = -Pi.T @ self.lagrange
        
        return T

class GconCD(Gcon):
    
    def __init__(self, body_i, sp_i_bar, body_j, sq_j_bar, c_vec, constraint_func = None):
        
        super().__init__(body_i, body_j, constraint_func) #pass arguments to 
                                                          #parent class
                                                          #constructor
        #store the local vectors ensuring they are numpy column vectors
        self.sp_i_bar = column(sp_i_bar)
        self.sq_j_bar = column(sq_j_bar)
        self.c_vec =    column(c_vec)
        
    def val(self, t):
        
        dij = self.dij_calc(self.sp_i_bar, self.sq_j_bar)
        
        return (self.c_vec.T @ dij - self.func.val(t)).item()
    
    def vel_rhs(self, t):
        
        return self.func.deriv1(t)
    
    def accel_rhs(self, t):
        
        #extract vectors into local variables for code readability
        pi_dot = self.body_i.p_dot
        pj_dot = self.body_j.p_dot
        
        t1 =  self.c_vec.T @ self.B(pi_dot, self.sp_i_bar) @ pi_dot
        t2 = -self.c_vec.T @ self.B(pj_dot, self.sq_j_bar) @ pj_dot
        
        gamma = t1 + t2 + self.func.deriv2(t)
        
        return gamma.item()
    
    def partial_p(self):
        
        pi = self.body_i.p
        pj = self.body_j.p
        
        c = self.c_vec
        
        dphi_dpi = -c.T @ self.B(pi, self.sp_i_bar)
        dphi_dpj =  c.T @ self.B(pj, self.sq_j_bar)
        
        return dphi_dpi.flatten(), dphi_dpj.flatten()
    
    def partial_r(self):
        
        dphi_dri = -self.c_vec.T
        dphi_drj =  self.c_vec.T
        
        return dphi_dri.flatten(), dphi_drj.flatten()
    
    def pi(self):
        
        dphi_dpi, dphi_dpj = self.partial_p()
        
        pi_i = 0.5*dphi_dpi @ self.body_i.E.T
        pi_j = 0.5*dphi_dpj @ self.body_j.E.T
        
        return pi_i, pi_j
    
    def reaction_force(self, body = 'i'):
        
        #note that reaction forces act at center of gravity, not joint
        phi_ri, phi_rj = self.partial_r()
        
        if body.lower() == 'j':
            partial = np.atleast_2d(phi_rj)
        else:
            partial = np.atleast_2d(phi_ri)
            
        F = -partial.T @ self.lagrange
        
        return F
    
    def reaction_torque(self, body = 'i'):
        
        #note that reaction torques act at center of gravity, not joint
        Pi_i, Pi_j = self.pi()
        
        if body.lower() == 'j':
            Pi = np.atleast_2d(Pi_j)
        else:
            Pi = np.atleast_2d(Pi_i)
            
        T = -Pi.T @ self.lagrange
        
        return T
        
