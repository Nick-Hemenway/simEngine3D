import numpy as np
from scipy.misc import derivative as deriv

def tilde(v):
    
    x,y,z = v
    v_tilde = np.array([[0, -z,  y],
                        [z,  0, -x],
                        [-y, x,  0]])
    
    return v_tilde

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
        
        if constraint_func == None:
            self.func = ConstraintFunc()
        else:
            self.func = constraint_func
    
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
    
class GconDP1(Gcon):
    
    def __init__(self, body_i, ai_bar, body_j, aj_bar, constraint_func = None):
        
        super().__init__(body_i, body_j, constraint_func) #pass arguments to 
                                                          #parent class
                                                          #constructor
        #store both of the local coordinate vectors ensuring they are both
        #numpy arrays
        self.ai_bar = np.atleast_2d(ai_bar).reshape(-1,1)
        self.aj_bar = np.atleast_2d(aj_bar).reshape(-1,1)
        
    def val(self, t):
        
        Ai = self.body_i.A
        Aj = self.body_j.A
        
        return self.ai_bar.T @ Ai.T @ Aj @ self.aj_bar - self.func.val(t)
    
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
        
        return gamma
    
class GconDP2(Gcon):
    
    def __init__(self):
        
        pass
    
    
class GconD(Gcon):
    
    def __init__(self):
        
        pass
    

class GconCD(Gcon):
    
    def __init__(self, body_i, sp_i_bar, body_j, sq_j_bar, c_vec, constraint_func = None):
        
        super().__init__(body_i, body_j, constraint_func) #pass arguments to 
                                                          #parent class
                                                          #constructor
        #store the local vectors ensuring they are numpy column vectors
        self.sp_i_bar = np.atleast_2d(sp_i_bar).reshape(-1,1)
        self.sq_j_bar = np.atleast_2d(sq_j_bar).reshape(-1,1)
        self.c_vec = np.atleast_2d(c_vec).reshape(-1,1)
        
    def val(self, t):
        
        Ai = self.body_i.A
        Aj = self.body_j.A
        
        rp = self.body_i.r + Ai @ self.sp_i_bar
        rq = self.body_j.r + Aj @ self.sp_j_bar
        
        dij = rq - rp
        
        return self.c_vec.T @ dij - self.func.val(t)
    
    def vel_rhs(self, t):
        
        return self.func.deriv1(t)
    
    def accel_rhs(self, t):
        
        #extract vectors into local variables for code readability
        pi_dot = self.body_i.p_dot
        pj_dot = self.body_j.p_dot
        
        t1 =  self.c_vec.T @ self.B(pi_dot, self.sp_i_bar) @ pi_dot
        t2 = -self.c_vec.T @ self.B(pj_dot, self.sq_j_bar) @ pj_dot
        
        gamma = t1 + t2 + self.func.deriv2(t)
        
        return gamma
    



