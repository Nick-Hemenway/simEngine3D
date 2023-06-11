import numpy as np
from simengine3D.geometric_constraints import primitives as prim
from simengine3D.utility import tilde, column

class DerivedConstraint():
    
    def __init__(self, gcon_list, body_i, body_j, sp_i_bar = None, sq_j_bar =  None):
        
        self.gcon_list = gcon_list
        
        self.body_i = body_i
        self.body_j = body_j
        
        self.i_idx = body_i.idx #index of body i in bodies list at system level
        self.j_idx = body_j.idx #index of body j in bodies list at system level
        
        #the following attributes are the location of the joint relative to each
        #body's center of gravity. It is used for computing reaction forces 
        #about the joint. If no locations are given, the vectors default to zero
        #vectors and the 'reaction()' method will simply return the reactons 
        #about the center of gravity.
        if sp_i_bar == None: 
            #if a location for the joint
            sp_i_bar = np.zeros((3,1))
            sq_j_bar = np.zeros((3,1))
        else:
            self.sp_i_bar = column(sp_i_bar) #location of joint on body i
            self.sq_j_bar = column(sq_j_bar) #location of joint on body j
        
        self.lagrange = None
        
    def val(self, t):
        
        vals = []
        for gcon in self.gcon_list:
            vals.append(gcon.val(t))
            
        vals = np.vstack(vals)
        
        return vals
    
    def vel_rhs(self, t):
        
        nu = []
        for gcon in self.gcon_list:
            nu.append(gcon.vel_rhs(t))
        nu = np.vstack(nu)
        
        return nu
    
    def accel_rhs(self, t):
        
        gamma =  []
        for gcon in self.gcon_list:
            gamma.append(gcon.accel_rhs(t))
        gamma = np.vstack(gamma)
        
        return gamma
        
    def partial_p(self):
        
        dphi_dpi = []
        dphi_dpj = []
        for gcon in self.gcon_list:
            t1, t2 = gcon.partial_p()
            dphi_dpi.append(t1)
            dphi_dpj.append(t2)
        
        dphi_dpi = np.vstack(dphi_dpi)
        dphi_dpj = np.vstack(dphi_dpj)
        
        return dphi_dpi, dphi_dpj
    
    def partial_r(self):
        
        dphi_dri = []
        dphi_drj = []
        for gcon in self.gcon_list:
            t1, t2 = gcon.partial_r()
            dphi_dri.append(t1)
            dphi_drj.append(t2)
        
        dphi_dri = np.vstack(dphi_dri)
        dphi_drj = np.vstack(dphi_drj)
        
        return dphi_dri, dphi_drj
    
    def pi(self):
        
        dphi_dpi, dphi_dpj = self.partial_p()
        
        #slide 16 of lecture 9 for going from phi to pi
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
    
    def reaction(self, body = 'i'):
        
        """Calculate reaction forces and torques at the joint location
        
        The reaction forces and torques are computed about the center of gravity
        automatically but it would be more useful to translate the torques to 
        to actually be about the joint so that we could use it to determine the
        torque required by a motor. To do this we do the following:
            
        Tg = rxF + T_joint --> T_joint = Tg - r x F
        
        Where Tg is the net torque that the contraint creates about the center 
        of gravity and is created by a combination of a moment about the joint
        plus the effects of the reaction for F being applied at the joint. To 
        determine the effects of just the moment at the joint, we must subtract
        off the contribution of the force at some lever arm from the net torque
        """
        
        #calculate reaction forces about center of gravity
        Tg = self.reaction_torque(body) #net torque about center of gravity
        Fg = self.reaction_force(body) #net force located at center of gravity
        
        if body.lower() == 'j':
            r = self.body_j.A @ self.sq_j_bar
        else:
            r = self.body_i.A @ self.sp_i_bar
        
        T_joint = Tg - tilde(r) @ Fg
        
        return T_joint, Fg
        
        
class GconPerp1(DerivedConstraint):
    
    """The perpendicular 1 constraint is discussed on slide 32 of lecture 8"""
    
    def __init__(self, body_i, ai_bar, bi_bar, body_j, cj_bar):
        
        DP1_1 = prim.GconDP1(body_i, ai_bar, body_j, cj_bar, constraint_func = None)
        DP1_2 = prim.GconDP1(body_i, bi_bar, body_j, cj_bar, constraint_func = None)
        
        super().__init__([DP1_1, DP1_2], body_i, body_j)
        
        self.DOF_constrained = 2 #constrains 2 DOF
        
    
class GconPerp2(DerivedConstraint):
    
    """The perpendicular 2 constraint is discussed on slide 34 of lecture 8"""
    
    def __init__(self, body_i, ai_bar, bi_bar, sp_i_bar, body_j, sq_j_bar):
        
        DP2_1 = prim.GconDP2(body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None)
        DP2_2 = prim.GconDP2(body_i, bi_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None)
        
        super().__init__([DP2_1, DP2_2], body_i, body_j)
        
        self.DOF_constrained = 2 #constrains 2 DOF
    
    
    