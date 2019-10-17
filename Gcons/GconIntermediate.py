#add the parent folder to the search path for importing modules
import sys
sys.path.append('..')

import numpy as np
from Gcons import GconPrimitives as prim

class DerivedConstraint():
    
    def __init__(self, gcon_list, i_idx, j_idx):
        
        self.gcon_list = gcon_list
        self.i_idx = i_idx
        self.j_idx = j_idx
        
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

class GconPerp1(DerivedConstraint):
    
    """The perpendicular 1 constraint is discussed on slide 32 of lecture 8"""
    
    def __init__(self, body_i, ai_bar, bi_bar, body_j, cj_bar):
        
        DP1_1 = prim.GconDP1(body_i, ai_bar, body_j, cj_bar, constraint_func = None)
        DP1_2 = prim.GconDP1(body_i, bi_bar, body_j, cj_bar, constraint_func = None)
        
        super().__init__([DP1_1, DP1_2], body_i.idx, body_j.idx)
        
        self.DOF_constrained = 2 #constrains 2 DOF
        
    
class GconPerp2(DerivedConstraint):
    
    """The perpendicular 2 constraint is discussed on slide 34 of lecture 8"""
    
    def __init__(self, body_i, ai_bar, bi_bar, sp_i_bar, body_j, sq_j_bar):
        
        DP2_1 = prim.GconDP2(body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None)
        DP2_2 = prim.GconDP2(body_i, bi_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None)
        
        super().__init__([DP2_1, DP2_2], body_i.idx, body_j.idx)
        
        self.DOF_constrained = 2 #constrains 2 DOF
    
    
    