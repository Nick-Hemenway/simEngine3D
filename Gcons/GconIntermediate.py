#add the parent folder to the search path for importing modules
import sys
sys.path.append('..')

import numpy as np
from utility import column
import GconPrimitives as prim

class Perp1():
    
    """The perpendicular 1 constraint is discussed on slide 32 of lecture 8"""
    
    def __init__(self, body_i, ai_bar, bi_bar, body_j, cj_bar):
        
        self.dp1_1 = prim.GconDP1(body_i, ai_bar, body_j, cj_bar, constraint_func = None)
        self.dp1_2 = prim.GconDP1(body_i, bi_bar, body_j, cj_bar, constraint_func = None)
        
    def val(self, t):
        
        vals = [self.dp1_1.val(t), self.dp1_2.val(t)]
        return column(vals)
    
    def vel_rhs(self, t):
        
        nu = [self.dp1_1.vel_rhs(t), self.dp1_2.vel_rhs(t)]
        return column(nu)
    
    def accel_rhs(self, t):
        
        gamma = [self.dp1_1.accel_rhs(t), self.dp1_2.accel_rhs(t)]
        return column(gamma)
    
    def partial_p(self):
        
        dphi_dpi_1, dphi_dpj_1 = self.dp1_1.partial_p()
        dphi_dpi_2, dphi_dpj_2 = self.dp1_2.partial_p()
        
        dphi_dpi = np.vstack( (dphi_dpi_1, dphi_dpi_2) )
        dphi_dpj = np.vstack( (dphi_dpj_1, dphi_dpj_2) )
        
        return dphi_dpi, dphi_dpj
    
    def partial_r(self):
        
        dphi_dri_1, dphi_drj_1 = self.dp1_1.partial_r()
        dphi_dri_2, dphi_drj_2 = self.dp1_2.partial_r()
        
        dphi_dri = np.vstack( (dphi_dri_1, dphi_dri_2) )
        dphi_drj = np.vstack( (dphi_drj_1, dphi_drj_2) )
        
        return dphi_dri, dphi_drj
    
class Perp2():
    
    """The perpendicular 2 constraint is discussed on slide 34 of lecture 8"""
    
    def __init__(self, body_i, ai_bar, bi_bar, sp_i_bar, body_j, sq_j_bar):
        
        self.dp2_1 = prim.GconDP2(body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None)
        self.dp2_2 = prim.GconDP2(body_i, bi_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None)
        
    def val(self, t):
        
        vals = [self.dp2_1.val(t), self.dp2_2.val(t)]
        return column(vals)
    
    def vel_rhs(self, t):
        
        nu = [self.dp2_1.vel_rhs(t), self.dp2_2.vel_rhs(t)]
        return column(nu)
    
    def accel_rhs(self, t):
        
        gamma = [self.dp2_1.accel_rhs(t), self.dp2_2.accel_rhs(t)]
        return column(gamma)
    
    def partial_p(self):
        
        dphi_dpi_1, dphi_dpj_1 = self.dp2_1.partial_p()
        dphi_dpi_2, dphi_dpj_2 = self.dp2_2.partial_p()
        
        dphi_dpi = np.vstack( (dphi_dpi_1, dphi_dpi_2) )
        dphi_dpj = np.vstack( (dphi_dpj_1, dphi_dpj_2) )
        
        return dphi_dpi, dphi_dpj
    
    def partial_r(self):
        
        dphi_dri_1, dphi_drj_1 = self.dp2_1.partial_r()
        dphi_dri_2, dphi_drj_2 = self.dp2_2.partial_r()
        
        dphi_dri = np.vstack( (dphi_dri_1, dphi_dri_2) )
        dphi_drj = np.vstack( (dphi_drj_1, dphi_drj_2) )
        
        return dphi_dri, dphi_drj    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    