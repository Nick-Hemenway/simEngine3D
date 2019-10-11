import sys
sys.path.append('..')

import numpy as np
from Gcons import GconPrimitives as prim
from Gcons import GconIntermediate as inter

class Joint():
    
    def __init__(self, gcon_list):
        
        self.gcon_list = gcon_list
        
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
    
class JointSpherical(Joint):
    
    def __init__(self, body_i, sp_i_bar, body_j, sq_j_bar):
        
        CD_1 = prim.GconCD(body_i, sp_i_bar, body_j, sq_j_bar, c_vec = [1,0,0]) 
        CD_2 = prim.GconCD(body_i, sp_i_bar, body_j, sq_j_bar, c_vec = [0,1,0]) 
        CD_3 = prim.GconCD(body_i, sp_i_bar, body_j, sq_j_bar, c_vec = [0,0,1]) 
        
        super().__init__([CD_1, CD_2, CD_3])

class JointUniversal(Joint):
    
    def __init__(self, body_i, sp_i_bar, ai_bar, body_j, sq_j_bar, aj_bar):
        
        SJ = JointSpherical(body_i, sp_i_bar, body_j, sq_j_bar)
        DP1 = prim.GconDP1(body_i, ai_bar, body_j, aj_bar)
        
        super().__init__([SJ, DP1])
        
        
class JointCylindrical(Joint):
    
    def __init__(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar):
        
        perp1 = inter.GconPerp1(body_i, ai_bar, bi_bar, body_j, cj_bar)
        perp2 = inter.GconPerp2(body_i, ai_bar, bi_bar, sp_i_bar, body_j, sq_j_bar)

        super().__init__([perp1, perp2])  
    
        
class JointRevolute(Joint):
    
    def __init__(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar):
        
        SJ = JointSpherical(body_i, sp_i_bar, body_j, sq_j_bar)
        perp1 = inter.GconPerp1(body_i, ai_bar, bi_bar, body_j, cj_bar)
        
        super().__init__([SJ, perp1])

class JointTranslational(Joint):
    
    def __init__(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, aj_bar, cj_bar):
        
        CJ = JointCylindrical(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar)
        DP1 = prim.GconDP1(body_i, ai_bar, body_j, aj_bar)
        
        super().__init__([CJ, DP1])
        
        
        
        
        
        
        
        
        
        
        
        
        