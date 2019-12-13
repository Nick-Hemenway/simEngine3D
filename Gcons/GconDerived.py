import sys
sys.path.append('..')

from Gcons import GconPrimitives as prim
from Gcons import GconIntermediate as inter
from Gcons.GconIntermediate import DerivedConstraint

class JointSpherical(DerivedConstraint):
    
    def __init__(self, body_i, sp_i_bar, body_j, sq_j_bar):
        
        CD_1 = prim.GconCD(body_i, sp_i_bar, body_j, sq_j_bar, c_vec = [1,0,0]) 
        CD_2 = prim.GconCD(body_i, sp_i_bar, body_j, sq_j_bar, c_vec = [0,1,0]) 
        CD_3 = prim.GconCD(body_i, sp_i_bar, body_j, sq_j_bar, c_vec = [0,0,1]) 
        self.DOF_constrained = 3 #constrains 3 DOF
        
        super().__init__([CD_1, CD_2, CD_3], body_i, body_j, sp_i_bar = sp_i_bar, sq_j_bar = sq_j_bar)

class JointUniversal(DerivedConstraint):
    
    def __init__(self, body_i, sp_i_bar, ai_bar, body_j, sq_j_bar, aj_bar):
        
        SJ = JointSpherical(body_i, sp_i_bar, body_j, sq_j_bar)
        DP1 = prim.GconDP1(body_i, ai_bar, body_j, aj_bar)
        self.DOF_constrained = 4 #constrains 4 DOF
        
        super().__init__([SJ, DP1], body_i, body_j, sp_i_bar = sp_i_bar, sq_j_bar = sq_j_bar)
        
        
class JointCylindrical(DerivedConstraint):
    
    def __init__(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar):
        
        perp1 = inter.GconPerp1(body_i, ai_bar, bi_bar, body_j, cj_bar)
        perp2 = inter.GconPerp2(body_i, ai_bar, bi_bar, sp_i_bar, body_j, sq_j_bar)
        self.DOF_constrained = 4 #constrains 4 DOF

        super().__init__([perp1, perp2], body_i, body_j, sp_i_bar = sp_i_bar, sq_j_bar = sq_j_bar)  
    
        
class JointRevolute(DerivedConstraint):
    
    def __init__(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar):
        
        SJ = JointSpherical(body_i, sp_i_bar, body_j, sq_j_bar)
        perp1 = inter.GconPerp1(body_i, ai_bar, bi_bar, body_j, cj_bar)
        self.DOF_constrained = 5 #constrains 5 DOF
        
        super().__init__([SJ, perp1], body_i, body_j, sp_i_bar = sp_i_bar, sq_j_bar = sq_j_bar)

class JointTranslational(DerivedConstraint):
    
    def __init__(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, aj_bar, cj_bar):
        
        CJ = JointCylindrical(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar)
        DP1 = prim.GconDP1(body_i, ai_bar, body_j, aj_bar)
        self.DOF_constrained = 5 #constrains 5 DOF
        
        super().__init__([CJ, DP1], body_i, body_j, sp_i_bar = sp_i_bar, sq_j_bar = sq_j_bar)
        
        
        
        
        
        
        
        
        
        
        
        
        