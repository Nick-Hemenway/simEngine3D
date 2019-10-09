#add folders to the path so python knows where to import self made modules from
import sys
sys.path.append('Gcons')

from Gcons import GconPrimitives
import rigidbody as body

class System():
    
    num_bodies = 0
    
    #note that there is a distinction between number of joints and number of 
    #constraints because some joints provide more than one constraint
    num_constraints = 0 #total number of constraints system
    num_joints = 0 #total number of joints in system
    
    def __init__(self):
        
        #always create a ground/global reference frame as the first body when
        #instantiating a system
        self.global_frame = body.RigidBody(idx = 0, name = 'global_frame')
        
        self.bodies = []
        self.constraints = []
        self.t = 0 #global time is initialized to zero
    
    def add_body(self, r = None, p = None, r_dot = None, p_dot = None, name = None):
        
        #increment the  number of bodies in the system
        System.num_bodies += 1
        
        new_body = body.RigidBody(r = r, p = p, r_dot = r_dot, p_dot = p_dot, idx = System.num_bodies, name = name)
        self.bodies.append(new_body)
        
        return self.bodies[new_body.idx - 1]
        
    def constraint_DP1(self, body_i, ai_bar, body_j, aj_bar, constraint_func = None):
        
        System.num_constraints += 1
        System.num_joints += 1
        
        con = GconPrimitives.GconDP1(body_i, ai_bar, body_j, aj_bar, constraint_func)
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
        
    def constraint_DP2(self, body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None):
        
        System.num_constraints += 1
        System.num_joints += 1
        
        con = GconPrimitives.GconDP2(body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func)
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def constraint_D(self, body_i, sp_i_bar, body_j, sq_j_bar, constraint_func = None):
        
        System.num_constraints += 1
        System.num_joints += 1
        
        con = GconPrimitives.GconD(body_i, sp_i_bar, body_j, sq_j_bar, constraint_func)
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def constraint_CD(self, body_i, sp_i_bar, body_j, sq_j_bar, c_vec, constraint_func = None):
        
        System.num_constraints += 1
        System.num_joints += 1
        
        con = GconPrimitives.GconCD(body_i, sp_i_bar, body_j, sq_j_bar, c_vec, constraint_func)
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def constraint_func(self, f = None, f_prime = None, f_pprime = None):
        
        return GconPrimitives.ConstraintFunc(f, f_prime, f_pprime)
    
    
    def jacobian(self):
        
        pass
    
    

def main():
    
    pass


if __name__ == '__main__': main()