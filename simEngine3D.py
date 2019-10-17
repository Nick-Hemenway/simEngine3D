#add folders to the path so python knows where to import self made modules from
import sys
import pathlib as pl
Gcon_folder = pl.Path('./Gcons/')
sys.path.append(str(Gcon_folder))

import numpy as np
import rigidbody
from Gcons import GconPrimitives, GconIntermediate, GconDerived
from utility import column

class System():
    
    num_bodies = 0
    
    #note that there is a distinction between number of joints and number of 
    #constraints because some joints provide more than one constraint
    num_constraints = 0 #total number of constraints in system
    num_joints = 0 #total number of joints in system
    
    def __init__(self):
        
        #always create a ground/global reference frame as the first body when
        #instantiating a system
        self.global_frame = rigidbody.RigidBody(idx = 0, name = 'global_frame')
        
        self.bodies = []
        self.constraints = []
        self.t = 0 #global time is initialized to zero
        
    def set_system_time(self, t):
        
        self.t = t
    
    def add_body(self, r = None, p = None, r_dot = None, p_dot = None, r_ddot = None, p_ddot = None, name = None):
        
        #increment the  number of bodies in the system
        System.num_bodies += 1
        
        new_body = rigidbody.RigidBody(r = r, p = p, r_dot = r_dot, p_dot = p_dot, r_ddot = r_ddot, p_ddot = p_ddot, idx = System.num_bodies, name = name)
        self.bodies.append(new_body)
        
        return self.bodies[new_body.idx - 1]
    
    
    ######################   CONSTRAINTS   ######################
        
    def constraint_DP1(self, body_i, ai_bar, body_j, aj_bar, constraint_func = None):
        
        System.num_joints += 1
        
        con = GconPrimitives.GconDP1(body_i, ai_bar, body_j, aj_bar, constraint_func)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
        
    def constraint_DP2(self, body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None):
        
        System.num_joints += 1
        
        con = GconPrimitives.GconDP2(body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def constraint_D(self, body_i, sp_i_bar, body_j, sq_j_bar, constraint_func = None):
        
        System.num_joints += 1
        
        con = GconPrimitives.GconD(body_i, sp_i_bar, body_j, sq_j_bar, constraint_func)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def constraint_CD(self, body_i, sp_i_bar, body_j, sq_j_bar, c_vec, constraint_func = None):
        
        System.num_joints += 1
        
        con = GconPrimitives.GconCD(body_i, sp_i_bar, body_j, sq_j_bar, c_vec, constraint_func)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def constraint_Perp1(self, body_i, ai_bar, bi_bar, body_j, cj_bar):
        
        System.num_joints += 1
        
        con = GconIntermediate.GconPerp1(body_i, ai_bar, bi_bar, body_j, cj_bar)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def constraint_Perp2(self, body_i, ai_bar, bi_bar, sp_i_bar, body_j, sq_j_bar):
        
        System.num_joints += 1
        
        con = GconIntermediate.GconPerp2(body_i, ai_bar, bi_bar, sp_i_bar, body_j, sq_j_bar)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def constraint_func(self, f = None, f_prime = None, f_pprime = None):
        
        return GconPrimitives.ConstraintFunc(f, f_prime, f_pprime)
    
    
    ######################   DERIVED JOINTS   ######################
    
    def joint_spherical(self, body_i, sp_i_bar, body_j, sq_j_bar):
        
        System.num_joints += 1
        
        con = GconDerived.JointSpherical(body_i, sp_i_bar, body_j, sq_j_bar)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def joint_universal(self, body_i, sp_i_bar, ai_bar, body_j, sq_j_bar, aj_bar):
        
        System.num_joints += 1
        
        con = GconDerived.JointUniversal(body_i, sp_i_bar, ai_bar, body_j, sq_j_bar, aj_bar)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def joint_cylindrical(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar):
        
        System.num_joints += 1
        
        con = GconDerived.JointCylindrical(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def joint_revolute(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar):
        
        System.num_joints += 1
        
        con = GconDerived.JointRevolute(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def joint_translational(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, aj_bar, cj_bar):
        
        System.num_joints += 1
        
        con = GconDerived.JointTranslational(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, aj_bar, cj_bar)
        System.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[System.num_joints - 1]
    
    def phi(self):
        
        #preallocate an array for phi
        phi_vec = np.zeros((System.num_constraints + System.num_bodies, 1))
        
        old_idx = 0
        #add all algebraic constraint equations to phi vector
        for constraint in self.constraints:
            
            new_idx = old_idx + constraint.DOF_constrained
            phi_vec[old_idx:new_idx] = constraint.val(self.t)
            
            old_idx = new_idx
        
        #append Euler parameter normalization constraints to end of vector
        for body in self.bodies:
            
            p = body.p
            phi_vec[new_idx] = p.T @ p - 1.0
            new_idx += 1
            
        return phi_vec
    
    def jacobian(self):
        
        #total number of constraints is the number of system constraints, plus
        #the number of bodies because each body has an Euler parameter constraint
        J = np.zeros((System.num_constraints + System.num_bodies, 7*System.num_bodies))
        offset = 3*System.num_bodies
        
        row = 0
        
        for constraint in self.constraints:
            
            N = constraint.DOF_constrained #number of rows in Jacobian
            #bodies start indexing at 1 because the global reference frame is at
            #index zero. However, the global reference frame doesn't end up in 
            #the jacobian so we must subtract one from the index for each body
            #because the array starts indexing from zero.
            i_idx = constraint.i_idx - 1 
            j_idx = constraint.j_idx - 1
            
            
            #compute partial derivatives
            dphi_dpi, dphi_dpj = constraint.partial_p()
            dphi_dri, dphi_drj = constraint.partial_r()
            
            #place partial derivatives in Jacobian but only if they are not ground
            if i_idx < 0: #body_i is ground
                                
                #assign position derivatives
                j = 3*j_idx
                J[row:row + N, j:j+3] = dphi_drj
                
                #assign orientation derivatives
                j = 4*j_idx + offset
                J[row:row + N, j:j+4] = dphi_dpj
                
            elif j_idx < 0: #body_j is ground
                
                #assign position derivatives
                i = 3*i_idx
                J[row:row + N, i:i+3] = dphi_dri
                
                #assign orientation derivatives
                i = 4*i_idx + offset
                J[row:row + N, i:i+4] = dphi_dpi   
                
            else: #neither body is ground
                
                #assign position derivatives
                i, j = 3*i_idx, 3*j_idx
                J[row:row + N, i:i+3] = dphi_dri
                J[row:row + N, j:j+3] = dphi_drj
                
                #assign orientation derivatives
                i = 4*i_idx + offset 
                j = 4*j_idx + offset
                J[row:row + N, i:i+4] = dphi_dpi
                J[row:row + N, j:j+4] = dphi_dpj
            
            row += N #move down in the Jacobian matrix
            
        #Add Euler Parameter constraint partial derivatives to bottom
        for i, body in enumerate(self.bodies):
            
            idx = i*4 + offset
            J[row, idx:idx+4] = 2*body.p.flatten()
            row += 1
            
        return J
                
    def vel_rhs(self):
        
        nu = np.zeros( (System.num_constraints + System.num_bodies, 1) )
        
        old_idx = 0
        #add all algebraic constraint equations to phi vector
        for constraint in self.constraints:
            
            new_idx = old_idx + constraint.DOF_constrained
            nu[old_idx:new_idx] = constraint.vel_rhs(self.t)
            
            old_idx = new_idx
        
        #nu for euler parameters equals zero and we have already preallocate
        #the array full of zeros so there's no need to worry about the bottom
        #of the array
        
        return column(nu)
        
    def accel_rhs(self):
    
        gamma = np.zeros( (System.num_constraints + System.num_bodies, 1) )
        
        old_idx = 0
        #add all algebraic constraint equations to phi vector
        for constraint in self.constraints:
            
            new_idx = old_idx + constraint.DOF_constrained
            gamma[old_idx:new_idx] = constraint.accel_rhs(self.t)
            
            old_idx = new_idx
        
        #append Euler parameter normalization constraints to end of vector
        for body in self.bodies:
            
            p_dot = body.p_dot
            gamma[new_idx] = -2*p_dot.T @ p_dot
            new_idx += 1
            
        return column(gamma)
    
    def _get_generalized_coords(self, level =  0):
        
        levels = [('r','p'), ('r_dot', 'p_dot'), ('r_ddot', 'p_ddot')]
        pos, orient = levels[level]
        
        q_vec = np.zeros( (7*self.num_bodies, 1))
        offset = 3*self.num_bodies
        
        for i, body in enumerate(self.bodies):
            
            #position coordinates
            idx = i*3
            q_vec[idx:idx+3] = getattr(body, pos)
            
            #orientation coordinates
            idx = i*4 + offset
            q_vec[idx:idx+4] = getattr(body, orient)
            
        return q_vec
            
    
    def _set_generalized_coords(self, q, level = 0):
        
        levels = [('set_position','set_orientation'), \
                  ('set_vel',     'set_ang_vel'),     \
                  ('set_accel',   'set_ang_accel')]
        
        pos, orient = levels[level]
        
        offset = 3*self.num_bodies
        
        for i, body in enumerate(self.bodies):
            
            #position coordinates
            idx = i*3
            getattr(body, pos)(q[idx:idx+3])
            
            #orientation coordinates
            idx = i*4 + offset
            getattr(body, orient)(q[idx:idx+4])
            
            
    def solve_position(self, tol = 1e-9):
        
        delta_mag = 2*tol
        
        while delta_mag > tol:
        
            J = self.jacobian()
            phi = self.phi()
            
            delta = np.linalg.solve(-J, phi)

            q_new = self._get_generalized_coords(level = 0) + delta
            self._set_generalized_coords(q_new, level = 0)
            
            delta_mag = np.linalg.norm(delta)
        
    def solve_velocity(self):

        J = self.jacobian()
        nu = self.vel_rhs()
        
        q_dot = np.linalg.solve(J, nu)
        
        self._set_generalized_coords(q_dot, level = 1)
        

    def solve_acceleration(self):

        J = self.jacobian()
        gamma = self.accel_rhs()
        
        q_ddot = np.linalg.solve(J, gamma)
        
        self._set_generalized_coords(q_ddot, level = 2)
        
    
def main():
    
    pass


if __name__ == '__main__': main()