#add folders to the path so python knows where to import self made modules from
import sys
import pathlib as pl
Gcon_folder = pl.Path('./Gcons/')
sys.path.append(str(Gcon_folder))

import numpy as np
import rigidbody
from Gcons import GconPrimitives, GconIntermediate, GconDerived
from utility import column, tilde

class System():
    
    num_bodies = 0
    
    #note that there is a distinction between number of joints and number of 
    #constraints because some joints provide more than one constraint
    num_constraints = 0 #total number of constraints in system
    num_joints = 0 #total number of joints in system
    
    def __init__(self, h = 0.01):
        
        #always create a ground/global reference frame as the first body when
        #instantiating a system
        self.global_frame = rigidbody.RigidBody(m = 0, J = [0,0,0], idx = 0, name = 'global_frame')
        self.gravity = self.global_frame.vector([0,0,-9.80665])
        
        self.bodies = []
        self.constraints = []
        self.t = 0 #global time is initialized to zero
        self.h = h
        
        self.history = {'r':[], 'p':[]}
        
    def set_system_time(self, t):
        
        self.t = t
        
    
    def add_body(self, m, J, r = None, p = None, r_dot = None, p_dot = None, r_ddot = None, p_ddot = None, name = None):
        
        #increment the  number of bodies in the system
        System.num_bodies += 1
        
        new_body = rigidbody.RigidBody(m = m, J = J, r = r, p = p, r_dot = r_dot, p_dot = p_dot, r_ddot = r_ddot, p_ddot = p_ddot, idx = System.num_bodies, name = name)
        new_body.add_force( force_vec = m*self.gravity, location = [0,0,0] ) #add gravity vector to body
        
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
    
    def M(self):
        
        M_mat = np.zeros( (3*self.num_bodies, 3*self.num_bodies) )
        
        for i, body in enumerate(self.bodies):
            
            idx = i*3
            M_mat[idx:idx+3, idx:idx+3] = body.m*np.eye(3)
            
        return M_mat
    
    def J(self):
        
        J_mat = np.zeros( (4*self.num_bodies, 4*self.num_bodies) )
        
        for i, body in enumerate(self.bodies):
            
            idx = i*4
            G = body.G
            J_mat[idx:idx+4, idx:idx+4] = 4 * G.T @ body.J @ G
            
        return J_mat
    
    def P(self):
        
        P_mat = np.zeros( (self.num_bodies, 4*self.num_bodies) )
        
        for i, body in enumerate(self.bodies):
            
            idx = i*4
            P_mat[i, idx:idx+4] = body.p.T
            
        return P_mat
    
    
    def partial_r(self):
        
        phi_r = np.zeros((System.num_constraints, 3*System.num_bodies))
        
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
            dphi_dri, dphi_drj = constraint.partial_r()
            
            #place partial derivatives in Jacobian but only if they are not ground
            if i_idx < 0: #body_i is ground
                                
                #assign position derivatives
                j = 3*j_idx
                phi_r[row:row + N, j:j+3] = dphi_drj
                
            elif j_idx < 0: #body_j is ground
                
                #assign position derivatives
                i = 3*i_idx
                phi_r[row:row + N, i:i+3] = dphi_dri
                
            else: #neither body is ground
                
                #assign position derivatives
                i, j = 3*i_idx, 3*j_idx
                phi_r[row:row + N, i:i+3] = dphi_dri
                phi_r[row:row + N, j:j+3] = dphi_drj
                
            row += N #move down in the Jacobian matrix
            
        return phi_r
    
    def partial_p(self):
        
        phi_p = np.zeros((System.num_constraints, 4*System.num_bodies))
        
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
            
            #place partial derivatives in Jacobian but only if they are not ground
            if i_idx < 0: #body_i is ground
                                
                #assign orientation derivatives
                j = 4*j_idx
                phi_p[row:row + N, j:j+4] = dphi_dpj
                
            elif j_idx < 0: #body_j is ground
                
                #assign orientation derivatives
                i = 4*i_idx
                phi_p[row:row + N, i:i+4] = dphi_dpi
                
            else: #neither body is ground
                
                #assign orientation derivatives
                i, j = 4*i_idx, 4*j_idx
                phi_p[row:row + N, i:i+4] = dphi_dpi
                phi_p[row:row + N, j:j+4] = dphi_dpj
                
            row += N #move down in the Jacobian matrix
            
        return phi_p
            
    def jacobian(self):
        
        #total number of constraints is the number of system constraints, plus
        #the number of bodies because each body has an Euler parameter constraint
        J = np.zeros((System.num_constraints + System.num_bodies, 7*System.num_bodies))
        offset = 3*System.num_bodies
        
        J[0:self.num_constraints, 0:3*self.num_bodies] = self.partial_r()
        J[0:self.num_constraints, 3*self.num_bodies::] = self.partial_p()
            
        #Add Euler Parameter constraint partial derivatives to bottom
        row = self.num_constraints
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
        
        r_vec = np.zeros( (3*self.num_bodies, 1))
        p_vec = np.zeros( (4*self.num_bodies, 1))
        
        for i, body in enumerate(self.bodies):
            
            #position coordinates
            idx = i*3
            r_vec[idx:idx+3] = getattr(body, pos)
            
            #orientation coordinates
            idx = i*4
            p_vec[idx:idx+4] = getattr(body, orient)
            
        return r_vec, p_vec
            
    
    def _set_generalized_coords(self, r, p, level = 0):
        
        levels = [('set_position','set_orientation'), \
                  ('set_vel',     'set_ang_vel'),     \
                  ('set_accel',   'set_ang_accel')]
        
        pos, orient = levels[level]
        
        for i, body in enumerate(self.bodies):
            
            #position coordinates
            idx = i*3
            getattr(body, pos)(r[idx:idx+3])
            
            #orientation coordinates
            idx = i*4
            getattr(body, orient)(p[idx:idx+4])
            
            
    def solve_position(self, tol = 1e-9, update_jacobian_every = 1):
        
        delta_mag = 2*tol #initial mag is twice tolerance to ensure loop starts
        iter_count = 0 #iteration count from last time jacobian was updated
        
        #compute the initial jacobian and q_vector
        J = self.jacobian()
        
        N = self.num_bodies
        q_old = np.zeros((7*N, 1))
        
        r_old, p_old = self._get_generalized_coords(level = 0)
        q_old[0:3*N] = r_old
        q_old[3*N::] = p_old
        
        q_new = q_old
        
        while delta_mag > tol:
        
            if iter_count >= update_jacobian_every:
                
                iter_count = 0 #reset the iteration count
                J = self.jacobian()
                
            phi = self.phi()
            delta = np.linalg.solve(-J, phi)
    
            q_new = q_old + delta
            self._set_generalized_coords(q_new[0:3*N], q_new[3*N::], level = 0)
            
            q_old = q_new
            
            delta_mag = np.linalg.norm(delta)
            iter_count += 1
        
    def solve_velocity(self):

        J = self.jacobian()
        nu = self.vel_rhs()
        
        q_dot = np.linalg.solve(J, nu)
        
        N = self.num_bodies
        
        self._set_generalized_coords(q_dot[0:3*N], q_dot[3*N::], level = 1)
        

    def solve_acceleration(self):

        J = self.jacobian()
        gamma = self.accel_rhs()
        
        q_ddot = np.linalg.solve(J, gamma)
        
        N = self.num_bodies
        
        self._set_generalized_coords(q_ddot[0:3*N], q_ddot[3*N::], level = 2)
        
    def solve_inverse_dynamics(self):
        
        #returns array lagrange multipliers from inverse solution
        
        #create LHS matrix
        lhs = np.zeros( (7*self.num_bodies, self.num_bodies + self.num_constraints) )
        
        offset1 = 3*self.num_bodies
        offset2 = self.num_constraints
        
        lhs[0:offset1, 0:offset2] = self.partial_r().T
        lhs[offset1::, 0:offset2] = self.partial_p().T
        lhs[offset1::, offset2::] = self.P().T
        
        #create RHS vector
        rhs = np.zeros( (7*self.num_bodies, 1) )
        
        offset = 3*self.num_bodies
        
        #preallocate generalized force and torque vectors
        F = np.zeros( (3*self.num_bodies, 1) )
        tau_hat = np.zeros((4*self.num_bodies, 1))

        #loop over each body and calculate net force and torque per body        
        for i, body in enumerate(self.bodies):
            
            #set net body force and torque to zero
            F_body = np.zeros((3,1))
            N_body = np.zeros((3,1))
            
            #calculate net force on body
            for force, loc in body.forces:
                
                fi = force.to_global()
                F_body += fi
                
            #determine how much torque each force creates about the CG
            #note that force is in global frame but torque vectors are in local frame
                N_body += tilde(loc) @ body.A.T @ fi
                
            for torque in body.torques:
                
                ti = torque.to_global()
                N_body += body.A.T @ ti #add torque in local reference frame to net torque
                
            #we must compute the generalized torques from the actual net torque
            G = body.G
            G_dot = body.G_dot
            tau_body = 2*G.T @ N_body + 8*G_dot.T @ body.J @ G_dot @ body.p
            
            #append body specific term to overall vector
            idx = i*3
            F[idx:idx+3] = F_body
            
            idx = i*4
            tau_hat[idx:idx+4] = tau_body
            
        #compute rhs
        r_ddot, p_ddot = self._get_generalized_coords(level = 2)
        rhs[0:offset] = F - self.M() @ r_ddot
        rhs[offset::] = tau_hat - self.J() @ p_ddot
        
        lagrange = np.linalg.solve(lhs, rhs)
        
        return lagrange
    
    def solve_kinematics(self, tol = 1e-9, update_jacobian_every = 1):
        
        self.solve_position(tol, update_jacobian_every)
        self.solve_velocity()
        self.solve_acceleration()
                
        
def main():
    
    pass


if __name__ == '__main__': main()