import numpy as np
from simengine3D import rigidbody
from simengine3D import solvers
from simengine3D import loads
import simengine3D.geometric_constraints as gcon
from simengine3D.utility import column

class System():
    
    def __init__(self, h = 0.01, order = 2, solver = solvers.BDF_Solver):
    
        self.num_bodies = 0
    
        #note that there is a distinction between number of joints and number of 
        #constraints because some joints provide more than one constraint
        self.num_constraints = 0 #total number of constraints in system
        self.num_joints = 0 #total number of joints in system
        
        #always create a ground/global reference frame as the first body when
        #instantiating a system
        self.global_frame = rigidbody.RigidBody(m = 0, J = [0,0,0], idx = 0, name = 'global_frame')
        self.gravity = self.global_frame.vector([0,0,-9.80665])
        
        self.bodies = []
        self.constraints = []
        self.t = 0 #global time is initialized to zero
        
        #solver settings
        self.order = order #order 
        self.h = h #step size
        
        self.step_num = 0 #num of steps in simulation
        
        #a dictionary filled with lists to track the history of the system. Note that
        #oldest entries come first in the list
        self.history = {'r':[], 'p':[], 'r_dot':[], 'p_dot':[], 'r_ddot':[], 'p_ddot':[]}
        
        self.solver = solver(h = h, order = order)
        self.is_initialized = False
        
        #lists to store system level forces and torques (e.g. TSDA, RSDA, etc.)
        self.sys_forces = []
        self.sys_torques = []
        
    def print_time(self, num_steps = 20):
        
        if self.step_num % num_steps == 0:
            print(f't: {self.t:.3f}')
        
    def set_gravity(self, vec):
        
        #vec is an array-like
        self.gravity = self.global_frame.vector(vec)
        
    def set_solver_order(self, order):
        
        self.order = order
        self.solver.set_order(order)
        
    def set_step_size(self, h):
        
        self.h = h
        self.solver.h = h
        
    def set_system_time(self, t):
        
        self.t = t
        
    def history_set(self):
        
        """This function appends the current state of the system to the history
        list attribute"""
        
        r, p = self._get_generalized_coords(level = 0)
        r_dot, p_dot = self._get_generalized_coords(level = 1)
        r_ddot, p_ddot = self._get_generalized_coords(level = 2)
        
        keys = ['r', 'p', 'r_dot', 'p_dot', 'r_ddot', 'p_ddot']
        attributes = [r, p, r_dot, p_dot, r_ddot, p_ddot]
        
        for key, attr in zip(keys, attributes):
            self.history[key].append(attr)
        
    def history_delete_oldest(self):
        
        for key in self.history.keys():
            
            self.history[key].pop(0)
        
    
    def add_body(self, m, J, r = None, p = None, r_dot = None, p_dot = None, r_ddot = None, p_ddot = None, name = None):
        
        #increment the  number of bodies in the system
        self.num_bodies += 1
        
        new_body = rigidbody.RigidBody(m = m, J = J, r = r, p = p, r_dot = r_dot, p_dot = p_dot, r_ddot = r_ddot, p_ddot = p_ddot, idx = self.num_bodies, name = name)
        Fg = loads.const_force( m*self.gravity, loc = [0,0,0] ) #returns force object for gravity
        new_body.add_force(force = Fg) #add gravity vector to body
        
        self.bodies.append(new_body)
        
        return self.bodies[new_body.idx - 1]
    
    
    ######################   CONSTRAINTS   ######################
        
    def constraint_DP1(self, body_i, ai_bar, body_j, aj_bar, constraint_func = None):
        
        self.num_joints += 1
        
        con = gcon.GconDP1(body_i, ai_bar, body_j, aj_bar, constraint_func)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
        
    def constraint_DP2(self, body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func = None):
        
        self.num_joints += 1
        
        con = gcon.GconDP2(body_i, ai_bar, sp_i_bar, body_j, sq_j_bar, constraint_func)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    def constraint_D(self, body_i, sp_i_bar, body_j, sq_j_bar, constraint_func = None):
        
        self.num_joints += 1
        
        con = gcon.GconD(body_i, sp_i_bar, body_j, sq_j_bar, constraint_func)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    def constraint_CD(self, body_i, sp_i_bar, body_j, sq_j_bar, c_vec, constraint_func = None):
        
        self.num_joints += 1
        
        con = gcon.GconCD(body_i, sp_i_bar, body_j, sq_j_bar, c_vec, constraint_func)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    def constraint_Perp1(self, body_i, ai_bar, bi_bar, body_j, cj_bar):
        
        self.num_joints += 1
        
        con = gcon.GconPerp1(body_i, ai_bar, bi_bar, body_j, cj_bar)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    def constraint_Perp2(self, body_i, ai_bar, bi_bar, sp_i_bar, body_j, sq_j_bar):
        
        self.num_joints += 1
        
        con = gcon.GconPerp2(body_i, ai_bar, bi_bar, sp_i_bar, body_j, sq_j_bar)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    def constraint_func(self, f = None, f_prime = None, f_pprime = None):
        
        return gcon.ConstraintFunc(f, f_prime, f_pprime)
    
    ######################   DERIVED JOINTS   ######################
    
    def joint_spherical(self, body_i, sp_i_bar, body_j, sq_j_bar):
        
        self.num_joints += 1
        
        con = gcon.JointSpherical(body_i, sp_i_bar, body_j, sq_j_bar)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    def joint_universal(self, body_i, sp_i_bar, ai_bar, body_j, sq_j_bar, aj_bar):
        
        self.num_joints += 1
        
        con = gcon.JointUniversal(body_i, sp_i_bar, ai_bar, body_j, sq_j_bar, aj_bar)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    def joint_cylindrical(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar):
        
        self.num_joints += 1
        
        con = gcon.JointCylindrical(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    def joint_revolute(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar):
        
        self.num_joints += 1
        
        con = gcon.JointRevolute(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    def joint_translational(self, body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, aj_bar, cj_bar):
        
        self.num_joints += 1
        
        con = gcon.JointTranslational(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, aj_bar, cj_bar)
        self.num_constraints += con.DOF_constrained
        
        self.constraints.append(con)
        
        return self.constraints[self.num_joints - 1]
    
    ######################   LOADS   ######################
    
    def force_constant(self, body, vec, location):
        
        """Creates a constant magnitude force object given the force vector and location
        
        body: rigid body object. Body to which the constant force should be applied to
        
        vec: a vector object created in the desired reference frame
        loc: an array-like specifying the point on the body at which the force is applied
        """
        
        F_obj = loads.const_force(vec, location)
        body.add_force(F_obj)
        
        
    def force_TSDA(self, body_i, sp_i_bar, body_j, sq_j_bar, k, l_0, c, h):
        
        tsda = loads.TSDA(body_i, sp_i_bar, body_j, sq_j_bar, k, l_0, c, h)
        self.sys_forces.append(tsda)
        
        Fi, Fj = tsda.create_force_objs()
        
        body_i.add_force(Fi)
        body_j.add_force(Fj)
        
        return self.sys_forces[-1]
        
    
    def torque_constant(self, body, vec):
        
        """Creates a constant magnitude torque object given the torque vector
        
        body: rigid body object. Body to which the constant torque should be applied to
        vec: a vector object created in the desired reference frame
        """
        
        T_obj = loads.const_torque(vec)
        body.add_torque(T_obj)
        

    def torque_RSDA(self, body_i, ai_bar, bi_bar, body_j, aj_bar, bj_bar, k, theta_0, c, h):
        
        rsda = loads.RSDA(body_i, ai_bar, bi_bar, body_j, aj_bar, bj_bar, k, theta_0, c, h)
        self.sys_torques.append(rsda)
        Ti, Tj = rsda.create_torque_objs()
        
        body_i.add_torque(Ti)
        body_j.add_torque(Tj)
        
        return self.sys_torques[-1]
        
    
    ######################   SYSTEM COMMANDS   ######################
    
    def phi(self):
        
        #preallocate an array for phi
        phi_vec = np.zeros((self.num_constraints, 1))
        
        old_idx = 0
        #add all algebraic constraint equations to phi vector
        for constraint in self.constraints:
            
            new_idx = old_idx + constraint.DOF_constrained
            phi_vec[old_idx:new_idx] = constraint.val(self.t)
            
            old_idx = new_idx
            
        return phi_vec
    
    def phi_euler(self):
        
        phi_vec_euler = np.zeros( (self.num_bodies, 1) )
        
        #append Euler parameter normalization constraints to end of vector
        for i, body in enumerate(self.bodies):
            
            p = body.p
            phi_vec_euler[i] = 0.5*p.T @ p - 0.5 #defined this way because of way
                                                 #we define Psi later
#            phi_vec_euler[i] = p.T @ p - 1.0

        return phi_vec_euler        
    
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
        
        #same as phi_r
        
        phi_r = np.zeros((self.num_constraints, 3*self.num_bodies))
        
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
        
        #same as phi_p
        phi_p = np.zeros((self.num_constraints, 4*self.num_bodies))
        
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
        J = np.zeros((self.num_constraints + self.num_bodies, 7*self.num_bodies))
        offset = 3*self.num_bodies
        
        J[0:self.num_constraints, 0:3*self.num_bodies] = self.partial_r()
        J[0:self.num_constraints, 3*self.num_bodies::] = self.partial_p()
            
        #Add Euler Parameter constraint partial derivatives to bottom
        row = self.num_constraints
        for i, body in enumerate(self.bodies):
            
            idx = i*4 + offset
#            J[row, idx:idx+4] = 2*body.p.flatten()
            J[row, idx:idx+4] = body.p.flatten() #defined this way because of defining
                                                 #phi^p as 0.5p^T p - 0.5 = 0 
            row += 1
            
        return J
    
    def nu(self):
        
        nu_arr = np.zeros( (self.num_constraints, 1) )
        
        old_idx = 0
        #add all algebraic constraint equations to phi vector
        for constraint in self.constraints:
            
            new_idx = old_idx + constraint.DOF_constrained
            nu_arr[old_idx:new_idx] = constraint.vel_rhs(self.t)
            
            old_idx = new_idx
        
        return column(nu_arr)    

    def nu_euler(self):
        
        return np.zeros( (self.num_bodies, 1) )
        
                
    def vel_rhs(self):
        
        rhs = np.zeros( (self.num_constraints + self.num_bodies, 1) )
        
        rhs[0:self.num_constraints] = self.nu() 
        
        #nu for euler parameters equals zero and we have already preallocate
        #the array full of zeros so there's no need to worry about the bottom
        #of the array
        
        return column(rhs)
    
    def gamma(self):
        
        gamma_arr = np.zeros( (self.num_constraints, 1) )
        
        old_idx = 0
        #add all algebraic constraint equations to phi vector
        for constraint in self.constraints:
            
            new_idx = old_idx + constraint.DOF_constrained
            gamma_arr[old_idx:new_idx] = constraint.accel_rhs(self.t)
            
            old_idx = new_idx
            
        return column(gamma_arr)
        
    def gamma_euler(self):
        
        gamma_arr_euler = np.zeros( (self.num_bodies, 1) )
        
        for i, body in enumerate(self.bodies):
            
            p_dot = body.p_dot
#            gamma_arr_euler[i] = -2*p_dot.T @ p_dot
            gamma_arr_euler[i] = -p_dot.T @ p_dot #defined this way because we defined
            #phi^p as 0.5p^T p - 0.5 = 0 

        return gamma_arr_euler
    
    def accel_rhs(self):
    
        gamma = np.zeros( (self.num_constraints + self.num_bodies, 1) )
        
        gamma[0:self.num_constraints] = self.gamma()
        gamma[self.num_constraints::] = self.gamma_euler()
            
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
            
            
    def _set_lagrange_params(self, lagrange, lagrange_euler):
        
        self.lagrange = lagrange
        self.lagrange_euler = lagrange_euler
        
        idx = 0
        
        for constraint in self.constraints:
            
            dof = constraint.DOF_constrained
            constraint.lagrange = column(lagrange[idx:idx+dof])
            idx += dof
            
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
        
        phi = np.zeros( (self.num_constraints + self.num_bodies, 1) )
        
        while delta_mag > tol:
        
            if iter_count >= update_jacobian_every:
                
                iter_count = 0 #reset the iteration count
                J = self.jacobian()
                
            phi[0:self.num_constraints] = self.phi()
            phi[self.num_constraints::] = self.phi_euler()
            
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
        
    def generalized_forces(self):
        
        F = np.zeros( (3*self.num_bodies, 1) )
        
        for i, body in enumerate(self.bodies):
            
            #set net body force and torque to zero
            F_body = np.zeros((3,1))
            
            for force in body.forces:
                
                fi = force.in_global(self.t)
                F_body += fi
                
                        #append body specific term to overall vector
            idx = i*3
            F[idx:idx+3] = F_body
            
        return F
    
    def generalized_torques(self):
        
        tau_hat = np.zeros((4*self.num_bodies, 1))
        
        #loop over each body and calculate net force and torque per body        
        for i, body in enumerate(self.bodies):
            
            #set net body force and torque to zero
            N_body = np.zeros((3,1))
            
            #calculate net force on body
            for force in body.forces:
                
                torque_from_force = force.resulting_torque(self.t, frame = body)
                N_body += torque_from_force.to_local(body)
                
            for torque in body.torques:
                
                ti = torque.to_local(self.t, body)
                N_body += ti #add torque in local reference frame to net torque
                
            #we must compute the generalized torques from the actual net torque
            G = body.G
            G_dot = body.G_dot
            tau_body = 2*G.T @ N_body + 8*G_dot.T @ body.J @ G_dot @ body.p
            
            #append body specific term to overall vector
            idx = i*4
            tau_hat[idx:idx+4] = tau_body
        
        return tau_hat
        
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
        
        F = self.generalized_forces()
        tau_hat = self.generalized_torques()
            
        #compute rhs
        r_ddot, p_ddot = self._get_generalized_coords(level = 2)
        rhs[0:offset] = F - self.M() @ r_ddot
        rhs[offset::] = tau_hat - self.J() @ p_ddot
        
        lam = np.linalg.solve(lhs, rhs)
        
        lagrange = lam[0:6*self.num_bodies]
        lagrange_euler = lam[6*self.num_bodies::]
        
        self._set_lagrange_params(lagrange, lagrange_euler)
        
    def solve_kinematics(self, tol = 1e-9, update_jacobian_every = 1):
        
        self.solve_position(tol, update_jacobian_every)
        self.solve_velocity()
        self.solve_acceleration()
        
    def initialize(self):
        
        """Initial conditions for position and velocity must be given. This 
        Function calculates the accelerations and lagrange multipliers for time 
        t = 0, sets the accelerations for each of the bodies, sets the lagrange 
        multipliers of the system for t = 0 and then takes a snapshot of the 
        history of the system at t = 0"""
        
        self.is_initialized = True
        
        #initialization comes from slide 16 of L15
        nb = self.num_bodies
        nc = self.num_constraints
        N = 8*nb + nc
        
        LHS = np.zeros( (N, N) )
        RHS = np.zeros( (N, 1) )
        
        #formulate LHS
        col_offset1 = 3*nb
        col_offset2 = col_offset1 + 4*nb
        col_offset3 = col_offset2 + nb
        
        row_offset1 = 3*nb
        row_offset2 = row_offset1 + 4*nb
        row_offset3 = row_offset2 + nb
        
        #row 1
        LHS[0:row_offset1, 0:col_offset1] = self.M()
        LHS[0:row_offset1, col_offset3::]  = self.partial_r().T
        
        #row 2
        LHS[row_offset1 : row_offset2, col_offset1 : col_offset2]  = self.J()
        LHS[row_offset1 : row_offset2, col_offset2 : col_offset3]  = self.P().T
        LHS[row_offset1 : row_offset2, col_offset3::]  = self.partial_p().T

        #row 3        
        LHS[row_offset2 : row_offset3, col_offset1 : col_offset2]  = self.P()

        #row 4
        LHS[row_offset3:: , 0 : col_offset1]  = self.partial_r()
        LHS[row_offset3:: , col_offset1 : col_offset2]  = self.partial_p()
        
        #formulate RHS
        RHS[0:row_offset1] = self.generalized_forces()
        RHS[row_offset1 : row_offset2] = self.generalized_torques()
        RHS[row_offset2 : row_offset3] = self.gamma_euler()
        RHS[row_offset3::] = self.gamma()
        
        sol_0 = np.linalg.solve(LHS, RHS)        
        
        r_ddot_0 = column(sol_0[0:row_offset1])
        p_ddot_0 = column(sol_0[row_offset1:row_offset2])
        
        lagrange_euler = column(sol_0[row_offset2:row_offset3])
        lagrange = column(sol_0[row_offset3::])
        
        self._set_generalized_coords(r_ddot_0, p_ddot_0, level = 2)

        self._set_lagrange_params(lagrange, lagrange_euler)
        
        self.history_set()
        
    def residual(self):
        
        """residual is shown on slide 5 of L15"""
        
        beta_0 = self.solver.beta_0
        h = self.h
        
        nb = self.num_bodies
        nc = self.num_constraints
        
        offset1 = 3*nb
        offset2 = offset1 + 4*nb
        offset3 = offset2 + nb
        
        r_ddot, p_ddot = self._get_generalized_coords(level = 2)
        
        g = np.zeros( (8*nb + nc, 1) )
        
        g[0:offset1] = self.M() @ r_ddot + self.partial_r().T @ self.lagrange - self.generalized_forces()
        g[offset1:offset2] = self.J() @ p_ddot + self.partial_p().T @ self.lagrange + self.P().T @ self.lagrange_euler - self.generalized_torques()
        
        coeff = 1/(beta_0**2 * h**2)
        
        g[offset2:offset3] = coeff* self.phi_euler()
        g[offset3::] = coeff* self.phi()
        
        return g
    
    def sensitivity(self):
        
        nc = self.num_constraints
        nb = self.num_bodies
        N = 8*nb + nc
        
        psi = np.zeros( (N, N) )
        
        #offsets
        col_offset1 = 3*nb
        col_offset2 = col_offset1 + 4*nb
        col_offset3 = col_offset2 + nb
        
        row_offset1 = 3*nb
        row_offset2 = row_offset1 + 4*nb
        row_offset3 = row_offset2 + nb
        
        #row 1
        psi[0:row_offset1, col_offset3::]  = self.partial_r().T
        
        psi[0:row_offset1, 0:col_offset1]  = self.M() # ~psi_11
        psi[0:row_offset1, col_offset1:col_offset2]  = 0 # ~psi_12
        
        #row 2
        psi[row_offset1 : row_offset2, col_offset2 : col_offset3]  = self.P().T
        psi[row_offset1 : row_offset2, col_offset3::]  = self.partial_p().T

        psi[row_offset1 : row_offset2, 0:col_offset1]  = 0 # ~psi_21
        psi[row_offset1 : row_offset2, col_offset1:col_offset2]  = self.J() # ~psi_22
        
        #row 3        
        psi[row_offset2 : row_offset3, col_offset1 : col_offset2]  = self.P()

        #row 4
        psi[row_offset3:: , 0 : col_offset1]  = self.partial_r()
        psi[row_offset3:: , col_offset1 : col_offset2]  = self.partial_p()
        
        return psi
        
    
    def _get_history_array(self, key, order):
        
        param = self.history[key][-order::] #list of numpy arrays
        arr = np.hstack(param)
        
        return arr                
    
    def step(self, tol = 1e-3):
       
        nb = self.num_bodies
        nc = self.num_constraints
        delta_mag = 2*tol
        
        #STEP 0: PRIME NEW STEP
        #initialize the problem if it has not been already
        if not self.is_initialized:
            self.initialize()
            
        #step time forward one step and increment the step counter
        self.t += self.h
        self.step_num += 1
        
        #set order of solver based on how much history is present
        history_cnt = len(self.history['r'])
        order = min(self.order, history_cnt)
        self.solver.set_order(order)
        
        #grab history 
        r_h = self._get_history_array(key = 'r', order = order)
        p_h = self._get_history_array(key = 'p', order = order)
        
        r_dot_h = self._get_history_array(key = 'r_dot', order = order)
        p_dot_h = self._get_history_array(key = 'p_dot', order = order)
        
        r_ddot, p_ddot = self._get_generalized_coords(level = 2)

        sol = np.zeros( (8*nb + nc, 1) )
        
        offset1 = 3*nb
        offset2 = offset1 + 4*nb
        offset3 = offset2 + nb

        sol[0:offset1] = r_ddot
        sol[offset1:offset2] = p_ddot
        sol[offset2 : offset3] = self.lagrange_euler
        sol[offset3::] = self.lagrange
        
        max_iter = 30
        iter_count = 0
#        psi = self.sensitivity()
        
        while delta_mag > tol and iter_count < max_iter:
        
            #STEP 1: COMPUTE POS AND VEL USING MOST RECENT ACCELERATIONS
    
            #step level 0 quantities
            r_n = self.solver.pos_step(accel = r_ddot, pos_history = r_h, vel_history = r_dot_h)        
            p_n = self.solver.pos_step(accel = p_ddot, pos_history = p_h, vel_history = p_dot_h)
            
            self._set_generalized_coords(r_n, p_n, level = 0)
            
            #step level 1 quantities
            r_dot_n = self.solver.vel_step(accel = r_ddot, vel_history = r_dot_h )        
            p_dot_n = self.solver.vel_step(accel = p_ddot, vel_history = p_dot_h)
    
            self._set_generalized_coords(r_dot_n, p_dot_n, level = 1)
            
            #STEP 2: COMPUTE RESIDUAL
            
            gn = self.residual()
            
            
            #STEP 3: SOLVE LINEAR SYSTEM FOR CORRECTION
            
            psi = self.sensitivity()
            delta = np.linalg.solve(psi, -gn)
            
            #STEP 4: IMPROVE APPROXIMATE SOLUTION
            
            sol = sol + delta
            
            #set system parameters from new solution
            r_ddot = sol[0:offset1]
            p_ddot = sol[offset1:offset2]
            self.lagrange_euler = column(sol[offset2 : offset3])
            self.lagrange = column(sol[offset3::])
            
            self._set_generalized_coords(r_ddot, p_ddot, level = 2)
            
            #STEP 5: CHECK STOPPING CRITERION. IF SMALL ENOUGH STEP IS DONE
        
            delta_mag = np.linalg.norm(delta)
            
#            if iter_count < 5:
#                print(delta_mag)
            
            iter_count += 1
        
        #use last acceleration data to set position and velocity level coordinates
        #step level 0 quantities
        r_n = self.solver.pos_step(accel = r_ddot, pos_history = r_h, vel_history = r_dot_h)        
        p_n = self.solver.pos_step(accel = p_ddot, pos_history = p_h, vel_history = p_dot_h)
        
        self._set_generalized_coords(r_n, p_n, level = 0)
        
        #step level 1 quantities
        r_dot_n = self.solver.vel_step(accel = r_ddot, vel_history = r_dot_h )        
        p_dot_n = self.solver.vel_step(accel = p_ddot, vel_history = p_dot_h)

        self._set_generalized_coords(r_dot_n, p_dot_n, level = 1)
        
        #snap-shot into system history
        self.history_set()
        
        if history_cnt >= self.order:
            #delete the oldest item out of history if we have enough data
            self.history_delete_oldest()
        
        #set system lagrange parameters
        lagrange_euler = sol[offset2 : offset3]
        lagrange = sol[offset3::]
        
        self._set_lagrange_params(lagrange, lagrange_euler)
        
        

def main():
    
    pass

if __name__ == '__main__': main()