###############   IMPORT MODULES   ###############
#add the parent folder to the search path for importing modules
import sys
sys.path.append('..')

import numpy as np
import simEngine3D as sim

############   CREATE SYSTEM   ############

sys1 = sim.System()

############   CREATE PENDULUM WITH REVOLUTE JOINT   ############

L = 2 #length to center of rod = 2 meters

rod = sys1.add_body(r = [0,0,-L]) #create rod body

A0 = np.array([[0,  0, 1],
               [0,  1, 0],
               [-1, 0, 0]])


rod.set_orientation_from_A(A0) #set initial orientation of rod to be straight down

#add a revolute joint between the pendulum rod and ground
rev = sys1.joint_revolute(body_i = sys1.global_frame, sp_i_bar = [0,0,0], ai_bar = [0,1,0], \
                    bi_bar = [0,0,1], body_j = rod, sq_j_bar = [-L, 0, 0], cj_bar = [0,0,1])

############   CREATE DRIVING CONSTRAINT   ############

theta =      lambda t:  (np.pi/4) * np.cos(2*t)
theta_dot =  lambda t: -(np.pi/2) * np.sin(2*t)
theta_ddot = lambda t:  -np.pi    * np.cos(2*t)

#we know theta as a function of time but we must convert that to the proper 
#to input for our DP1 constraint. the dot product of two vectors is: 
# a_vec * b_vec  = |a||b| cos(theta). If a and b are both unit vectors then 
# f(t) = cos(theta)

offset = np.pi/2 #angle between two vectors of driving constraint at theta = 0

f =      lambda t:  np.cos(theta(t) + offset) 
f_dot =  lambda t: -np.sin(theta(t) + offset) * theta_dot(t)
f_ddot = lambda t: -np.cos(theta(t) + offset) * theta_dot(t)**2 - np.sin(theta(t) + offset)*theta_ddot(t)

func = sys1.constraint_func(f = f, f_prime = f_dot, f_pprime = f_ddot) #create constraint function object for DP1 constraint

dp1 = sys1.constraint_DP1(body_i = sys1.global_frame, ai_bar = [0,0,-1], body_j = rod, aj_bar = [0,1,0], constraint_func = func)

############   PROBLEM 2   ############

sys1.set_system_time(t = 0)
sys1.solve_position() #make sure the initial configuration is correct

phi = sys1.phi()
J = sys1.jacobian()
nu = sys1.vel_rhs()
gamma = sys1.accel_rhs()

print(f'\nPhi: \n\n {phi}')
print(f'\nJacobian: \n\n {J}')
print(f'\nNu: \n\n {nu}')
print(f'\nGamma: \n\n {gamma}')


