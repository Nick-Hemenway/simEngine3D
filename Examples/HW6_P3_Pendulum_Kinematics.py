###############   IMPORT MODULES   ###############
#add the parent folder to the search path for importing modules
import sys
sys.path.append('..')

import numpy as np
import simEngine3D as sim
import matplotlib.pyplot as plt

############   CREATE SYSTEM   ############

sys1 = sim.System()

############   CREATE PENDULUM WITH REVOLUTE JOINT   ############

m = 1
J = [1,1,1]

L = 2 #length to center of rod = 2 meters

rod = sys1.add_body(m = m, J = J, r = [0,0,-L]) #create rod body

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


############   PROBLEM 3   ############

dt = 1e-3
t_end = 10

time = np.arange(0, t_end + dt/2, dt) #create array of all of the solution times

N = time.shape[0]

#preallocate solution arrays for all of the quantities we want time histories of
r_o = np.zeros( (N, 3) )
r_o_dot = np.zeros( (N, 3) )
r_o_ddot = np.zeros( (N, 3) )

r_q = np.zeros( (N, 3) )
r_q_dot = np.zeros( (N, 3) )
r_q_ddot = np.zeros( (N, 3) )

pQ = rod.create_point([-L,0,0])

for i, t in enumerate(time):
    
    #update the time of the system
    sys1.set_system_time(t)
    
    #solve kinematics
    sys1.solve_position(tol = 1e-6)
    sys1.solve_velocity()
    sys1.solve_acceleration()
    
    #append solutions for each time step to the overall solution arrays to track 
    #the history of the system
    r_o[i,:] =      rod.r.flatten()
    r_o_dot[i,:] =  rod.r_dot.flatten()
    r_o_ddot[i,:] = rod.r_ddot.flatten()
    
    r_q[i,:] =      pQ.position.flatten()
    r_q_dot[i,:] =  pQ.velocity.flatten()
    r_q_ddot[i,:] = pQ.acceleration.flatten()
    
###########   CREATE PLOTS   ############
    
plt.close('all') #close all currently open plots
   
#we need to make a similar plot multiple times so we can create a general plotting 
#function that takes in the values to plot as well as the level of the derivative
#level = 0 --> position
#level = 1 --> velocity
#level = 2 --> acceleration

def make_plot(t, M, level = 0):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    y_labels = ['Position [m]', 'Velocity [m/s]', 'Acceleration [m/$\mathrm{s^2}$]']
    y_label = y_labels[level]
    
    label_padding = [('$', '$'), ('$\dot{', '}$'), ('$\ddot{', '}$')]
    lpad, rpad = label_padding[level][0], label_padding[level][1]
    
    ax.plot(t, M[:,0], label = lpad + 'x' + rpad)
    ax.plot(t, M[:,1], label = lpad + 'y' + rpad)
    ax.plot(t, M[:,2], label = lpad + 'z' + rpad)
    
    ax.set_xlabel('Time, t [s]', fontsize = 14)
    ax.set_ylabel(y_label, fontsize = 14)
    ax.legend(bbox_to_anchor = (0.0, 1.05, 1, 0.1), loc = 'center', ncol = 3)
    
    fig.tight_layout()
    
    return fig

#plots for point o
fig1 = make_plot(time, r_o, level = 0)
fig2 = make_plot(time, r_o_dot, level = 1)
fig3 = make_plot(time, r_o_ddot, level = 2)

#plots for point q
fig4 = make_plot(time, r_q, level = 0)
fig5 = make_plot(time, r_q_dot, level = 1)
fig6 = make_plot(time, r_q_ddot, level = 2)

save_figs = False
if save_figs:
    fig1.savefig('Position_O.svg')
    fig2.savefig('Velocity_O.svg')
    fig3.savefig('Acceleration_O.svg')
    fig4.savefig('Position_Q.svg')
    fig5.savefig('Velocity_Q.svg')
    fig6.savefig('Acceleration_Q.svg')

