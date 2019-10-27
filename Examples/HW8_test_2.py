###############   IMPORT MODULES   ###############
#add the parent folder to the search path for importing modules
import sys
sys.path.append('..')

import numpy as np
import simEngine3D as sim
import matplotlib.pyplot as plt

import time as timer

############   CREATE SYSTEM   ############

sys1 = sim.System()

############   CREATE PENDULUM WITH REVOLUTE JOINT   ############


L = 2 #length to center of rod = 2 meters
rho = 7800 #kg/m^3
w = 0.05 #m

m = rho*2*L*w**2

Ixx = (1/12)*m*(w**2 + w**2)
Iyy = (1/12)*m*((2*L)**2 + w**2)
Izz = Iyy

J = [Ixx, Iyy, Izz]

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

driving_con = sys1.constraint_DP1(body_i = sys1.global_frame, ai_bar = [0,0,-1], body_j = rod, aj_bar = [0,1,0], constraint_func = func)

#############   PROBLEM 3   ############

dt = 1e-3
t_end = np.pi

time = np.arange(0, t_end + dt/2, dt) #create array of all of the solution times

N = time.shape[0]

Torque = np.zeros((N,3))
Force_reaction = np.zeros((N,3))

#preallocate solution arrays for all of the quantities we want time histories of
r_o = np.zeros( (N, 3) )
r_o_dot = np.zeros( (N, 3) )
r_o_ddot = np.zeros( (N, 3) )

for i, t in enumerate(time):
    
    #update the time of the system
    sys1.set_system_time(t)
    
    #solve kinematics for time step
    sys1.solve_kinematics(tol = 1e-6, update_jacobian_every = 3)
    
    #solve inverse dynamics for time step for lagrange multipliers
    lam = sys1.solve_inverse_dynamics()
    
    #solve for joint reaction forces
    phi_ri, phi_rj = rev.partial_r()
    F = -phi_rj.T @ lam[0:5]
    
    #solve for generalized driving torque
    phi_pi, phi_pj = driving_con.partial_p()
    E = rod.E
    
    T = -0.5*phi_pj @ E.T *lam[5] #slide 16 of lecture 9 for going from phi to pi
    
    #append reaction forces/torques to array
    Torque[i,:] = T.flatten()
    Force_reaction[i,:] = F.flatten()
    
    #append kinematic solutions for each time step to the overall solution arrays 
    #to track the history of the system
    r_o[i,:] =      rod.r.flatten()
    r_o_dot[i,:] =  rod.r_dot.flatten()
    r_o_ddot[i,:] = rod.r_ddot.flatten()
    

###########   CREATE PLOTS   ############
    
plt.rc('font', family='serif') 
plt.rc('font', size = 14)
plt.rc('font', serif='Times New Roman') 
plt.rc('mathtext', fontset = 'cm')

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
#
#fig1.savefig('Plots/Position_O.svg')
#fig2.savefig('Plots/Velocity_O.svg')
#fig3.savefig('Plots/Acceleration_O.svg')
    
fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
ax7.set_xlabel('Time, t [s]', fontsize = 14)
ax7.set_ylabel(r'Torque, $\tau$ [Nm]', fontsize = 14)

ax7.plot(time, Torque[:,0])
fig7.tight_layout()
#fig7.savefig('Plots/Driving_Torque.svg')
#fig7.savefig('Plots/Driving_Torque.png')

fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
ax8.set_xlabel('Time, t [s]', fontsize = 14)
ax8.set_ylabel(r'Reaction Force, F [N]', fontsize = 14)

ax8.plot(time, Force_reaction[:,1], label = '$F_y$')
ax8.plot(time, Force_reaction[:,2], label = '$F_z$')
ax8.legend(fontsize = 12)
fig8.tight_layout()
#fig8.savefig('Plots/Reaction_Forces.svg')
#fig8.savefig('Plots/Reaction_Forces.png')

