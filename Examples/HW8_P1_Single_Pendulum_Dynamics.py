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


L = 2 #length to center of rod = 2 meters
rho = 7800 #kg/m^3
w = 0.05 #m

m = rho*2*L*w**2

Ixx = (1/12)*m*(w**2 + w**2)
Iyy = (1/12)*m*((2*L)**2 + w**2)
Izz = Iyy

J = [Ixx, Iyy, Izz]

theta_0 = np.pi/4
c = np.cos(theta_0)
s = np.sin(theta_0)

rod = sys1.add_body(m = m, J = J) #create rod body

A1 = np.array([[0,  0, 1],
               [0,  1, 0],
               [-1, 0, 0]])

A2 = np.array([[c,  -s, 0],
               [s,   c, 0],
               [0,   0, 1]])

A0 = A1 @ A2

rod.set_orientation_from_A(A0) #set initial orientation of rod to be straight down
rod.set_position(r = [0, L*s, -L*c])

#add a revolute joint between the pendulum rod and ground
rev = sys1.joint_revolute(body_i = sys1.global_frame, sp_i_bar = [0,0,0], ai_bar = [0,1,0], \
                    bi_bar = [0,0,1], body_j = rod, sq_j_bar = [-L, 0, 0], cj_bar = [0,0,1])

############   START STEPPING THROUGH TIME   ############

#solver parameters
h = 0.01
order = 2

#set up solver
sys1.set_solver_order(order)
sys1.set_step_size(h)

#initialize problem by solving for initial accelerations
sys1.initialize()

r_o = [rod.r.flatten()]
v_o = [rod.r_dot.flatten()]

t_stop = 10
time = [0]

while sys1.t < t_stop:
    
    sys1.step(tol = 1e-3)

    #append desired results
    r_o.append(rod.r.flatten())
    v_o.append(rod.r_dot.flatten())
    
    time.append(sys1.t)

t = np.array(time)
pos = np.vstack(r_o)
vel = np.vstack(v_o)


############   PLOTTING   ############

plt.close('all')

#plot position
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(t, pos[:,0], label = '$x$')
ax.plot(t, pos[:,1], label = '$y$')
ax.plot(t, pos[:,2], label = '$z$')

ax.set_xlabel('Time, t [s]',  fontsize = 14)
ax.set_ylabel('Position [m]', fontsize = 14)

ax.legend(bbox_to_anchor = (0.0, 1.05, 1, 0.1), loc = 'center', ncol = 3)
fig.tight_layout()

#plot velocity
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.plot(t, vel[:,0], label = '$\dot{x}$')
ax1.plot(t, vel[:,1], label = '$\dot{y}$')
ax1.plot(t, vel[:,2], label = '$\dot{z}$')

ax1.set_xlabel('Time, t [s]',  fontsize = 14)
ax1.set_ylabel('Velocity [m/s]', fontsize = 14)

ax1.legend(bbox_to_anchor = (0.0, 1.05, 1, 0.1), loc = 'center', ncol = 3)
fig1.tight_layout()

