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

############   ROD 1   ############

m1 = rho*2*L*w**2

Ixx_1 = (1/12)*m1*(w**2 + w**2)
Iyy_1 = (1/12)*m1*((2*L)**2 + w**2)
Izz_1 = Iyy_1

J1 = [Ixx_1, Iyy_1, Izz_1]

theta_0 = np.pi/2
c = np.cos(theta_0)
s = np.sin(theta_0)

rod1 = sys1.add_body(m = m1, J = J1) #create rod body

A_down = np.array([[0,  0, 1],
               [0,  1, 0],
               [-1, 0, 0]])

A_theta = np.array([[c,  -s, 0],
               [s,   c, 0],
               [0,   0, 1]])

A1 = A_down @ A_theta

rod1.set_orientation_from_A(A1) #set initial orientation of rod to be straight down
rod1.set_position(r = [0, L*s, -L*c])

############   ROD 2   ############

L2 = L/2

m2 = rho*2*(L2)*w**2

Ixx_2 = (1/12)*m2*(w**2 + w**2)
Iyy_2 = (1/12)*m2*((2*L2)**2 + w**2)
Izz_2 = Iyy_2

J2 = [Ixx_2, Iyy_2, Izz_2]

rod2 = sys1.add_body(m = m2, J = J2) #create rod body

rod2.set_orientation_from_A(A_down) #set initial orientation of rod to be straight down
rod2.set_position(r = [0, 2*L, -L2])


#add a revolute joint between the pendulum rod and ground
rev1 = sys1.joint_revolute(body_i = sys1.global_frame, sp_i_bar = [0,0,0], ai_bar = [0,1,0], \
                    bi_bar = [0,0,1], body_j = rod1, sq_j_bar = [-L, 0, 0], cj_bar = [0,0,1])

#add a revolute joint between two pendulums
rev2 = sys1.joint_revolute(body_i = rod1, sp_i_bar = [L,0,0], ai_bar = [1,0,0], \
                    bi_bar = [0,1,0], body_j = rod2, sq_j_bar = [-L2, 0, 0], cj_bar = [0,0,1])

sys1.initialize()

############   START STEPPING THROUGH TIME   ############

#solver parameters
h = 0.01
order = 2

#set up solver
sys1.set_solver_order(order)
sys1.set_step_size(h)

#initialize problem by solving for initial accelerations
sys1.initialize()

r_o1 = [rod1.r.flatten()]
v_o1 = [rod1.r_dot.flatten()]

r_o2 = [rod2.r.flatten()]
v_o2 = [rod2.r_dot.flatten()]

t_stop = 10
time = [0]

while sys1.t < t_stop:
    
    sys1.step(tol = 1e-3)

    #append desired results
    r_o1.append(rod1.r.flatten())
    v_o1.append(rod1.r_dot.flatten())
    
    r_o2.append(rod2.r.flatten())
    v_o2.append(rod2.r_dot.flatten())
    
    time.append(sys1.t)
    
    if sys1.step_num % 20 == 0:
        print(sys1.t)

t = np.array(time)
pos1 = np.vstack(r_o1)
vel1 = np.vstack(v_o1)
pos2 = np.vstack(r_o2)
vel2 = np.vstack(v_o2)


############   PLOTTING   ############

plt.close('all')

#body 1 position
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(t, pos1[:,0], label = '$x_1$')
ax.plot(t, pos1[:,1], label = '$y_1$')
ax.plot(t, pos1[:,2], label = '$z_1$')

ax.set_xlabel('Time, t [s]',  fontsize = 14)
ax.set_ylabel('Position [m]', fontsize = 14)

ax.legend(bbox_to_anchor = (0.0, 1.05, 1, 0.1), loc = 'center', ncol = 3)
fig.tight_layout()

#body 2 position
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.plot(t, pos2[:,0], label = '$x_2$')
ax1.plot(t, pos2[:,1], label = '$y_2$')
ax1.plot(t, pos2[:,2], label = '$z_2$')

ax1.set_xlabel('Time, t [s]',  fontsize = 14)
ax1.set_ylabel('Position [m]', fontsize = 14)

ax1.legend(bbox_to_anchor = (0.0, 1.05, 1, 0.1), loc = 'center', ncol = 3)
fig1.tight_layout()

#body 1 velocity
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

ax2.plot(t, vel1[:,0], label = '$\dot{x}_1$')
ax2.plot(t, vel1[:,1], label = '$\dot{y}_1$')
ax2.plot(t, vel1[:,2], label = '$\dot{z}_1$')

ax2.set_xlabel('Time, t [s]',  fontsize = 14)
ax2.set_ylabel('Velocity [m/s]', fontsize = 14)

ax2.legend(bbox_to_anchor = (0.0, 1.05, 1, 0.1), loc = 'center', ncol = 3)
fig2.tight_layout()

#body 2 velocity
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

ax3.plot(t, vel2[:,0], label = '$\dot{x}_2$')
ax3.plot(t, vel2[:,1], label = '$\dot{y}_2$')
ax3.plot(t, vel2[:,2], label = '$\dot{z}_2$')

ax3.set_xlabel('Time, t [s]',  fontsize = 14)
ax3.set_ylabel('Velocity [m/s]', fontsize = 14)

ax3.legend(bbox_to_anchor = (0.0, 1.05, 1, 0.1), loc = 'center', ncol = 3)
fig3.tight_layout()

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)

#theta = np.
#ax4.plot(t, L_mag)

