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

A0 = np.array([[0,  0, 1],
               [0,  1, 0],
               [-1, 0, 0]])

A_theta = np.array([[c,  -s, 0],
               [s,   c, 0],
               [0,   0, 1]])

A1 = A0 @ A_theta

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

rod2.set_orientation_from_A(A0) #set initial orientation of rod to be straight down
rod2.set_position(r = [0, 2*L, -L2])

############   REVOLUTE JOINTS   ############

#add a revolute joint between the pendulum rod and ground
rev1 = sys1.joint_revolute(body_i = sys1.global_frame, sp_i_bar = [0,0,0], ai_bar = [0,1,0], \
                    bi_bar = [0,0,1], body_j = rod1, sq_j_bar = [-L, 0, 0], cj_bar = [0,0,1])

#add a revolute joint between two pendulums
rev2 = sys1.joint_revolute(body_i = rod1, sp_i_bar = [L,0,0], ai_bar = [1,0,0], \
                    bi_bar = [0,1,0], body_j = rod2, sq_j_bar = [-L2, 0, 0], cj_bar = [0,0,1])

############   SET UP SIMULATION PARAMETERS   ############

#solver parameters
h = 0.05
order = 2

#set up solver
sys1.set_solver_order(order)
sys1.set_step_size(h)

t_stop = 10

#parameters to keep track of in simulation
r1 = [rod1.r.flatten()]
omega1 = [rod1.omega.flatten()]

r2 = [rod2.r.flatten()]
omega2 = [rod2.omega.flatten()]

rod1_tip = rod1.create_point([L,0,0])
tip = [rod1_tip.position.flatten()]

time = [0]

############   START STEPPING THROUGH TIME   ############

#initialize problem by solving for initial accelerations
sys1.initialize()

while sys1.t < t_stop:
    
    sys1.step(tol = 1e-3)

    #append desired results
    r1.append(rod1.r.flatten())
    omega1.append(rod1.omega.flatten())
    
    r2.append(rod2.r.flatten())
    omega2.append(rod2.omega.flatten())
    
    tip.append(rod1_tip.position.flatten())
    
    time.append(sys1.t)
    
    if sys1.step_num % 20 == 0:
        print(sys1.t)

########################   POST PROCESSING   ########################

#convert lists to numpy arrays for plotting
t = np.array(time)

#rod 1
r1 = np.vstack(r1)
tip = np.vstack(tip)
theta1 = np.rad2deg(np.arctan2(r1[:,1], -r1[:,2]))
omega1 = np.vstack(omega1)

#rod 2
r2 = np.vstack(r2)
dir_vec = r2 - tip
theta2 = np.rad2deg(np.arctan2(dir_vec[:,1], -dir_vec[:,2]))
omega2 = np.vstack(omega2)

########################   PLOTTING   ########################

def plot(t, y, ylabel, labels, xlabel = 'Time, t [s]'):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(t, y[:,0], label = labels[0])
    ax.plot(t, y[:,1], label = labels[1])
    ax.plot(t, y[:,2], label = labels[2])
    
    ax.set_xlabel(xlabel,  fontsize = 14)
    ax.set_ylabel(ylabel, fontsize = 14)
    
    ax.legend(bbox_to_anchor = (0.0, 1.05, 1, 0.1), loc = 'center', ncol = 3)
    fig.tight_layout()
    
    return fig, ax

plt.close('all')

#plot angles
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t, theta1, label = r'$\theta_1$')
ax.plot(t, theta2, label = r'$\theta_2$')
ax.legend()
ax.set_xlabel('Time, t [s]', fontsize = 14)
ax.set_ylabel(r'Angle, $\theta$ [$\degree$]', fontsize = 14)
fig.tight_layout()

#plot positions
fig1, ax1 = plot(t, r1, ylabel = 'Position [m]', labels = '$x_1$ $y_1$ $z_1$'.split())
fig2, ax2 = plot(t, r2, ylabel = 'Position [m]', labels = '$x_2$ $y_2$ $z_2$'.split())

#plot angular velocities
fig3, ax3 = plot(t, omega1, ylabel = 'Angular Velocity $\omega$ [rad/s]', \
                 labels = '$\omega_{x,1}$ $\omega_{y,1}$ $\omega_{z,1}$'.split())

fig4, ax4 = plot(t, omega2, ylabel = 'Angular Velocity $\omega$ [rad/s]', \
                 labels = '$\omega_{x,2}$ $\omega_{y,2}$ $\omega_{z,2}$'.split())

save_figs = False
if save_figs:
    fig1.savefig('Body1_Pos.svg')
    fig2.savefig('Body2_Pos.svg')
    fig3.savefig('Body1_Omega.svg')
    fig4.savefig('Body2_Omega.svg')