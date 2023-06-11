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
h = 0.05
order = 2

#set up solver
sys1.set_solver_order(order)
sys1.set_step_size(h)

#initialize problem by solving for initial accelerations
sys1.initialize()

r_o = [rod.r.flatten()]
omega = [rod.omega.flatten()]
torque = [rev.reaction_torque(body = 'j').flatten()]
T_pin = [rev.reaction(body = 'j')[0].flatten()]

#velocity constraint violation
q_dot = np.vstack( (rod.r_dot, rod.p_dot) )
phi_q = np.hstack( (rev.partial_r()[1], rev.partial_p()[1]) )
nu = rev.vel_rhs(sys1.t)

violation = [np.linalg.norm( (phi_q @ q_dot - nu).flatten() )]

t_stop = 10
time = [0]

start = timer.time()

while sys1.t < t_stop:
    
    sys1.step(tol = 1e-3)

    #append desired results
    r_o.append(rod.r.flatten())
    omega.append(rod.omega.flatten())
    torque.append(rev.reaction_torque(body = 'j').flatten())
    T_pin.append(rev.reaction(body = 'j')[0].flatten())
    
    q_dot = np.vstack( (rod.r_dot, rod.p_dot) )
    phi_q = np.hstack( (rev.partial_r()[1], rev.partial_p()[1]) )
    nu = rev.vel_rhs(sys1.t)
   
    violation.append(np.linalg.norm( (phi_q @ q_dot - nu).flatten() ))
    
    time.append(sys1.t)

    if sys1.step_num % 20 == 0:
        print(sys1.t)
        
stop = timer.time()

print(f'Total Time {stop - start}')
        
t = np.array(time)

pos = np.vstack(r_o)
vel = np.vstack(omega)
T = np.vstack(torque)
T_pin = np.vstack(T_pin)

############   PLOTTING   ############

plt.close('all')

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
    
    
#plot position
fig1, ax1 = plot(t, pos, ylabel = 'Position [m]', labels = '$x$ $y$ $z$'.split())

#plot velocity
fig2, ax2 = plot(t, vel, ylabel = 'Angular Velocity, \omega [rad/s]', \
                 labels = '$\omega_{x}$ $\omega_{y}$ $\omega_{z}$'.split())
fig3, ax3 = plot(t, T, ylabel = r'Torque, $\tau$ [N-m]', \
                 labels = '$T_{x}$ $T_{y}$ $T_{z}$'.split())
    
fig4, ax4 = plot(t, T_pin, ylabel = r'Pin Torque, $\tau$ [N-m]', \
                 labels = '$T_{x}$ $T_{y}$ $T_{z}$'.split())
    
#plot velocity constraint violation magnitude

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(t, violation)

ax.set_xlabel('Time, t [s]',  fontsize = 14)
ax.set_ylabel('Velocity Violation', fontsize = 14)

fig.tight_layout()

save_figs = False
if save_figs:
    fig1.savefig('Position.svg')
    fig2.savefig('Omega.svg')
    fig3.savefig('Torque.svg')
    fig.savefig('Part1_Velocity_Violation.svg')