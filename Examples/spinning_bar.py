###############   IMPORT MODULES   ###############
#add the parent folder to the search path for importing modules
import sys
sys.path.append('..')

import numpy as np
import matplotlib.pyplot as plt
import simEngine3D as sim

#%% ############   CREATE SYSTEM   ############

sys1 = sim.System()

#%% ############   CREATE ROD WITH PIVOT AT CENTER   ############

L = 1 #Total length of rod
rho = 7800 #kg/m^3
w = 0.05 #m

############   ROD 1   ############

m = rho*L*w**2

Ixx = (1/12)*m*(L**2 + w**2)
Iyy = (1/12)*m*(w**2 + w**2)
Izz = Ixx

J = [Ixx, Iyy, Izz]

rod = sys1.add_body(m = m, J = J) #create rod body

#initial conditions
theta_init = 90
theta_dot_init = 0

A0 = sim.utility.rotate_x(theta_init)
rod.set_orientation_from_A(A0)
rod.set_ang_vel_from_omega([theta_dot_init, 0, 0])

#%% ############   ADD RSDA ELEMENT   ############

theta_0 = 0 #initial spring postion in degrees
k = 10 #N/m
c = 2 #damping N/(m/s)
T_act = 0
actuator = lambda theta, theta_dot, t: -T_act

RSDA = sys1.torque_RSDA(body_i = sys1.global_frame, ai_bar = [1,0,0], bi_bar = [0,1,0], 
                 body_j = rod, aj_bar = [1,0,0], bj_bar = [0,1,0], k=k, 
                 theta_0 = theta_0, c = c, h = actuator)

#%% ############   START STEPPING THROUGH TIME   ############

#solver parameters
h = 0.01 #time step
order = 2 #solver order

#set up solver
sys1.set_solver_order(order)
sys1.set_step_size(h)

#initialize problem by solving for initial accelerations
sys1.initialize()

#variables of interest
omega = [rod.omega.flatten()]
theta = [np.rad2deg(RSDA.calc_theta())]

t_stop = 10
time = [0]

while sys1.t < t_stop:
    
    sys1.step(tol = 1e-3)
    omega.append(rod.omega.flatten())
    theta.append(np.rad2deg(RSDA.calc_theta()))
    
    time.append(sys1.t)
    sys1.print_time(num_steps=20)
    
        
#%% ############   ANALYTICAL SOLUTION   ############     
        
t_analytic = np.linspace(0, t_stop, 100)

omega_n = np.sqrt(k/Ixx)
zeta = c/(2*Ixx*omega_n)

#case 1: Underdamped
if zeta < 1:
    omega_d = omega_n*np.sqrt(1-zeta)
    
    theta_analytic = np.exp(-zeta*omega_n*t_analytic)*(theta_init*np.cos(omega_d*t_analytic) + (theta_dot_init + zeta*omega_n*theta_init)/omega_d*np.sin(omega_d*t_analytic))
        
# omega_analytic = T_act/Ixx*t_analytic
        
#%% ############   POST PROCESSING AND PLOTTING   ############        
        
omega = np.vstack(omega)
theta = np.array(theta)

fig = plt.figure()
ax = fig.add_subplot(111)

# ax.plot(time, omega[:,0], label = 'Simulated')
ax.plot(time, theta, label = 'Simulated')
ax.plot(t_analytic, theta_analytic, ls = '--', label = 'Analytic')

ax.set_xlabel('Time, t [s]')
ax.set_ylabel('Angular Velocity, $\omega$ [rad/s]')
# ax.legend()
fig.tight_layout()