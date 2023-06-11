###############   IMPORT MODULES   ###############
#add the parent folder to the search path for importing modules
import sys
sys.path.append('..')

import numpy as np
import matplotlib.pyplot as plt
import simEngine3D as sim

#%% ############   CREATE SYSTEM   ############

sys1 = sim.System()
sys1.set_gravity([0,0,0])

#%% ############   CREATE SPHERICAL MASS   ############

r = 0.050 #50 mm
V = (4/3)*np.pi*r**3 #volume of a sphere
rho = 7800 #density of steel [kg/m^3]

m = rho*V #mass in kg
Ixx = (2/5)*m*r**2 #inertia for solid sphere

J = [Ixx, Ixx, Ixx] #inertia is same in all directions

bob = sys1.add_body(m = m, J = J) #create spherical bob

#%% ############   ADD SPRING AND DAMPER   ############

l0 = 0 #unstretched spring length
k = 100 #N/m
c = 5 #damping N/(m/s)
# c = 0
actuator = lambda l, l_dot, t: 0

resting_pos = l0 + m*9.80665/k #position of mass at rest accounting for gravity

#initial conditions
h0 = 0.1
v0 = 0

bob.set_position([0, 0, h0]) #-resting_pos + h0]) #bob starts sitting at springs unstretched position
bob.set_vel([0, 0, v0])

sys1.force_TSDA(sys1.global_frame, [0,0,0], bob, [0,0,0], k, l0, c, actuator)

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
r = [bob.r.flatten()]
v = [bob.r_dot.flatten()]

t_stop = 6
time = [0]

while sys1.t < t_stop:
    
    sys1.step(tol = 1e-3)
    r.append(bob.r.flatten())
    v.append(bob.r_dot.flatten())
    
    time.append(sys1.t)
    sys1.print_time(num_steps = 20)
        
#%% ############   ANALYTICAL SOLUTION   ############     
        
def MSD(m,k,c,x0,v0,t_arr):
    
    wn = np.sqrt(k/m)
    zeta = c/(2*m*wn)
    
    #underdamped
    if zeta < 1:
        
        wd = wn*np.sqrt(1-zeta)
        
        sin_terms = ((v0 + zeta*wn*x0)/wd)*np.sin(wd*t_arr) + x0*np.cos(wd*t_arr)
        
        y = np.exp(-zeta*wn*t_arr)*sin_terms
        
        return y    
    
t_analytic = np.linspace(0, t_stop, 100)

z_analytic = MSD(m,k,c,h0,v0,t_analytic)
        
# omega_n = np.sqrt(k/m)
# zeta = c/(2*m*omega_n)

# #case 1: Underdamped
# if zeta < 1:
#     omega_d = omega_n*np.sqrt(1-zeta)
    
#     z_analytic = np.exp(-zeta*omega_n*t_analytic)*(h0*np.cos(omega_d*t_analytic) + (v0 + zeta*omega_n*h0)/omega_d*np.sin(omega_d*t_analytic))

        
#%% ############   POST PROCESSING AND PLOTTING   ############        
        
r = np.vstack(r)
# r += resting_pos #add resting pos so that it oscillates about zero

v = np.vstack(v)
v_max = v.max()

fig = plt.figure()
ax = fig.add_subplot(111)

# ax.axhline(-resting_pos, ls = '--', color = 'k', zorder = -1)
# ax.axhline(0, ls = '-', lw = 1, color = 'k', zorder = -1)
ax.plot(time, r[:,2], label = 'Simulated')
ax.plot(t_analytic, z_analytic, ls = '--', label = 'Analytic')

ax.set_xlabel('Time, t [s]')
ax.set_ylabel('Bob Position, $z$ [m]')
ax.legend()
fig.tight_layout()