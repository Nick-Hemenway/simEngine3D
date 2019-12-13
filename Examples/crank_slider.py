import numpy as np
import simEngine3D as sim
import matplotlib.pyplot as plt

############   CREATE SYSTEM   ############

sys1 = sim.System()

############   INPUT VARIABLES   ############

motor_torque = 0  #N-m
k = 0             #spring constant N-m/rad
c = 0.6           #damping coeff N-m/(rad/s)

############   COMMON PROPERTIES   ############

rho = 7800 #kg/m^3
w = 0.05 #thickness into page of all of links

############   CREATE CRANK   ############

L_crank = 0.2 #length to center of rod = 2 meters
m_crank = rho*L_crank*w**2

Ixx_crank = (1/12)*m_crank*(L_crank**2 + w**2)
Iyy_crank = (1/12)*m_crank*(w**2 + w**2)
Izz_crank = Ixx_crank

J_crank = [Ixx_crank, Iyy_crank, Izz_crank]


crank = sys1.add_body(m = m_crank, J = J_crank) #create crank body

#set position so that crank is aligned with y-axis
crank.set_position([0, L_crank/2, 0])

############   ROD   ############

L_rod = 3*L_crank

m_rod = rho*L_rod*w**2

Ixx_rod = (1/12)*m_rod*(L_rod**2 + w**2)
Iyy_rod = (1/12)*m_rod*(w**2 + w**2)
Izz_rod = Ixx_rod

J_rod = [Ixx_rod, Iyy_rod, Izz_rod]

rod = sys1.add_body(m = m_rod, J = J_rod) #create rod body

#set rod to be aligned with y-axis
rod.set_position([0, L_crank + L_rod/2, 0])

############   SLIDER   ############

#create slider as a cylinder

L_slider = 0.1
r_slider = 0.05

m_slider = rho*L_slider*np.pi*r_slider**2

Ixx_slider = (1/4)*m_slider*r_slider**2 + (1/12)*m_slider*L_slider**2
Iyy_slider = 0.5*m_slider*r_slider**2
Izz_slider = Ixx_slider

J_slider = [Ixx_slider, Iyy_slider, Izz_slider]

slider = sys1.add_body(m = m_slider, J = J_slider) #create rod body

#set rod to be aligned with y-axis
slider.set_position([0, L_crank + L_rod, 0])


############   CREATE JOINTS BETWEEN LIKAGES   ############

#note that there are three bodies with 6 DOF each for a total of 18 DOF. To 
#simulate the crank slider there can only be one DOF. simEngine3D can not handle
#redundant constraints and thus we must choose our joints that we connect each 
#link with carefully.

# 2 revolute = 2*5 = 10
# 1 spherical = 3
# 1 cylindrical = 4

#using this combination leaves one degree of freedom which is what we need to 
#solve the dynamics

#add a revolute joint between the crank and ground
rev1 = sys1.joint_revolute(body_i=sys1.global_frame, sp_i_bar =[0,0,0], 
                           ai_bar = [0,1,0], bi_bar = [0,0,1], 
                           body_j = crank, sq_j_bar = [0, -L_crank/2, 0],
                           cj_bar = [1,0,0])

#add a spherical joint between the crank and rod
sphere1 = sys1.joint_spherical(body_i=crank, sp_i_bar = [0, L_crank/2, 0], 
                               body_j=rod, sq_j_bar=[0, -L_rod/2, 0])
    
#add a revolute joint between rod and slider
rev2 = sys1.joint_revolute(body_i=rod, sp_i_bar = [0, L_rod/2, 0], 
                           ai_bar=[0,1,0], bi_bar=[0,0,1],
                           body_j=slider, sq_j_bar=[0,0,0], cj_bar=[1,0,0])


#add a cylindrical joint between the slider and ground
cyl = sys1.joint_cylindrical(body_i=sys1.global_frame, sp_i_bar=[0,0,0],
                             ai_bar=[1,0,0], bi_bar=[0,0,1], body_j=slider,
                             sq_j_bar=[0,0,0], cj_bar=[0,1,0])
    

############   ADD RSDA ELEMENT TO DRIVE CRANK-SLIDER   ############

theta_0 = 0 #initial spring postion in degrees

actuator = lambda theta, theta_dot, t: motor_torque

RSDA = sys1.torque_RSDA(body_i = sys1.global_frame, ai_bar = [1,0,0], bi_bar = [0,1,0], 
                 body_j = crank, aj_bar = [1,0,0], bj_bar = [0,1,0], k=k, 
                 theta_0 = theta_0, c = c, h = actuator) 

############   SET UP SIMULATION PARAMETERS   ############

#solver parameters
h = 0.02
order = 2

#set up solver
sys1.set_solver_order(order)
sys1.set_step_size(h)

t_stop = 10

#parameters to keep track of in simulation
r_crank = [crank.r.flatten()]
omega_crank = [crank.omega.flatten()]
r_rod = [rod.r.flatten()]
r_slider = [slider.r.flatten()]

p = crank.create_point([0,L_crank/2,0])
point = [p.position.flatten()]


time = [0]

############   START STEPPING THROUGH TIME   ############

#initialize problem by solving for initial accelerations
sys1.initialize()

while sys1.t < t_stop:
    
    sys1.step(tol = 1e-3)

    #append desired results
    r_crank.append(crank.r.flatten())
    omega_crank.append(crank.omega.flatten())
    r_rod.append(rod.r.flatten())
    r_slider.append(slider.r.flatten())
    point.append(p.position.flatten())
    
    time.append(sys1.t)
    
    sys1.print_time(num_steps = 20)

########################   POST PROCESSING   ########################

#convert lists to numpy arrays for plotting
t = np.array(time)

r_crank = np.vstack(r_crank)
omega_crank = np.vstack(omega_crank)
r_rod = np.vstack(r_rod)
r_slider = np.vstack(r_slider)
point = np.vstack(point)

#save data to text file for making animation
out_arr = np.hstack((t.reshape(-1,1), point, r_slider))
np.savetxt('Crank_Slider_Data.txt', out_arr, header = 't, py, py, pz, sx, sy, sz')


########################   PLOTTING   ########################

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel('Time, $t$ [s]')
ax.set_ylabel('Crank Angular Velocity, $\omega$ [rad/s]')

ax.plot(t, omega_crank[:,0])
fig.tight_layout()

########################   ANIMATION   ########################

animate = False

if animate:

    print('\nDone with simulation.... Creating animation\n')

    import crank_slider_animation
    crank_slider_animation.animate(fname = 'Crank_Slider_Damping.mp4')

    print('\nDone')    


























