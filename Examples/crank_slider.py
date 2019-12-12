import numpy as np
import simEngine3D as sim
import matplotlib.pyplot as plt

############   CREATE SYSTEM   ############

sys1 = sim.System()

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


crank = sys1.add_body(m = m_crank, J = J_crank) #create rod body

#set position so that crank is aligned with y-axis
crank.set_position([0, L_crank/2, 0])

############   rod   ############

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
rod.set_position([0, L_crank + L_rod + L_slider/2, 0])


############   CREATE JOINTS BETWEEN LIKAGES   ############

#add a revolute joint between the crank and ground
rev1 = sys1.joint_revolute(body_i=sys1.global_frame, sp_i_bar =[0,0,0], 
                           ai_bar = [0,1,0], bi_bar = [0,0,1], 
                           body_j = crank, sq_j_bar = [0, -L_crank/2, 0],
                           cj_bar = [1,0,0])

#add a spherical joint between the crank and rod
sphere1 = sys1.joint_spherical(body_i=crank, sp_i_bar = [0, L_crank/2, 0], 
                               body_j=rod, sq_j_bar=[0, -L_rod/2, 0])
    
#add a revolute joint between rod and slider
rev2 = sys1.joint_revolute(body_i=rod, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar)
    
#add a translational joint between the slider and ground


    
############   ADD RSDA ELEMENT TO DRIVE CRANK-SLIDER   ############

#run with all of RSDA set to zero and observe motion of slider just under gravity    

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

