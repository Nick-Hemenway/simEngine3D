import numpy as np
import simEngine3D as sim
import appliedLoads as loads
import rigidbody
import matplotlib.pyplot as plt

body_i = rigidbody.RigidBody(m=4, J=[1,1,1], r = [0,0,0], r_dot = [-1,0,0])
body_j = rigidbody.RigidBody(m=2, J=[1,1,1], r = [1,0,0], r_dot = [1,0,0])

########################   TEST TSDA   ########################

sp_i_bar = [1,1,1]
sq_j_bar = [1,1,1]

k = 10
l_0 = 1.5 #set l_0 to 1 to see effects of only damping
c = 2
h = lambda lij, lij_dot, t: 0

F = loads.TSDA(body_i, sp_i_bar, body_j, sq_j_bar, k, l_0, c, h)

eij, lij, lij_dot = F._calc_info()
Fp, Fq = F(10)

# print(f'\nUnit vector pointing from P to Q:\n\n {eij}')
# print(f'\nDistance between P and Q:\n\n {lij}')
# print(f'\nTime rate of change of length between P and Q:\n\n {lij_dot}')
# print(f'\nForce acting on point P:\n\n {Fp}')
# print(f'\nForce acting on point Q:\n\n {Fq}')


########################   TEST RSDA   ########################

body_i = rigidbody.RigidBody(m=4, J=[1,1,1], r = [0,0,0], r_dot = [-1,0,0])
body_j = rigidbody.RigidBody(m=2, J=[1,1,1], r = [1,0,0], r_dot = [1,0,0])

body_i.set_ang_vel_from_omega_bar([0,0,-1])
body_j.set_ang_vel_from_omega_bar([0,0,1])

k = 10
theta_0 = 0
c = 2
h = 0

ai_bar = [0,0,1]
aj_bar = [0,0,1]

bi_bar = [1,0,0]
bj_bar = [1,0,0]

T = loads.RSDA(body_i, ai_bar, bi_bar, body_j, aj_bar, bj_bar, k, theta_0, c, h)

theta = np.rad2deg(T._calc_theta_ij())
theta_dot = T._calc_theta_ij_dot()

def calc_A(angle):
    
    angle = np.deg2rad(angle)
    c = np.cos(angle)
    s = np.sin(angle)
    
    A = np.array([[ c, -s, 0 ],
                   [ s,  c, 0 ],
                   [ 0,  0, 1 ]])
    
    return A

#generate profile in time of what the angle theta should be
offset = 0
t1 = np.linspace(0, -700, 100) + offset
t2 = np.linspace(-700, 600, 100) + offset
t3 = np.linspace(600, -200, 100) + offset
theta_prof = np.hstack((t1,t2,t3))

measured_angle = []

for angle in theta_prof:
    
    A = calc_A(angle)
    body_j.set_orientation_from_A(A)
    
    measured_angle.append(np.rad2deg(T._calc_theta_ij()))
    
    
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(measured_angle)
ax.set_xlabel('Step')
ax.set_ylabel('Measured Angle')
    
    
    

