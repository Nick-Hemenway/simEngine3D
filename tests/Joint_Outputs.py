###############   IMPORT MODULES   ###############
#add the parent folder to the search path for importing modules
import sys
sys.path.append('../..')

import numpy as np
import simEngine3D as sim

###############   CREATE A SYSTEM   ###############

sys1 = sim.System()

###############   INPUT PARAMETERS FOR EACH BODY   ###############
ri = [8,6,-3]
ri_dot = [7,8,9]

temp = np.array([4,3,-5,1])
pi = temp/np.linalg.norm(temp) 

pi_dot = np.array([-0.2,1.3,3.4,0])
pi_dot[-1] = -np.dot(pi, pi_dot)/pi[-1]
pi_dot = pi_dot/np.linalg.norm(pi_dot)

rj = [-0.5,1.6,-6.3]
rj_dot = [11, 12, 13]

temp = np.array([3.3, -4, 5.1, 6])
pj = temp/np.linalg.norm(temp) 

pj_dot = np.array([0.6, -3.7, 5.1, 0])
pj_dot[-1] = -np.dot(pj, pj_dot)/pj[-1]
pj_dot = pj_dot/np.linalg.norm(pj_dot)

################   CREATE BODIES   ###############

body_i = sys1.add_body(ri, pi, ri_dot, pi_dot)
body_j = sys1.add_body(rj, pj, rj_dot, pj_dot)

###############   VECTORS FOR CONSTRAINTS   ###############

ai_bar = [-1.2, 1 ,0.3]
aj_bar = [1.2, 4.5, 3.1]
bi_bar = [-2, 4.2, -1.6]
cj_bar = [0.3, 0.4, -6]
sp_i_bar = [0.1, -0.3, 6.0]
sq_j_bar = [0.2, -1.0, 1.5]

###############   CREATE JOINTS   ###############

SJ = sys1.joint_spherical(body_i, sp_i_bar, body_j, sq_j_bar)
UJ = sys1.joint_universal(body_i, sp_i_bar, ai_bar, body_j, sq_j_bar, aj_bar)
CJ = sys1.joint_cylindrical(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar)
RJ = sys1.joint_revolute(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, cj_bar)
TJ = sys1.joint_translational(body_i, sp_i_bar, ai_bar, bi_bar, body_j, sq_j_bar, aj_bar, cj_bar)

###############   PRINT OUTPUT FROM JOINTS   ###############

joints = [SJ, UJ, CJ, RJ, TJ]
names =  ['Spherical', 'Universal', 'Cylindrical', 'Revolute', 'Translational']

t = 0

for joint, name in zip(joints, names):
    
    N = 30
    title = '\n' + '#'*N + f'{name}'.center(20) + '#'*N + '\n'
    print(title)
    
    print(f'Phi:\n\n {joint.val(t)}')
    print(f'\nNu:\n\n {joint.vel_rhs(t)}')
    print(f'\nGamma:\n\n {joint.accel_rhs(t)}')
    
    #partial derivatives    
    dphi_dpi, dphi_dpj = joint.partial_p()
    print(f'\nPhi_pi:\n\n {dphi_dpi}')
    print(f'\nPhi_pj:\n\n {dphi_dpj}')
    
    dphi_dri, dphi_drj = joint.partial_r()
    print(f'\nPhi_ri:\n\n {dphi_dri}')
    print(f'\nPhi_rj:\n\n {dphi_drj}')
    
