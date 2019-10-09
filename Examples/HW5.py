###############   IMPORT MODULES   ###############
#add the parent folder to the search path for importing modules
import sys
sys.path.append('..')

import numpy as np
import simEngine3D as sim


#%% ###############   DEFINE FUNCTION TO PRINT OUTPUTS   ###############
def constraint_output(constraint, name):
    
    """This function gets passed a constraint object, calls its various methods
       and prints it to the console in a formatted manner"""
    
    value = constraint.val(t = 0)
    nu = constraint.vel_rhs(t=0)
    gamma = constraint.accel_rhs(t = 0)
    dphi_dpi, dphi_dpj = constraint.partial_p()
    dphi_dri, dphi_drj = constraint.partial_r()
    
    N = 50
    title = '\n' + f'{name}'.center(10).center(N, '#') + '\n'
    print(title)
    print(f'Value: {value}')
    print(f'Nu: {nu}')
    print(f'Gamma: {gamma}')
    print(f'dphi_dri: {dphi_dri}')
    print(f'dphi_drj: {dphi_drj}')
    print(f'dphi_dpi: {dphi_dpi}')
    print(f'dphi_dpj: {dphi_dpj}')


#%% ##############   CREATE A SYSTEM   ##############

sys1 = sim.System()

#%% ###############   CREATE BODY 1   ###############

r1 = [8,6,-3]
r1_dot = [7,8,9]


#to get a valid set of Euler parameters we can normalize an arbitary vector
temp = np.array([4,3,-5,1])
p1 = temp/np.linalg.norm(temp) 

#to be valid, the derivative p_dot must be orthogonal to p
#to achieve this, we can set the first three values to an arbitrary value
#and then determine what the last element of p_dot needs to be to make the 
#two vectors orthogonal. We set the last element to zero to begin with
p1_dot = np.array([-0.2,1.3,3.4,0])

#now we determine what the last element needs to be as:
p1_dot[-1] = -np.dot(p1, p1_dot)/p1[-1]
p1_dot = p1_dot/np.linalg.norm(p1_dot)

body1 = sys1.add_body(r1, p1, r1_dot, p1_dot)

#%% ###############   CREATE BODY 2   ###############

r2 = [-0.5,1.6,-6.3]
r2_dot = [11, 12, 13]

#to get a valid set of Euler parameters we can normalize an arbitary vector
temp = np.array([3.3, -4, 5.1, 6])
p2 = temp/np.linalg.norm(temp) 

#again determining a valid set of Euler parameter derivatives
p2_dot = np.array([0.6, -3.7, 5.1, 0])
p2_dot[-1] = -np.dot(p2, p2_dot)/p2[-1]
p2_dot = p2_dot/np.linalg.norm(p2_dot)

body2 = sys1.add_body(r2, p2, r2_dot, p2_dot)

#%% ###############   CREATE CONSTRAINT FUNC   ###############

f = lambda t: 1.2
f_prime = lambda t: 2.5
f_pprime = lambda t: 0.2

func = sys1.constraint_func(f, f_prime, f_pprime)

#%% ###############   CREATE DP1 CONSTRAINT   ###############

ai_bar = [-1.2, 1 ,0.3]
aj_bar = [1.2, 4.5, 3.1]

DP1 = sys1.constraint_DP1(body1, ai_bar, body2, aj_bar, constraint_func = func)

constraint_output(DP1, name = 'DP1')

#%% ###############   CREATE DP2 CONSTRAINT   ###############
    
ai_bar = [-1.2, 1 ,0.3]
sp_i_bar = [0.1, -0.3, 6.0]
sq_j_bar = [0.2, -1.0, 1.5]

DP2 = sys1.constraint_DP2(body1, ai_bar, sp_i_bar, body2, sq_j_bar, constraint_func = func)

constraint_output(DP2, name = 'DP2')

#%% ###############   CREATE D CONSTRAINT   ###############

sp_i_bar = [0.1, -0.3, 6.0]
sq_j_bar = [0.2, -1.0, 1.5]

D = sys1.constraint_D(body1, sp_i_bar, body2, sq_j_bar, constraint_func = func)

constraint_output(D, name = 'D')

#%% ###############   CREATE CD CONSTRAINT   ###############

c_vec = [0.3, 0.4, -6]
sp_i_bar = [0.1, -0.3, 6.0]
sq_j_bar = [0.2, -1.0, 1.5]

CD = sys1.constraint_CD(body1, sp_i_bar, body2, sq_j_bar, c_vec, constraint_func = func)

constraint_output(CD, name = 'CD')