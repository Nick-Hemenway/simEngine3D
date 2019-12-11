import simEngine3D as sim
import appliedLoads as loads
import rigidbody

body_i = rigidbody.RigidBody(m=4, J=[1,1,1], r = [0,0,0], r_dot = [-1,0,0])
body_j = rigidbody.RigidBody(m=2, J=[1,1,1], r = [1,0,0], r_dot = [1,0,0])

sp_i_bar = [1,1,1]
sq_j_bar = [1,1,1]

k = 10
l_0 = 1.5 #set l_0 to 1 to see effects of only damping
c = 2
h = lambda lij, lij_dot, t: 0

F = loads.TSDA(body_i, sp_i_bar, body_j, sq_j_bar, k, l_0, c, h)

eij, lij, lij_dot = F._calc_info()

print(f'\nUnit vector pointing from P to Q:\n\n {eij}')
print(f'\nDistance between P and Q:\n\n {lij}')
print(f'\nTime rate of change of length between P and Q:\n\n {lij_dot}')

Fp, Fq = F(10)

print(f'\nForce acting on point P:\n\n {Fp}')
print(f'\nForce acting on point Q:\n\n {Fq}')

