import sys
import pathlib as pl
root = pl.Path('../')
sys.path.append(str(root))

import simEngine3D as sim

r1 = [1,0,0]
r1_dot = [0,0,0]
p1 = [1, 0, 0, 0]
p1_dot = [0,0,0,0]

r2 = [1,0,0]
r2_dot = [0,0,0]
p2 = [1, 0, 0, 0]
p2_dot = [0,0,0,0]


body1 = sim.RigidBody(r1, r1_dot, p1, p1_dot)
body2 = sim.RigidBody(r2, r2_dot, p2, p2_dot)

a1_bar = [1,0,0]
a2_bar = [2,0,0]

dp1 = sim.GconDP1(body1, a1_bar, body2, a2_bar)

print(dp1.accel_rhs(0))

