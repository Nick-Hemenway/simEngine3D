import numpy as np
import matplotlib.pyplot as plt

x_ddot, y_ddot, z_ddot = np.loadtxt('Acceleration.txt', unpack = True)

dt = 1e-3
t_end = 10

g = 9.81

time = np.arange(0, t_end + dt/2, dt)
theta = (np.pi/4) * np.cos(2*time)
theta_ddot = -np.pi* np.cos(2*time)


#sum T = I*alpha

L = 2 #length to center of rod = 2 meters
rho = 7800 #kg/m^3
w = 0.05 #m

m = rho*2*L*w**2
I = (1/12)*m*((2*L)**2 + w**2) + m*L**2 #inertia about pin

Fz = m*z_ddot + m*g
Fy = m*y_ddot

Tg = -m*g*L*np.sin(theta)

T = I*theta_ddot - Tg

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(time, Fz)
ax.plot(time, Fy)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.plot(time, T)
