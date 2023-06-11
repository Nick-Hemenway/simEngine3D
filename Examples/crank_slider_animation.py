import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.style.use('seaborn-darkgrid')
plt.rcParams['font.family'] = 'Times New Roman'

def update(i, fig, ax, params):
    
    t, py, pz, sy, sz = params

    ax.cla()
    
    L_crank = py[0]
    L_rod = sy[0] - L_crank
    
    x1 = [0, py[i]]
    y1 = [0, pz[i]]
    
    x2 = [py[i], sy[i]]
    y2 = [pz[i], sz[i]]
    
    xmax = L_crank + L_rod
    ax.set_xlim(-1.5*L_crank, 1.1*xmax)
    ax.set_ylim(-1.5*L_crank, 1.5*L_crank)
    
    ax.set_xlabel('$y$ [m]')
    ax.set_ylabel('$z$ [m]')
    ax.set_aspect('equal')
    # ax.grid()
    
    ax.plot(x1, y1, lw = 2)
    ax.plot(x2, y2, lw = 2)
    ax.plot(sy[i], sz[i], marker = 's', markersize = 10)
    
    fig.tight_layout()
    
def animate(fname):   
    
    t, px, py, pz, sx, sy, sz = np.loadtxt('output/Crank_Slider_Data.txt', unpack=True)
    
    fps = 1/(t[1] - t[0])
    
    params = (t, py, pz, sy, sz)
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # update(10, fig, ax)
    
    frames = range(len(py))
    
    ani = animation.FuncAnimation(fig, update, frames = frames, fargs = (fig, ax, params))
    ani.save(fname, writer = 'ffmpeg', fps = fps, dpi = 150)
    
if __name__ == '__main__': 
    
    animate()
    