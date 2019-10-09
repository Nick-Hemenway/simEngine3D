import numpy as np
from utility import tilde, column

class Orientation():
    
    def __init__(self):
        
        pass

class RigidBody():
    
    def __init__(self, r = None, p = None, r_dot = None, p_dot = None, idx = None, name = None):
        
        self.idx = idx
        
        if name == None:
            name = str(idx)
        
        #set default values to zero if none are provided
        if r is None:
            r = [0,0,0]
        if p is None:
            p = [0,0,0,0]
        if r_dot is None:
            r_dot = [0,0,0]
        if p_dot is None:
            p_dot = [0,0,0,0]
        
        #set attributes ensuring they are all column vectors
        self.r =      column(r)
        self.r_dot =  column(r_dot)
        self.p =      column(p)
        self.p_dot =  column(p_dot)
        
    def set_position(self, r):
        self.r = column(r)
    
    def set_orientation(self, q):
        self.q = column(q)
        
    def set_vel(self, r_dot):
        self.r_dot = column(r_dot)
        
    def set_ang_vel(self, p_dot):
        self.p_dot = column(p_dot)
    
    @property
    def A(self):
        
        e0 = self.p[0]
        e = self.p[1::]
        e_tilde = tilde(e)
        
        A_mat = (2*e0**2 - 1)*np.eye(3) + 2*(e @ e.T + e0*e_tilde)
        
        return A_mat
    

def main():
    
    r = [1,0,0]
    r_dot = [0,0,0]
    p = [1, 0, 0, 0]
    p_dot = [0,0,0]
    
    b1 = RigidBody(r, r_dot, p, p_dot)
    print(b1.A)


if __name__ == '__main__': main()
       
        
        
        
        
        