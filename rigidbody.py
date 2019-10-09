import numpy as np
from utility import tilde

class RigidBody():
    
    def __init__(self, r, r_dot, p, p_dot):
        
        self.r = np.atleast_2d(r).reshape(-1,1)
        self.r_dot = np.atleast_2d(r_dot).reshape(-1,1)
        self.p = np.atleast_2d(p).reshape(-1,1)
        self.p_dot = np.atleast_2d(p_dot).reshape(-1,1)
    
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
       
        
        
        
        
        