import numpy as np

def tilde(v):
    
    v = v.flatten()
    x,y,z = v
    v_tilde = np.array([[0, -z,  y],
                        [z,  0, -x],
                        [-y, x,  0]])
    
    return v_tilde
