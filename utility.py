import numpy as np

def tilde(v):
    
    v = v.flatten()
    x,y,z = v
    v_tilde = np.array([[0, -z,  y],
                        [z,  0, -x],
                        [-y, x,  0]])
    
    return v_tilde

def column(v):
    """converts any array like into a numpy column vector"""
    return np.atleast_2d(v).reshape(-1,1)

def main():
    
    v_test = [1,2,3,4]
    v_col = column(v_test)
    print(v_col)
    
if __name__ == '__main__': main()