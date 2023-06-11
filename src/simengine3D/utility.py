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

def A_to_p(A):
    """converts a valid orientation matrix to the corresponding set of Euler
    parameters"""
    
    trace = A.trace()
    p = np.zeros(4)
    p[0] = np.sqrt((trace + 1)/4)
    
    temp = np.zeros(3)
    for i in range(3):
        temp[i] = (1 + 2*A[i,i] - trace)/4
    
    #case 1: e0 != 0
    if p[0] != 0:
        sign_arr = np.sign([A[2,1] - A[1,2], A[0,2] - A[2,0], A[1,0] - A[0,1]])
        p[1::] = sign_arr * np.sqrt(temp)
    
    #case 2: e0 == 0
    else:   
        e_vec = p[1::]
        
        idx = np.argmax(temp)
        val = np.sqrt(temp[idx])
        p[idx + 1] = val
        
        if idx == 0:
            e_vec[1] = (A[1,0] + A[0,1])/(4*val)
            e_vec[2] = (A[2,0] + A[0,2])/(4*val)
                
        elif idx == 1:
            e_vec[0] = (A[1,0] + A[0,1])/(4*val)
            e_vec[2] = (A[2,1] + A[1,2])/(4*val)
            
        else:
            e_vec[0] = (A[2,0] + A[2,0])/(4*val)
            e_vec[1] = (A[2,1] + A[1,2])/(4*val)

    return p

def normalize(vec):
    
    return vec/np.linalg.norm(vec)

def rotate_x(theta):
    
    """Return rotation matrix A for rotation about x-axis by theta"""
    
    theta = np.deg2rad(theta)
    c = np.cos(theta)
    s = np.sin(theta)
    
    A = np.array([[  1, 0,  0 ],
                   [ 0, c, -s ],
                   [ 0, s,  c ]])
    
    return A

def rotate_y(theta):
    
    """Return rotation matrix A for rotation about y-axis by theta"""
    
    theta = np.deg2rad(theta)
    c = np.cos(theta)
    s = np.sin(theta)
    
    A = np.array([[  c,  0, s ],
                   [ 0,  1, 0 ],
                   [-s,  0, c ]])
    
    return A
    
def rotate_z(theta):
    
    """Return rotation matrix A for rotation about z-axis by theta"""
    
    theta = np.deg2rad(theta)
    c = np.cos(theta)
    s = np.sin(theta)
    
    A = np.array([[ c, -s, 0 ],
                   [ s,  c, 0 ],
                   [ 0,  0, 1 ]])
    
    return A
    
def main():
    
    v_test = [1,2,3,4]
    v_col = column(v_test)
    print(v_col)
    
if __name__ == '__main__': main()