import numpy as np
from simengine3D.utility import tilde, column, A_to_p

class ReferenceFrame():
    
    def __init__(self, p = None, p_dot = None, p_ddot = None):
        
        if p is None:
            p = [1,0,0,0] #A = I_3
        
        if p_dot is None:
            p_dot = [0,0,0,0]

        if p_ddot is None:
            p_ddot = [0,0,0,0]
            
        #set attributes ensuring they are all column vectors
        self.set_orientation(p)
        self.set_ang_vel(p_dot)
        self.set_ang_accel(p_ddot)
        
    def set_orientation(self, p):

        if len(p) != 4:
            raise ValueError('p must be of length 4')
        
        else:
            self.p = column(p)
            
    def set_orientation_from_A(self, A):
        
        p = A_to_p(A)
        self.set_orientation(p)
        
    def set_ang_vel(self, p_dot):
        
        if len(p_dot) != 4:
            raise ValueError('p_dot must be of length 4')
        
        else:
            self.p_dot = column(p_dot)
            
    def set_ang_accel(self, p_ddot):
        
        if len(p_ddot) != 4:
            raise ValueError('p_ddot must be of length 4')
        
        else:
            self.p_ddot = column(p_ddot)
        
    def set_ang_vel_from_omega(self, omega):
        
        #comes from slide 16 of lecture 7
        omega = column(omega)
        
        p_dot = 0.5 * self.E.T @ omega
        
        self.set_ang_vel(p_dot)
    
    def set_ang_vel_from_omega_bar(self, omega_bar):
        
        #comes from slide 16 of lecture 7
        omega_bar = column(omega_bar)
        
        p_dot = 0.5 * self.G.T @ omega_bar
        
        self.set_ang_vel(p_dot)
        
    @property
    def omega(self):
       
        omega_vec = 2*self.E @ self.p_dot
        return omega_vec

    @property
    def omega_bar(self):

        omega_bar_vec = 2* self.G @ self.p_dot
        return omega_bar_vec
    
    @property
    def omega_dot(self):
        
        #returns angular acceleration (alpha/omega_dot) in the global reference
        #frame
        omega_dot_vec = 2*self.E @ self.p_ddot
        return omega_dot_vec
        
    @property
    def A(self):
        
        e0 = self.p[0]
        e = self.p[1::]
        e_tilde = tilde(e)
        
        A_mat = (2*e0**2 - 1)*np.eye(3) + 2*(e @ e.T + e0*e_tilde)
        
        return A_mat
    
    @property
    def E(self):
        
        #comes from slide 12 of lecture 7
        e0 = self.p[0]
        e = self.p[1::]
        e_tilde = tilde(e)
        
        E_mat = np.zeros((3,4))
        E_mat[:,0] = -e.flatten()
        E_mat[:,1::] = e_tilde + e0*np.eye(3)
        
        return E_mat
    
    @property
    def G(self):
        
        #comes from slide 12 of lecture 7
        e0 = self.p[0]
        e = self.p[1::]
        e_tilde = tilde(e)
        
        G_mat = np.zeros((3,4))
        G_mat[:,0] = -e.flatten()
        G_mat[:,1::] = -e_tilde + e0*np.eye(3)
        
        return G_mat
    
    @property
    def G_dot(self):
        
        #comes from slide 12 of lecture 7
        e0_dot = self.p_dot[0]
        e_dot = self.p_dot[1::]
        e_dot_tilde = tilde(e_dot)
        
        G_dot_mat = np.zeros((3,4))
        G_dot_mat[:,0] = -e_dot.flatten()
        G_dot_mat[:,1::] = -e_dot_tilde + e0_dot*np.eye(3)
        
        return G_dot_mat
    
    def vector(self, local):
        
        return Vector(local, self)
    
class Vector():
    
    """
    Class for describing vectors in general
    
    All vectors have an associated reference frame that they are created in
    """
    
    def __init__(self, vec, frame):
        
        """
        vec: the vector represented in the reference frame "frame"
        frame: the local frame that the input vector is represented in (can be
               the global frame)
        """
        
        self.frame = frame
        self.vec = column(vec)
        
    def in_global(self):
        
        #returns the local vector in the global frame
        return self.frame.A @ self.vec
    
    def in_local(self):
        
        #returns the vector in the local frame that it was created
        return self.vec
    
    def to_local(self, frame):
        
        #frame is the frame to which you want the vector to be expressed in
        
        glob_vec = self.in_global() #convert vector object to local frame
        local_vec = frame.A.T @ glob_vec #convert global vector to frame of interest
        
        return local_vec
    
    def __mul__(self, scalar):
        
        return Vector(vec = self.vec*scalar, frame = self.frame)
    
    def __rmul__(self, scalar):
        
        return Vector(vec = self.vec*scalar, frame = self.frame)
    
class Point():
    
    """
    A point is different from a vector in that it has a body associated with it
    and not just a reference frame. The point describes a point on a rigid body
    and allows users to do things like determine the position, velocicity, and 
    acceleration of any point on a rigid body.
    """
    
    def __init__(self, vec, body):
        
        self.vec = column(vec)
        self.body = body
    
    @property    
    def position(self):
        
        pos = self.body.r + self.body.A @ self.vec
        
        return pos
    
    @property
    def velocity(self):
        
        vel = self.body.r_dot + tilde(self.body.omega) @ self.body.A @ self.vec
        
        return vel

    @property        
    def acceleration(self):
        
        omega_tilde = tilde(self.body.omega)
        omega_dot_tilde = tilde(self.body.omega_dot)
        
        t1 = self.body.r_ddot
        t2 = omega_tilde @ omega_tilde @ self.body.A @ self.vec
        t3 = omega_dot_tilde @ self.body.A @ self.vec
        
        accel = t1 + t2 + t3
        
        return accel
    
class RigidBody(ReferenceFrame):
    
    def __init__(self, m, J, r = None, p = None, r_dot = None, p_dot = None, r_ddot = None, p_ddot = None, idx = None, name = None):
        
        self.m = m
        self.J = np.diag(J)
        self.idx = idx
        
        if name == None:
            name = str(idx)
        
        #set default values to zero if none are provided
        if r is None:
            r = [0,0,0]
        
        if r_dot is None:
            r_dot = [0,0,0]
            
        if r_ddot is None:
            r_ddot = [0,0,0]

        #set attributes ensuring they are all column vectors
        self.r =      column(r)
        self.r_dot =  column(r_dot)
        self.r_ddot = column(r_ddot)
        
        super().__init__(p, p_dot, p_ddot)
        
        self.forces = []
        self.torques = []
        
    def set_position(self, r):
        self.r = column(r)
    
    def set_vel(self, r_dot):
        self.r_dot = column(r_dot)
        
    def set_accel(self, r_ddot):
        self.r_ddot = column(r_ddot)
        
    def create_point(self, vec):
        
        p = Point(vec, self)
        
        return p
    
    def add_force(self, force):
        
        """force_vec must be a vector object"""
        
        self.forces.append(force)
        
    def add_torque(self, torque):
        
        """torque must be a vector object"""
        
        self.torques.append(torque)
        
    

def main():
    
    r = [1,0,0]
    r_dot = [0,0,0]
    p = [1, 0, 0, 0]
    p_dot = [0,0,0,0]
    
    b1 = RigidBody(m = 1, J = [1,1,1], r = r, r_dot = r_dot, p = p, p_dot = p_dot)
#    print(b1.p_dot)
    
    b1.set_ang_vel([1,2,3,4])
#    print(b1.p_dot)
    
    print(b1.A)
    v = b1.vector([1,2,3])
    
    v2 = 2*v
    
    print(v2.to_global())
    print(v2.to_local())


if __name__ == '__main__': main()
       
        
        
        
        
        