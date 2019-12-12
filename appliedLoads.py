import numpy as np
import rigidbody
from utility import column, normalize, tilde

"""
Note: all specific types of forces and torques return a "Force" or "Torque" object
"""

class Force():
    
    def __init__(self, f, loc):
        
        """
        f: a function that takes time as an input and returns a vector object
           representing the force vector
           
        loc: an array like marking the location of where the force is being 
             applied on the rigid body
        """

        self.f = f
        self.loc = column(loc)
        
    def in_global(self, t):
        
        F = self.f(t)
        return F.in_global()
    
    def in_local(self, t):
        
        #returns force in frame it was created in
        
        F = self.f(t)
        return F.in_local()
    
    def to_local(self, t, frame):
        
        #returns force in arbitary reference frame "frame"
        
        F = self.f(t)
        return F.to_local(frame)
    
    def resulting_torque(self, t, frame):
        
        #frame is rigid_bodies frame that the force is being applied to
        
        F = self.f(t) #vector object
        force = F.to_local(frame) #bring it into the bodies local frame
        
        T = tilde(self.loc) @ force #torque in local frame as np array
        
        #make torque into actual vector object
        T = rigidbody.Vector(T, frame)
        
        return T
    
class Torque():
    
    def __init__(self, f):
        
        """
        f: a function that takes time as an input and returns a vector object
           representing the torque vector
        """

        self.f = f
        
    def in_global(self, t):
        
        T = self.f(t)
        return T.in_global()
    
    def in_local(self, t):
        
        #returns torque in frame it was created in
        
        T = self.f(t)
        return T.in_local()
    
    def to_local(self, t, frame):
        
        #returns torque in arbitary reference frame "frame"
        
        T = self.f(t)
        return T.to_local(frame)
        
        
def const_force(vec, loc):
    
    #vec is a vector object from within rigidbody.py
    #loc is location of force application
    f = lambda t: vec
    
    return Force(f, loc)

def const_torque(vec):
    
    #vec is a vector object from within rigidbody.py
    f = lambda t: vec
    
    return Torque(f)
    
class TSDA():
    
    """Translational Spring-Damper-Actuator"""
    
    def __init__(self, body_i, sp_i_bar, body_j, sq_j_bar, k, l_0, c, h):
        
        """
        sp_i_bar: location of point on body_i at which point force acts
        sq_j_bar: location of point on body_j at which point force acts

        k: spring constant of actuator
        l_0: unstretched length of spring
        c: damping coefficient
        h: function that describes the effects of an actuator (hydraulic, electric, etc.)
           must be of the form h(lij, lij_dot, t)
        t: time         
        
        """
        
        self.body_i = body_i
        self.body_j = body_j
        
        self.sp_i_bar = sp_i_bar
        self.sq_j_bar = sq_j_bar
        self.P = rigidbody.Point(sp_i_bar, body_i)
        self.Q = rigidbody.Point(sq_j_bar, body_j)
        
        self.k = k
        self.l_0 = l_0
        self.c = c
        self.h = h
        
        
    def _calc_info(self):
        
        """compute the distance between points P and Q as well as unit vector
        and time rate of change of length between both points
        """
        
        #get position and velocity of P and Q
        rP = self.P.position
        rQ = self.Q.position
        
        vP = self.P.velocity
        vQ = self.Q.velocity
        
        #position level calcs
        dij = rQ-rP
        
        lij = np.linalg.norm(dij)
        
        eij = dij/lij #unit vector pointing from point P to point Q
        
        #velocity level calcs
        dij_dot = vQ - vP
        
        lij_dot = ((dij_dot.T @ dij)/(np.sqrt(dij.T @ dij)) ).item()
        
        return eij, lij, lij_dot
    
    def _force_info(self, t):
        
        eij, lij, lij_dot = self._calc_info()
        
        F_mag = self.k*(lij - self.l_0) + self.c*lij_dot + self.h(lij, lij_dot, t)
        
        #eij points from P to Q
        return F_mag, eij
    
    def Fi(self, t):
        
        #returns force object for force on body i
        F_mag, eij = self._force_info(t)
        eij_local = self.body_i.A.T @ eij
        
        F_vec = self.body_i.vector(F_mag*eij_local)
        
        return F_vec
    
    def Fj(self, t):
        
        #returns force object for force on body j
        F_mag, eij = self._force_info(t)
        eij_local = self.body_j.A.T @ -eij
        
        F_vec = self.body_j.vector(F_mag*eij_local)
        
        return F_vec
    
    def create_force_objs(self):
        
        #creates force objects by passing force methods in as functions to the
        #force obj constructor
        Fi = Force(self.Fi, loc = self.sp_i_bar)
        Fj = Force(self.Fj, loc = self.sq_j_bar)
        
        return Fi, Fj
        
class RSDA():
    
    """Translational Spring-Damper-Actuator"""
    
    def __init__(self, body_i, ai_bar, bi_bar, body_j, aj_bar, bj_bar, k, theta_0, c, h):
        
        """
        ai_bar: axis about which the torque is applied to body i (unit vector)
        bi_bar: perpendicular to ai_bar and used with bj_bar to define rotation
                angle theta_ij
        
        aj_bar: axis about which the torque is applied to body j (unit vector)
        bj_bar: perpendicular to aj_bar and used with bi_bar to define rotation
                angle theta_ij
                
        k: spring constant of actuator
        theta_0: zero_tension angle
        c: damping coefficient
        h: function that describes the effects of an actuator (hydraulic, electric, etc.)
           must be of the form h(theta, theta_dot, t)
        t: time         
        
        """
        
        self.body_i = body_i
        self.body_j = body_j
        
        self.ai_bar = normalize(column(ai_bar))
        self.bi_bar = normalize(column(bi_bar))

        self.aj_bar = normalize(column(aj_bar))
        self.bj_bar = normalize(column(bj_bar))

        self.k = k
        self.theta_0 = theta_0
        self.c = c
        self.h = h
        
        self.theta_old = 0
        self.n = 0 #number of full revolutions
        
    def set_num_revolutions(self, n):
        
        """Manually set the number of revolutions for the revolution counter
        
        This function can be useful if you want to set an offset in the 
        revolution count of the RSDA joint. One particularly important use case
        is if you want the initial angle at the start of the simulation to be 
        negative instead of positive. In this case the initial revolutions should
        be set to be n = -1
        
        """
        self.n = n
    
    def _calc_theta_ij(self):
        
        #convert needed vectors to global frame
        bi = self.body_i.A @ self.bi_bar
        bj = self.body_j.A @ self.bj_bar
        ai = self.body_i.A @ self.ai_bar
        
        #flatten arrays to be used with numpys dot and cross functions
        bi = bi.flatten()
        bj = bj.flatten()
        ai = ai.flatten()

        #see pages 327 - 329 of Haugs book for explanation        
        gi = np.cross(ai, bi)
        
        c = bi @ bj
        s = gi @ bj
        
        arc = np.arcsin(s)
        
        if (s>=0) and (c>=0):
            theta = arc
        elif (s>=0) and (c<0):
            theta = np.pi - arc
        elif (s<0) and (c<0):
            theta = np.pi - arc
        else:
            theta = 2*np.pi + arc
            
        diff = theta - self.theta_old
        self.theta_old = theta
        
        #increment revolution counter depending on which way the zero degree 
        #line is crossed
        if diff < -0.9*2*np.pi:
            self.n += 1
        if diff > 0.9*2*np.pi:
            self.n -= 1
            
        #NOTE: The revolution counter feature works, but the initial angle at 
        #the beginning of the simulation must be positive (theta_init >= 0)
        
        return theta + 2*self.n*np.pi
    
    def _calc_theta_ij_dot(self):
        
        #see eq. 9.2.62 on page 335 of Haug
        Ai = self.body_i.A
        Aj = self.body_j.A
        
        omega_i_bar = self.body_i.omega_bar
        omega_j_bar = self.body_j.omega_bar
        
        theta_dot = self.ai_bar.T @ (Ai.T @ Aj @ omega_j_bar - omega_i_bar)
        theta_dot = theta_dot.item() #extract scalar from np array
        
        return theta_dot
        
    def Ti(self, t):
        
        theta_ij = self._calc_theta_ij()
        theta_ij_dot = self._calc_theta_ij_dot()
        
        T_mag = (self.k*(theta_ij - self.theta_0) + self.c*theta_ij_dot
                 + self.h(theta_ij, theta_ij_dot, t))
    
        #Note that this function calculates the torques in the global frame first
        #and then converts them to each of the local frames
        
        ai = self.body_i.A @ self.ai_bar
        Ti = T_mag*ai
        Ti_bar = self.body_i.A.T @ Ti
        
        Ti_bar = rigidbody.Vector(Ti_bar, self.body_i)  #turn Ti_bar into vector object      
        
        return Ti_bar
    
    
    def Tj(self, t):
        
        Ti_bar = self.Ti(t)
        Ti = Ti_bar.in_global()
        Tj = -Ti
        
        Tj_bar = self.body_j.A.T @ Tj
        Tj_bar = rigidbody.Vector(Tj_bar, self.body_j)
        
        return Tj_bar
    
    def create_torque_objs(self):
        
        #creates force objects by passing force methods in as functions to the
        #force obj constructor
        Ti = Torque(self.Ti)
        Tj = Torque(self.Tj)
        
        return Ti, Tj
    
