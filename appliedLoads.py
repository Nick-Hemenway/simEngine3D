import numpy as np
import rigidbody
from utility import column, tilde, normalize
        
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
        Vpq = vP - vQ #relative velocity of P w.r.t. Q
        
        #project relative velocity along unit vector in direction of two points
        lij_dot = (Vpq.T @ eij)[0,0] 
        
        return eij, lij, lij_dot
    
    def __call__(self, t):
        
        eij, lij, lij_dot = self._calc_info()
        
        F_mag = self.k*(lij - self.l_0) + self.c*lij_dot + self.h(lij, lij_dot, t)
        
        #eij points from P to Q
        Fp = F_mag*eij
        Fq = -F_mag*eij
        
        return Fp, Fq

        
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
        
        #need to implement logic for revolutions greater than one
        self.n = 0 #number of full revolutions
    
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
        
    def __call__(self, t):
        
        theta_ij = self.calc_theta_ij()
        theta_ij_dot = self.calc_theta_ij_dot()
        
        T_mag = (self.k*(theta_ij - self.theta_0) + self.c*theta_ij_dot
                 + self.h(theta_ij, theta_ij_dot, t))
    
        #Note that this function calculates the torques in the global frame first
        #and then converts them to each of the local frames
        
        ai = self.body_i.A @ self.ai_bar
        Ti = T_mag*ai
        Tj = -Ti
        
        Ti_bar = self.body_i.A.T @ Ti
        Tj_bar = self.body_j.A.T @ Tj
        
        return Ti_bar, Tj_bar
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    