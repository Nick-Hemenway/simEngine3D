import numpy as np
from simengine3D.utility import column

class BDF_Solver():
    
    def __init__(self, order = 1, h = 0.01):
        
        self.h = h
        self.order = order
        self.set_order(order)
        
    def set_order(self, order):
        
        self.order = order
        
        beta_arr = [1, 2/3, 6/11, 12/25, 60/137, 60/147]
        self.beta_0 = beta_arr[order - 1]
        
        #note that earlier indices corresponds to most recent points
        alpha_lst = [[-1],
                     [-4/3, 1/3], 
                     [-18/11, 9/11, -2/11], 
                     [-48/25, 36/25, -16/25, 3/25],
                     [-300/137, 300/137, -200/137, 75/137, -12/137],
                     [-360/147, 450/147, -400/147, 225/147, -72/147, 10/147]]
        
        #set the alpha array and reverse the order so that earlier indices
        #correspond to the oldest points
        alpha_coeff = alpha_lst[order-1]
        alpha_coeff.reverse()

        self.alpha_arr = column( alpha_coeff )
        
        
    def C_pos(self, pos_history, vel_history):
        
        C = - pos_history @ self.alpha_arr + self.beta_0*self.h*self.C_vel(vel_history)
        
        return column( C )
    
    def C_vel(self, vel_history):
        
        """vel_history is an array where the rows correspond to a variable and the 
        columns correspond to a time step. Earlier (lower index) columns correspond
        to older data"""
        
        C = - vel_history @ self.alpha_arr
        
        return column( C )
    
    def pos_step(self, accel, pos_history, vel_history):
        
        pos_n = self.C_pos(pos_history, vel_history) + self.beta_0**2* self.h**2 * accel
        
        return pos_n
        
    def vel_step(self, accel, vel_history):
        
        vel_n = self.C_vel(vel_history) + self.beta_0*self.h * accel
        
        return vel_n
        
        
def main():
    
    h = 0.1
    
    pos_h = np.array([[1,2],[1,2]])
    vel_h = pos_h*4
    
    bdf = BDF_Solver(order = 2, h = h)
    
    print(bdf.C_vel(vel_h))
    Cv = (4/3)*vel_h[:,1] - (1/3)*vel_h[:,0]
    print(Cv)
    
    print(bdf.C_pos(pos_h, vel_h))
    Cp = (4/3)*pos_h[:,1] - (1/3)*pos_h[:,0] + (8/9)*h*vel_h[:,1] - (2/9)*h*vel_h[:,0]
    print(Cp)
    
if __name__ == '__main__': main()