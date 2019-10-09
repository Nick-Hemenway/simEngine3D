#add folders to the path so python knows where to import self made modules from
import sys
import pathlib as pl
gcon_folder = pl.Path('./Gcons/')
sys.path.append(gcon_folder)

import numpy as np
import Gcons
import rigidbody as body

class System():
    
    def __init__(self):
        
        pass
    
    def create_body(self):
        
        pass
    
    def constraint_DP1(self):
        
        pass
    
    def constraint_DP2(self):
        
        pass
    
    
    def constraint_D(self):
        
        pass
    
    def constraint_CD(self):
        
        pass
    
    def jacobian(self):
        
        pass
    
    

def main():
    
    pass


if __name__ == '__main__': main()