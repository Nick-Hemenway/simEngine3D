#add folders to the path so python knows where to import self made modules from
import sys
import pathlib as pl
gcon_folder = pl.Path('./Gcons')
sys.path.append(str(gcon_folder))

import numpy as np
from GconPrimitives import *
from rigidbody import *

class System():
    
    def __init__(self):
        
        pass
    
    

def main():
    
#    dp1 = GconPrimitives.
    pass


if __name__ == '__main__': main()