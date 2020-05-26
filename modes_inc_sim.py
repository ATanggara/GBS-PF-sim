"""
This will run simscript-gen.py one or more times according to maxmodes
    modes will be incremented by 2 starting from 2 modes.

to run this:
    $ python modes_inc_sim.py ndata maxmodes ns r s

- ndata: number of data per simulation
- maxmodes: maximum number of modes
- ns: (maximum) output pattern. e.g: for 6 modes "000000"
- r: squeezing parameter
- s: save all if "1", if "0" only save pfavgs and errs

@author: andrewtanggara
"""

import os
import sys
from GBSPF.util import *


s = "python3 simscript-gen.py"

max_modes = int(sys.argv[2])
r = sys.argv[4]
ndata = sys.argv[1]
ns = sys.argv[3]
time = timestamp() #fix a timestamp for all simulations

for i in range(2,max_modes+1,2):
    print("\n---------- "+str(i)+" modes ----------")
    nsi = ns[:i]
    print("Output: "+nsi)
    
    bs_arr = sqbs_arr(i)
        
    print("Beamsplitter arrangement: "+bs_arr)
    com = s+" "+ndata+" "+str(nsi)+" "+r+" "+bs_arr+" "+time+" "+sys.argv[5]
    print("Running command: "+com)
    os.system(com)



