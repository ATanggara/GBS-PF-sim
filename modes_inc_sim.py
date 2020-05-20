"""
This will run simscript-gen.py one or more times according to maxmodes
    modes will be incremented by 2 starting from 2 modes.

to run this:
    $ python modes_inc_sim.py ndata maxmodes ns r

- ndata: number of data per simulation
- maxmodes: maximum number of modes
- ns: output pattern
- r: squeezing parameter

@author: andrewtanggara
"""

import os
import sys

def sqbs_arr(m):
    """
    return str of Beamsplitter arrangement in "1,4,7,12,14" format
    """
    barr = ""
    for j in range(1,m):
        bs_arrj = ""
        if j%2==1:
            for k in range(1,m,2):
                if j!=1 or k!=1:
                    bs_arrj = bs_arrj + ","
                bs_arrj = bs_arrj + str(k)
        else:
            for k in range(2,m-1,2):
                if j!=1:
                    bs_arrj = bs_arrj + ","
                bs_arrj = bs_arrj + str(k)
        barr = barr + bs_arrj
    return barr

s = "python3 simscript-gen.py "

max_modes = int(sys.argv[2])
r = sys.argv[4]
ndata = sys.argv[1]
ns = sys.argv[3]

for i in range(2,max_modes+1,2):
    print("\n---------- "+str(i)+" modes ----------")
    nsi = ns[:i]
    print("Output: "+nsi)
    
    bs_arr = sqbs_arr(i)
        
    print("Beamsplitter arrangement: "+bs_arr)
    com = s+ndata+" "+str(nsi)+" "+r+" "+bs_arr
    print("Running command: "+com)
    os.system(com)



