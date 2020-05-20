"""
This will run simscript-gen.py one or more times according to maxmodes
    modes will be incremented by 2 starting from 2 modes.

to run this:
    $ python modes_inc_sim.py ndata modes ns r_max r_inc

- ndata: numebr of data per simulation
- modes: number of modes
- ns: output pattern
- r_max: maximum squeezing param r (in x.xx format. e.g: "4.0")
- r_inc: increment of squeezing r starting from 0 (e.g: "0.25")

@author: andrewtanggara
"""

import os
import sys
import numpy as np

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

modes = int(sys.argv[2])
ns= sys.argv[3]
r_max = float(sys.argv[4])
r_inc = float(sys.argv[5])
ndata = sys.argv[1]

rinc = np.arange(0,r_max+r_inc,r_inc)

for r in rinc:
    bs_arr = sqbs_arr(modes)
    com = s+ndata+" "+str(ns)+" "+str(r)+" "+bs_arr
    print("Running command: "+com)
    os.system(com)
    




