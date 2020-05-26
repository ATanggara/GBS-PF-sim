"""
This will run simscript-gen.py one or more times according to maxmodes
    modes will be incremented by 2 starting from 2 modes.

to run this:
    $ python modes_inc_sim.py ndata modes ns r_max r_inc

- ndata: number of data per simulation
- modes: number of modes
- ns: output pattern
- r_max: maximum squeezing param r (in x.xx format. e.g: "4.0")
- r_inc: increment of squeezing r starting from 0 (e.g: "0.25")

@author: andrewtanggara
"""

import os
import sys
import numpy as np
from GBSPF.util import *


s = "python3 simscript-gen.py"

modes = int(sys.argv[2])
ns= sys.argv[3]
r_max = float(sys.argv[4])
r_inc = float(sys.argv[5])
ndata = sys.argv[1]

rinc = np.arange(0,r_max+r_inc,r_inc)

for r in rinc:
    bs_arr = sqbs_arr(modes)
    com = s+" "+ndata+" "+str(ns)+" "+str(r)+" "+bs_arr+" "+timestamp()
    print("Running command: "+com)
    os.system(com)
    



