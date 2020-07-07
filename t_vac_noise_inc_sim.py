"""
Simulate vacuum noise GBS with different noise intensity (quantified by noise BS transmissivity)

This will run simscript-gen.py one or more times according to
transmissivity t_noi of noise beamsplitter: 
    value of t_noi between 0 and 1: 1 no noise, 0 all noise
    for all values increment of t_inc deom 0 to 1 

to run this:
    $ python noise_inc_sim.py ndata m ns r t_inc sv

- ndata: number of data per simulation
- m: maximum number of modes
- ns: output pattern. e.g: "000000"
- r: input squeezing parameter
- t_inc: increment of noise BS transmissivity from 0 to 1
- sv: save all if "1", if "0" only save pfavgs and errs

@author: andrewtanggara
"""

import os
import sys
import numpy as np
from GBSPF.util import *


s = "python3 simscript-gen.py"

ndata = sys.argv[1]
m = int(sys.argv[2])
ns= sys.argv[3]
r = sys.argv[4]
t_inc = float(sys.argv[5])
sv = sys.argv[6]
n_bar = 0 #0 avg photon thermal st. (vacuum noise)


t_nois = np.arange(0, 1+t_inc, t_inc) #array of noise BS transmissivities

ts = timestamp() #one timestamp for all data from this simulation

### run simulations
for j in range(t_nois.shape[0]):
    bs_arr = sqbs_arr(m)
    com = (s+" "+ndata+" "+str(ns)+" "+str(r)+" "+bs_arr+" "+
           ts+" "+str(sv)+" "+str(n_bar)+" "+str(t_nois[j]))
    print("Running command: "+com)
    os.system(com)



