"""
This will run simscript-gen.py one or more times according to maximum noise n_bar from: 
    -1 (no noise), 0 (vacuum), and increments of inc_n_bar up to max_n_bar

to run this:
    $ python noise_inc_sim.py ndata m ns r max_n_bar inc_n_bar sv

- ndata: number of data per simulation
- m: number of modes
- ns: output pattern. e.g: "000000"
- r: input squeezing parameter
- max_n_bar: maximum noise (avg photons in the noise modes)
- inc_n_bar: increment of noise
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
max_n_bar = float(sys.argv[5])
inc_n_bar = float(sys.argv[6])
sv = sys.argv[7]
t_noi = 0.5

n_bars = np.arange(0, max_n_bar+inc_n_bar, inc_n_bar) #array of noise n_bar
nb = []
nb.append(-1)
for i in range(n_bars.shape[0]):
    nb.append(n_bars[i])
n_bars = np.array(nb)

ts = timestamp() #one timestamp for all data from this simulation

for j in range(n_bars.shape[0]):
    bs_arr = sqbs_arr(m)
    n_bar = n_bars[j]
    com = (s+" "+ndata+" "+str(ns)+" "+str(r)+" "+bs_arr+" "+
           ts+" "+str(sv)+" "+str(n_bar)+" "+str(t_noi))
    print("Running command: "+com)
    os.system(com)



