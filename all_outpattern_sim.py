"""
Simulate for all possible output patterns

This will run simscript-gen-v1.py once for a specific setting 
according to number of modes m and number of input squeezed state k.

to run this:
    $ python all_outpattern_sim.py m k r sv xi gam maxN stepN

- m:  number of modes
- k: number of input squeezed st (k<m must be true)
- r: input squeezing parameter
- sv: save all if "1", if "0" only save pfavgs and errs
- xi: value in (0,1) - lower-bound for PF-GBS approximation error
- gam: value in (0,1) - upper-bound for probability measure of k with error larger than xi
- maxN: maximum number of homodyne data
- stepN: number of homodyne data addition each step (while the sum of stdev is
            larger than the bound)

@author: andrewtanggara
"""


import os
import sys
from GBSPF.util import *

s = "python3 simscript-gen-v1.py"


m = int(sys.argv[1])
k = int(sys.argv[2])
r = sys.argv[3]
sv = sys.argv[4]
xi = float(sys.argv[5])
gam = float(sys.argv[6])
maxN = int(sys.argv[7])
stepN = int(sys.argv[8])


ts = timestamp() #one timestamp for all data from this simulation
n_bar = -1 #no noise
t_noi = 1 #no noise


bs_arr = sqbs_arr(m) #generate interferometer
com = (s+" "+str(m)+" "+str(k)+" "+str(r)+" "+bs_arr+" "+
       ts+" "+str(sv)+" "+str(n_bar)+" "+str(t_noi)+" "+
       str(xi)+" "+str(gam)+" "+str(maxN)+" "+str(stepN))
print("Running command: "+com)
os.system(com)

