"""
This will run simscript-gen.py one or more times according to maxmodes
    modes will be incremented by 2 starting from 2 modes.

to run this:
    $ python modes_inc_sim.py ndata maxmodes ns r

- ndata: numebr of data per simulation
- maxmodes: maximum number of modes
- ns: output pattern
- r: squeezing parameter

@author: andrewtanggara
"""

import os
import sys

s = "python3 simscript-gen.py "

max_modes = int(sys.argv[2])
r = sys.argv[4]
ndata = sys.argv[1]
ns = sys.argv[3]

for i in range(2,max_modes+1,2):
    print("\n---------- "+str(i)+" modes ----------")
    nsi = ns[:i]
    print("Output: "+ns)
    bs_arr = ""
    for j in range(1,i):
        bs_arrj = ""
        if j%2==1:
            for k in range(1,i,2):
                bs_arrj = bs_arrj + str(k)
        else:
            for k in range(2,i-1,2):
                bs_arrj = bs_arrj + str(k)
        bs_arr = bs_arr + bs_arrj
    print("Beamsplitter arrangement: "+bs_arr)
    com = s+ndata+" "+str(ns)+" "+r+" "+bs_arr
    print("Running command: "+com)
    os.system(com)




