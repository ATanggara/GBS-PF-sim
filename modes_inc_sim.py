"""
This will run simscript-gen.py one or more times according to maxmodes
    modes will be incremented by 2 starting from 2 modes.

to run this:
    $ python modes_inc_sim.py ndata maxmodes s

- ndata: numebr of data per simulation
- maxmodes: maximum number of modes
- s: squeezing parameter

@author: andrewtanggara
"""

import os
import sys

s = "python3 simscript-gen.py "

max_modes = int(sys.argv[2])
r = sys.argv[3]
ndata = sys.argv[1]
ns = ""

for i in range(2,max_modes+1,2):
    print("\n---------- "+str(i)+" modes ----------")
    ns = ns + "00"
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
    print("Running command: "+"python3 simscript-gen.py "+ndata+" "+str(ns)+" "+r+" "+bs_arr)
    os.system(s+ndata+" "+str(ns)+" "+r+" "+bs_arr)




