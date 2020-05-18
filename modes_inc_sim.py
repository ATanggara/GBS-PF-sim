"""
Created on Mon May 18 11:11:19 2020

python modes_inc_sim.py ndata maxmodes

@author: andrewtanggara
"""

import os
import sys

s = "python3 simscript-gen.py "

max_modes = int(sys.argv[2])
r = "0.25"
ns = ""
ndata = sys.argv[1]

for i in range(2,max_modes+1,2):
    print("--"+str(i))
    ns = ns + "00"
    print(ns)
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
#        print(bs_arrj)
    print(bs_arr)
    print("running command: "+"python simscript-gen.py "+ndata+" "+str(ns)+" "+r+" "+bs_arr)
    os.system(s+ndata+" "+str(ns)+" "+r+" "+bs_arr)




