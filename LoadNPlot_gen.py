"""
Load a result saved in ./EPR_res/
and plot it
To run: $ python LoanNpPlot.py x t
    - x is number of data
    - t is timestamp
"""

import numpy as np
from matplotlib import pyplot as plt
import hafnian
import sys

from GBSPF.opers import *
from GBSPF.PF import *


# Define GBS output pattern and interferometer params
r = 0.25
rs = np.array([-r,r,-r,r]) #squeezing param
t = 0.5 #bs transmissivity
m = rs.shape[0] #number of modes
bs_arr = np.array([1,3,2,1,3])
ns = np.array([1,1,1,1])

# calculate hafnian prob
Pn_haf = prob_haf(rs,ns,bs_arr)


### load result

time = "2020-05-15-1118" #file timestamp
reldir = "gen_res/" #relative save directory

# create output pattern string
pattern = "("
for i in range(ns.shape[0]):
    if i!=0:
        pattern = pattern + ","
    pattern = pattern + str(ns[i])
pattern = pattern + ")"

ndata = 50000
resn = ["pfss", "pfprods", "pfsums", "pfavgs", "errs", "errsums"]

lress = []
for i in range(len(resn)):
    lress.append(np.load(reldir+resn[i]+pattern+"_"+str(ndata)+"_"+str(r)+" - "+str(time)+".npy"))

pfssl, pfprodsl, pfsumsl, pfavgsl, errsl, errsumsl = (lress[0], lress[1], 
                                                      lress[2], lress[3], 
                                                      lress[4], lress[5])


#### Plot

#print("p(1,1) from PF: "+str(pfavgsl[-1]))
#print("from hafnian: "+str(Pn_haf))
#plt.figure()
#plt.plot(pfavgsl, label="p(1,1) from PF")
#plt.plot(np.array([0,pfavgsl.shape[0]]), Pn_haf*np.array([1,1]), label="p(1,1) from hafnian")
#plt.plot(pfavgsl+errsl, label="up error")
#plt.plot(pfavgsl-errsl, label="lower error")
#plt.ylim(-0.1,0.2)
#plt.ylabel("PF avg.")
#plt.xlabel("# of homodyne data")
#plt.legend()
#plt.savefig(reldir+"prob_"+pattern+"_"+str(ndata)+" - "+str(time)+".png")
#plt.show()
#
#print("error: " + str(errsl[-1]))
#plt.figure()
#plt.plot(errsl)
#plt.ylabel("error")
#plt.xlabel("# of homodyne data")
#plt.ylim(0,0.15)
#plt.savefig(reldir+"err_"+pattern+"_"+str(ndata)+" - "+str(time)+".png")
#plt.show()