import numpy as np
from matplotlib import pyplot as plt
from itertools import product
from scipy import special
import hafnian
import datetime
import sys

from GBSPF.opers import *
from GBSPF.PF import *


# load result
time = "2020-05-07-1346" #file timestamp
reldir = "EPR_res/" #relative save directory
pattern = "(1,1)"

resn = ["pfss", "pfprods", "pfsums", "pfavgs", "errs"]
lress = []
for i in range(len(resn)):
    lress.append(np.load(np.load(reldir+resn[i]+pattern+","+str(ndata)+" - "+str(time)+".npy")))


pfssl, pfprodsl, pfsumsl, pfavgsl, errsl = lress[0], lress[1], lress[2], lress[3], lress[4]


plt.figure()
plt.plot(pfavgs)
plt.plot(np.array([0,ndata]), Pn_haf*np.array([1,1]))
plt.plot(pfavgs+stdevs)
plt.plot(pfavgs-stdevs)
plt.ylim(-0.1,0.2)
plt.ylabel("PF avg.")
plt.xlabel("# of homodyne data")
plt.show()

plt.figure()
plt.plot(errs)
# plt.ylim(0,0.02)
plt.show()
