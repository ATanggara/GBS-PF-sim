import numpy as np
from matplotlib import pyplot as plt
from itertools import product
from scipy import special
import hafnian
import datetime
import sys

from GBSPF.opers import *
from GBSPF.PF import *


# # Define GBS output pattern and interferometer params
# ns = np.array([1,1])
# rs = np.array([-0.25,0.25]) #squeezing param
# t = 0.5 #bs transmissivity
# ## calculate hafnian prob
# s = sq_haf(rs)@sq_haf(rs).T
# s = (1/2) * beamsplitter_haf(s, t, 1)
# B = np.diag(np.tanh(rs))
# Bs = np.array([[np.sqrt(t), -np.sqrt(1-t)],
#               [np.sqrt(1-t), np.sqrt(t)]])
# B = Bs.T@B@Bs
# s_Q = s + np.eye(s.shape[0])*(1/2)
# det_s_Q = np.linalg.det(s_Q)
# #only for the case n=(1,1), therefore B_S=B 
# haf = hafnian.hafnian(B)
# #probability
# Pn_haf = (1/np.sqrt(det_s_Q)) * haf**2


# load result
time = "2020-05-07-1346" #file timestamp
reldir = "EPR_res/" #relative save directory
pattern = "(1,1)"
ndata = 50000

resn = ["pfss", "pfprods", "pfsums", "pfavgs", "errs"]
lress = []
for i in range(len(resn)):
    lress.append(np.load(np.load(reldir+resn[i]+pattern+"_"+str(ndata)+" - "+str(time)+".npy")))

pfssl, pfprodsl, pfsumsl, pfavgsl, errsl = lress[0], lress[1], lress[2], lress[3], lress[4]


plt.figure()
plt.plot(pfavgsl)
# plt.plot(np.array([0,ndata]), Pn_haf*np.array([1,1]))
plt.plot(pfavgsl+errsl)
plt.plot(pfavgsl-errsl)
plt.ylim(-0.1,0.2)
plt.ylabel("PF avg.")
plt.xlabel("# of homodyne data")
plt.show()

plt.figure()
plt.plot(errsl)
# plt.ylim(0,0.02)
plt.show()
