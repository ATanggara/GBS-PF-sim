"""
Plot result of multiple multi-mode simulations

Created on Tue May 19 11:59:15 2020

@author: andrewtanggara
"""

import numpy as np
from matplotlib import pyplot as plt
from itertools import product
from scipy import special
import hafnian
import datetime

from GBSPF.opers import *
from GBSPF.PF import *

def haf_pars(r,ns,t):
    """
    Create parameters for Hafnian Probability:
    - rs: squeezing parameters
    - nsh: output pattern
    - bs_arrl: BS arrangement
    """
    bs_arr = ""
    m = len(ns)
    for j in range(1,m):
        bs_arrj = ""
        if j%2==1:
            for k in range(1,m,2):
                bs_arrj = bs_arrj + str(k)
        else:
            for k in range(2,m-1,2):
                bs_arrj = bs_arrj + str(k)
        bs_arr = bs_arr + bs_arrj

    nsh = [] #create array from string
    for i in ns:
        nsh.append(int(i))
    nsh = np.array(nsh)
    
    rs = [] #create array of squeezing params
    for i in range(m):
        rs.append(-1 * (-1)**i * r)
    rs = np.array(rs)
    
    bs_arrl = [] #create array from string for BS arrangement
    for i in bs_arr:
        bs_arrl.append(int(i))
    bs_arrl = np.array(bs_arrl)

    return rs, nsh, bs_arrl

    
### put parameters here:
max_modes = 10
r = 0.25
ns = ""
ndata = 50000
t = 0.5 #bs transmissivity - same for all simulations so far


# write timestamps for each simulation

stamps = np.array(["2020-05-19-1022", "2020-05-19-1022", "2020-05-19-1022",
                  "2020-05-19-1022", "2020-05-19-1023"]) #10k r=0.25
stamps = np.array(["2020-05-19-1025", "2020-05-19-1026", "2020-05-19-1026",
                  "2020-05-19-1027", "2020-05-19-1027"]) #20k r=0.25
stamps = np.array(["2020-05-19-1223", "2020-05-19-1223", "2020-05-19-1224",
                  "2020-05-19-1224", "2020-05-19-1225"]) #20k r=0.25
stamps = np.array(["2020-05-19-1226", "2020-05-19-1227", "2020-05-19-1228",
                  "2020-05-19-1230", "2020-05-19-1232"]) #50k r=0.25


############ Load files ############
    
pfssls = []
pfprodsls = []
pfsumsls = []
pfavgsls = []
errsls = []
errsumsls = []

Pn_hafs = []
nss = []

for j in range(int(max_modes/2)):
    ns = ns + "00"
    nss.append(ns)
    
    time = stamps[j] #file timestamp
    reldir = "gen_res/" #relative save directory

    # create output pattern string
    pattern = "("
    for i in range(len(ns)):
        if i!=0:
            pattern = pattern + ","
        pattern = pattern + str(ns[i])
    pattern = pattern + ")"

    resn = ["pfss", "pfprods", "pfsums", "pfavgs", "errs", "errsums"]
    resn = ["pfavgs", "errs"]
    
    lress = []
    for i in range(len(resn)):
        lress.append(np.load(reldir+resn[i]+pattern+"_"+str(ndata)+"_"+str(r)+" - "+stamps[j]+".npy"))

    pfavgsl, errsl = (lress[0], lress[1])
    
#     pfssls.append(pfssl)
#     pfprodsls.append(pfprodsl)
#     pfsumsls.append(pfsumsl)
    pfavgsls.append(pfavgsl)
    errsls.append(errsl)
#     errsumsls.append(errsumsls)
    
    rsh, nsh, bs_arrh = haf_pars(r,ns,t)
    
    Pn_haf = prob_haf(rsh,nsh,bs_arrh,t)
    Pn_hafs.append(Pn_haf)



######### plot probabilities #########
    
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
         "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]

plt.figure(figsize=(12,7))
plt.title("prob for $r="+str(r)+"$ up to "+str(max_modes)+" modes")
for i in range(int(max_modes/2)):
    print("p"+str(nss[i])+" from PF: "+str(pfavgsls[i][-1]))
    print("from hafnian: "+str(Pn_hafs[i]))
    plt.plot(pfavgsls[i], label="prob PF "+str(nss[i]), color=colors[i])
    plt.plot(np.array([0,pfavgsls[i].shape[0]]), Pn_hafs[i]*np.array([1,1]), label="prob hafnian "+str(nss[i]),
             linestyle=':', color=colors[i])
#     plt.plot(pfavgsls[i]+errsls[i], label="up error")
#     plt.plot(pfavgsls[i]-errsls[i], label="lower error")
#plt.ylim(0.5, 1)
i = int(max_modes/4)
plt.ylim(pfavgsls[i][-1]-errsls[i][-1]*13, pfavgsls[i][-1]+errsls[i][-1]*13)
plt.ylabel("PF avg.")
plt.xlabel("# of homodyne data")
# plt.xscale("log")
# plt.yscale("log")
plt.legend()
plt.show()

##plot error
plt.figure(figsize=(12,7))
plt.title("true and average error for $r="+str(r)+"$ up to "+str(max_modes)+" modes")
for i in range(int(max_modes/2)):
    print("error "+nss[i]+": " + str(errsls[i][-1]))
    print("|hafnian p - PF p| "+str(nss[i][-1])+": "+
          str(np.abs(pfavgsls[i] - Pn_hafs[i]*np.ones(pfavgsls[i].shape[0]))[-1]))
    plt.plot(errsls[i], label="error "+str(nss[i]), color=colors[i])
#     plt.plot(errsls[i]*2, label="error$\\times2$ "+str(nss[i]))
    plt.plot(np.abs(pfavgsls[i] - Pn_hafs[i]*np.ones(pfavgsls[i].shape[0])), label="|hafnian p - PF p| "+str(nss[i]),
             linestyle=':', color=colors[i])
plt.legend()
plt.ylabel("error")
plt.xlabel("# of homodyne data")
# plt.xscale("log")
# plt.yscale("log")
plt.ylim(0,errsls[i][-1]*10)
plt.show()
