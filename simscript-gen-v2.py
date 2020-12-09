"""
To run: $ python simscript-gen-v2.py m k r ns ts maxN nInterfs nbar t_noi
    - m number of modes
    - k number of squeezing inputs
    - r is squeezing param grater than or equal to 0
    - ns is output pattern separated by comma. e.g: "1,0,5,12,3"
    - maxN is number of Homodyne data
    - nInterfs is number of random interferometers
    - nbar is noise amount (equal to all modes), -1 being no noise, 0 being vacuum
    - t_noi is transmissivity of noise BS
    - sa is save mode {0,1,2}

@author: andrewtanggara
"""

import numpy as np
import datetime
import sys

from GBSPF.opers import *
from GBSPF.opers_haf_v2 import prob_haf_gen_interf
from GBSPF.PFv2 import *
from GBSPF.util import *
from GBSPF.randMtr import randomUnitary


m = int(sys.argv[1])
k = int(sys.argv[2])
r = float(sys.argv[3])
ns = str(sys.argv[4])
maxN = int(sys.argv[5])
nInterfs = int(sys.argv[6])
n_bar = float(sys.argv[7]) #noise amount
t_noi = float(sys.argv[8]) #transmissivity of noise beamsplitter
sa = int(sys.argv[9])
#b = str(sys.argv[4])
#xi = float(sys.argv[9])
#gam = float(sys.argv[10])
#stepN = int(sys.argv[12])

n_bars = np.ones(nInterfs)*n_bar
t_nois = np.ones(nInterfs)*t_noi

#import random
#random.seed(0)

#squeezing params
rs = []
for i in range(k): #only first k modes are squeezed
    if (i%2)==0:
        rs.append(-r)
    else:
        rs.append(r)
for i in range(m-k):
    rs.append(0)
rs = np.array(rs)


### generate array of random unitary interferometers
interfs = []
for i in range(nInterfs):    
    interfs.append(randomUnitary(m, re=False))

#calculate prob
print("Calculating probabilities using PF...")
pfsss, pfprodss, pfsumss, pfavgss, errss, errsumss = pf_avg_prod_sum_rand_interferometer(
        m,k,ns,interfs,rs,n_bars,t_nois,maxN=maxN)


#### CALCULATE HAFNIAN PROBS for all collision free output states
print("\n** Calculating Hafnian Probabilities...")
Pn_hafs = []
for i in range(nInterfs):  
    Pn_haf = prob_haf_gen_interf(rs=rs, ns=ns, T=interfs[i], n_bar=n_bar, t_noi=t_noi)
    Pn_hafs.append(Pn_haf)
Pn_hafs = np.array(Pn_hafs)
print("Probabilities from hafnian: "+str(Pn_hafs))


#### save result
date = str(datetime.date.today())
hr = str(datetime.datetime.now().hour)
mint = str(datetime.datetime.now().minute)
time = date+"-"+hr+""+mint #timestamp file

reldir = "gen_res/" #relative save directory

inputs = "("+str(m)+","+str(k)+")" #number of total modes and SMSS
    
ress = [pfsss, pfprodss, pfsumss, pfavgss, errss, errsumss]
resn = ["pfsss", "pfprodss", "pfsumss", "pfavgss", "errss", "errsumss"]


if sa=="1": #save all data
    for i in range(len(ress)):
        np.save(reldir+resn[i]+inputs+"_"+str(k)+"_"+str(r)+"_"+str(ns)+
            "_"+str(maxN)+"_"+str(n_bar)+" - "+str(time), ress[i])
elif sa=="0": #only save pfavgs and errs
    sv = [3,4]
    for i in sv:
        np.save(reldir+resn[i]+inputs+"_"+str(k)+"_"+str(r)+"_"+str(ns)+
            "_"+str(maxN)+"_"+str(n_bar)+" - "+str(time), ress[i])
    np.save(reldir+"interfs"+inputs+"_"+str(k)+"_"+str(r)+"_"+str(ns)+
            "_"+str(maxN)+"_"+str(n_bar)+" - "+str(time), interfs)
    np.save(reldir+"hafs"+inputs+"_"+str(k)+"_"+str(r)+"_"+str(ns)+
            "_"+str(maxN)+"_"+str(n_bar)+" - "+str(time), Pn_hafs)
else: 
    sv = [3,4]
    for i in sv:
        np.save(reldir+resn[i]+inputs+"_"+str(k)+"_"+str(r)+"_"+str(ns)+
            "_"+str(maxN)+"_"+str(n_bar)+" - "+str(time), ress[i])
    np.save(reldir+"hafs"+inputs+"_"+str(k)+"_"+str(r)+"_"+str(ns)+
            "_"+str(maxN)+"_"+str(n_bar)+" - "+str(time), Pn_hafs)

