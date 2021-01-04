#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate for all possible output patterns according to the specific setting 
according to parameters described below.

to run this:
    $ python simscript_accuracy.py m k r sv xi gam

- m: number of modes
- k: number of input squeezed st (k<m must be true)
- r: input squeezing parameter
- sa: save all if "1", if "0" only save pfavgs and errs
- xi: value in (0,1) - lower-bound for PF-GBS approximation error
- gam: value in (0,1) - upper-bound for probability measure of k with error larger than xi
- maxN: maximum number of homodyne data (0 for no limit)
- stepN: number of homodyne data generated before each error bound check

@author: andrewtanggara
"""

from GBSPF.PFv1 import *
from GBSPF.util import *

m = int(sys.argv[1])
k = int(sys.argv[2])
r = float(sys.argv[3])
sa = sys.argv[4]
xi = float(sys.argv[5])
gam = float(sys.argv[6])
maxN = int(sys.argv[7])
stepN = int(sys.argv[8])


### generate stdev coefficients & its bound
nss_mk = no_cols(m,k)
c_mk = nss_mk.shape[0]
vs = sd_coef(c_mk, gam) #coefficients in the stdev sum
vs_erf = prod_erf(vs)
beta = ((xi * c_mk)/2) * ((vs_erf-1+gam)/vs_erf)


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

### Generate circuit description
bs_arr = bs_str_to_array(sqbs_arr(m)) #auto-generate a "square" interferometer
t = 0.5
n_bar = -1 #no noise
t_noi = 1 #no noise
cir = (rs,bs_arr,t,n_bar,t_noi) #Homodyne GBS circuit description


# Calculate prob from PF
print("Calculating probabilities using PF...")
pfsss, pfprodss, pfsumss, pfavgss, errss, errsumss = pf_avg_prod_sum_bound_error(
        m,k,beta,vs,cir,stepN=stepN,maxN=maxN)


#### CALCULATE HAFNIAN PROBS for all collision free output states
print("\n** Calculating Hafnian Probability...")
Pn_hafs = []
for ns in nss_mk:
    Pn_haf = prob_haf_gen(rs, ns, bs_arr, t, n_bar, t_noi)
    Pn_hafs.append(Pn_haf)
Pn_hafs = np.array(Pn_hafs)
print("** Probability from hafnian: "+str(Pn_hafs))


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



