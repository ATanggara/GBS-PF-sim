"""
To run: $ python simscript-gen-v1.py m k r b ts s nbar t_noi xi gam
    - m number of modes
    - k number of squeezing inputs
    - r is squeezing param grater than or equal to 0
    - b is beamsplitter arrangement 1,3,2 means means mix mode 1&2, 3&4, then 2&3
    - ts is specified save timestamp. "0" to not specify
    - s is save all if "1", if "0" only save pfavgs and errs
    - nbar is noise amount (equal to all modes), -1 being no noise, 0 being vacuum
    - t_noi is transmissivity of noise BS
    - xi: value in (0,1) - lower-bound for PF-GBS approximation error
    - gam: value in (0,1) - upper-bound for probability measure of k with error larger than xi

@author: andrewtanggara
"""

import numpy as np
import hafnian
import datetime
import sys

from GBSPF.opers import *
from GBSPF.opers_haf import prob_haf_gen
from GBSPF.PFv1 import *
from GBSPF.util import *


m = int(sys.argv[1])
k = int(sys.argv[2])
r = float(sys.argv[3])
b = str(sys.argv[4])
n_bar = float(sys.argv[7]) #noise amount
t_noi = float(sys.argv[8]) #transmissivity of noise beamsplitter
xi = float(sys.argv[9])
gam = float(sys.argv[10])


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

t = 0.5 #bs transmissivity

### generate stdev coefficients its bound
nss_mk = no_cols(m,k)
c_mk = nss_mk.shape[0]
vs = sd_coef(c_mk, gam) #coefficients in the stdev sum
vs_erf = prod_erf(vs)
beta = ((xi * c_mk)/2) * ((vs_erf-1+gam)/vs_erf)

### Generate circuit description
bs_arr = bs_str_to_array(b)
t = 0.5
cir = (rs,bs_arr,t,n_bar,t_noi) #Homodyne GBS circuit description


#calculate prob
print("Calculating probabilities using PF...")
pfss, pfprods, pfsums, pfavgs, errs, errsums = pf_avg_prod_sum_bound_error(
        m,k,beta,vs,cir,stepN=500,maxN=50000)


#### CALCULATE HAFNIAN PROBS for all collision free output states
Pn_hafs = []
for ns in nss_mk:
    Pn_haf = prob_haf_gen(rs, ns, bs_arr, t, n_bar, t_noi)
    Pn_hafs.append(Pn_haf)
Pn_hafs = np.array(Pn_hafs)
print("Probability from hafnian: "+str(Pn_hafs))


#### save result
time = sys.argv[5]
if sys.argv[5] == "0":
    date = str(datetime.date.today())
    hr = str(datetime.datetime.now().hour)
    mint = str(datetime.datetime.now().minute)
    time = date+"-"+hr+""+mint #timestamp file

reldir = "gen_res/" #relative save directory

inputs = "("+str(m)+","+str(k)+")" #number of total modes and SMSS
    
ress = [pfsss, pfprodss, pfsumss, pfavgss, errss, errsumss]
resn = ["pfsss", "pfprodss", "pfsumss", "pfavgss", "errss", "errsumss"]

sa = sys.argv[6]
if sa=="1":
    for i in range(len(ress)):
        np.save(reldir+resn[i]+inputs+"_"+str(ndata)+"_"+str(r)+"_"+
                str(n_bar)+"_"+str(t_noi)+" - "+str(time), ress[i])
else: #only save pfavgs and errs
    sv = [3,4]
    for i in sv:
        np.save(reldir+resn[i]+inputs+"_"+str(ndata)+"_"+str(r)+"_"+
                str(n_bar)+"_"+str(t_noi)+" - "+str(time), ress[i])



