#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 12:17:45 2020

@author: andrewtanggara
"""

from GBSPF.PFv1 import *
from GBSPF.util import *

m = int(sys.argv[1])
k = int(sys.argv[2])
r = float(sys.argv[3])
sv = sys.argv[4]
xi = float(sys.argv[5])
gam = float(sys.argv[6])


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
        m,k,beta,vs,cir,stepN=200,maxN=20000)


#### CALCULATE HAFNIAN PROBS for all collision free output states
print("\n** Calculating Hafnian Probability...")
Pn_hafs = []
for ns in nss_mk:
    Pn_haf = prob_haf_gen(rs, ns, bs_arr, t, n_bar, t_noi)
    Pn_hafs.append(Pn_haf)
Pn_hafs = np.array(Pn_hafs)
print("** Probability from hafnian: "+str(Pn_hafs))


