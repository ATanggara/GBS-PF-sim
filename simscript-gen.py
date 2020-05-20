"""
To run: $ python simscript-gen.py x y r b
    - x is number of data
    - y is output pattern '0011' or '11' or '1001', etc.
    - r is squeezing param grater than or equal to 0
    - b is beamsplitter arrangement 132 means means mix mode 1&2, 3&4, then 2&3
"""

import numpy as np
import hafnian
import datetime
import sys

from GBSPF.opers import *
from GBSPF.PF import *


#### Define GBS output pattern and interferometer params

ns = []
for i in range(len(sys.argv[2])):
    ns.append(int(sys.argv[2][i]))
ns = np.array(ns)

#squeezing param
r = float(sys.argv[3])
rs = []
for i in range(len(sys.argv[2])):
    if (i%2)==0:
        rs.append(-r)
    else:
        rs.append(r)
rs = np.array(rs)

t = 0.5 #bs transmissivity

#beamsplitter arrangement
bs_arr = []
for i in range(len(sys.argv[4])):
    bs_arr.append(int(sys.argv[4][i]))
bs_arr = np.array(bs_arr)


#### CALCULATE HAFNIAN PROB

Pn_haf = prob_haf(rs, ns, bs_arr, t)
print("Probability from hafnian: "+str(Pn_haf))


#### CALCULATE PF PROB

#compute first and second moments
mu = np.zeros(ns.shape[0]*2)
s = msqueezer(rs).T@msqueezer(rs)
for i in range(bs_arr.shape[0]):
    s = beamsplitter(s,t,bs_arr[i])
    
#create data
ndata = int(sys.argv[1])
print("\n========================\nCreating "+str(ndata)+" homodyne data...\n")
qp = []
for i in range(ndata):
    #rand phase shift
    ts = np.ones(ns.shape[0])*np.random.rand(1)
    pi = np.ones(ns.shape[0])*np.pi
    s_rand = mphasor(np.multiply(2*pi,ts)).T@s@mphasor(np.multiply(2*pi,ts))
    qp.append(np.random.multivariate_normal(mu, s_rand*0.5))
qp = np.array(qp)

#calculate prob
print("Calculating probability using PF...")
pfss, pfprods, pfsums, pfavgs, errs, errsums = pf_avg_prod_sum(qp, ns)


#### save result

date = str(datetime.date.today())
hr = str(datetime.datetime.now().hour)
mint = str(datetime.datetime.now().minute)
time = date+"-"+hr+""+mint #timestamp file
reldir = "gen_res/" #relative save directory
pattern = "("
for i in range(ns.shape[0]):
    if i!=0:
        pattern = pattern + ","
    pattern = pattern + str(ns[i])
pattern = pattern + ")"

ress = [pfss, pfprods, pfsums, pfavgs, errs, errsums]
resn = ["pfss", "pfprods", "pfsums", "pfavgs", "errs", "errsums"]
for i in range(len(ress)):
    np.save(reldir+resn[i]+pattern+"_"+str(ndata)+"_"+str(r)+" - "+str(time), ress[i])

