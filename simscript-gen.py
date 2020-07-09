"""
To run: $ python simscript-gen.py n_data out r b t s n
    - n_data is number of data
    - out is output pattern '0011' or '11' or '1001', etc.
    - r is squeezing param grater than or equal to 0
    - b is beamsplitter arrangement 1,3,2 means means mix mode 1&2, 3&4, then 2&3
    - t is specified save timestamp. "0" to not specify
    - s is save all if "1", if "0" only save pfavgs and errs
    - n is noise amount (equal to all modes), -1 being no noise, 0 being vacuum
    - t_noi is transmissivity of noise BS
"""

import numpy as np
import hafnian
import datetime
import sys

from GBSPF.opers import *
from GBSPF.opers_haf import prob_haf_gen
from GBSPF.PF import *


#### Define GBS output pattern and interferometer params

ns = []
for i in range(len(sys.argv[2])):
    ns.append(int(sys.argv[2][i]))
ns = np.array(ns) #output pattern

m = len(sys.argv[2]) #number of modes

n_bar = float(sys.argv[7]) #noise amount
t_noi = float(sys.argv[8]) #transmissivity of noise beamsplitter


#squeezing param
r = float(sys.argv[3])
rs = []
for i in range(m):
    if (i%2)==0:
        rs.append(-r)
    else:
        rs.append(r)
rs = np.array(rs)

t = 0.5 #bs transmissivity

#beamsplitter arrangement
bs_arr_str = str.split(sys.argv[4], sep=',')
bs_arr = []
for i in range(len(bs_arr_str)):
    bs_arr.append(int(bs_arr_str[i]))
bs_arr = np.array(bs_arr)


#### CALCULATE HAFNIAN PROB
Pn_haf = prob_haf_gen(rs, ns, bs_arr, t, n_bar, t_noi)

print("Probability from hafnian: "+str(Pn_haf))


####### CALCULATE PF PROB

#compute first and second moments
mu = np.zeros(ns.shape[0]*2)
s = msqueezer(rs).T@msqueezer(rs)
for i in range(bs_arr.shape[0]):
    s = beamsplitter(s,t,bs_arr[i])

#put noise into cov matrix
for i in range(1,m+1):
    s,_ = therm_vac_noise(st=s, fm=i, n_bar=n_bar, t=t_noi)
np.set_printoptions(precision=3, linewidth=200)
print("quad cov mtr:")
print(s)

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

time = sys.argv[5]
if sys.argv[5] == "0":
    date = str(datetime.date.today())
    hr = str(datetime.datetime.now().hour)
    mint = str(datetime.datetime.now().minute)
    time = date+"-"+hr+""+mint #timestamp file

reldir = "gen_res/" #relative save directory
pattern = "(" #output pattern
for i in range(ns.shape[0]):
    if i!=0:
        pattern = pattern + ","
    pattern = pattern + str(ns[i])
pattern = pattern + ")"

ress = [pfss, pfprods, pfsums, pfavgs, errs, errsums]
resn = ["pfss", "pfprods", "pfsums", "pfavgs", "errs", "errsums"]

sa = sys.argv[6]
if sa=="1":
    for i in range(len(ress)):
        np.save(reldir+resn[i]+pattern+"_"+str(ndata)+"_"+str(r)+"_"+str(n_bar)+"_"+str(t_noi)+" - "+str(time), ress[i])
else: #only save pfavgs and errs
    sv = [3,4]
    for i in sv:
        np.save(reldir+resn[i]+pattern+"_"+str(ndata)+"_"+str(r)+"_"+str(n_bar)+"_"+str(t_noi)+" - "+str(time), ress[i])



