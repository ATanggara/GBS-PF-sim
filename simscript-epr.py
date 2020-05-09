import numpy as np
import hafnian
import datetime
import sys

from GBSPF.opers import *
from GBSPF.PF import *

"""
To run: $ python simscript-epr.py x y r
    - x is number of data
    - y is pattern '00' or '11'
    - r is squeezing param grater than or equal to 0
"""


#### Define GBS output pattern and interferometer params

r = float(sys.argv[3])
ns = np.array([int(sys.argv[2][0]), int(sys.argv[2][1])])
rs = np.array([-r,r]) #squeezing param
t = 0.5 #bs transmissivity


#### CALCULATE HAFNIAN PROB

s = sq_haf(rs)@sq_haf(rs).T
s = (1/2) * beamsplitter_haf(s, t, 1)
print("cov matrix:")
print(s)

B = np.diag(np.tanh(rs))
Bs = np.array([[np.sqrt(t), -np.sqrt(1-t)],
              [np.sqrt(1-t), np.sqrt(t)]])
B = Bs.T@B@Bs
print("B matrix:")
print(B)

s_Q = s + np.eye(s.shape[0])*(1/2)
det_s_Q = np.linalg.det(s_Q)
print("s_Q:")
print(s_Q)
print("det of s_Q: "+str(det_s_Q))

#only for the case n=(1,1), therefore B_S=B 
haf = hafnian.hafnian(B)
print("\nhafnian: "+str(haf))

#probability
Pn_haf = (1/np.sqrt(det_s_Q)) * haf**2
print("\nP(1,1) = "+str(Pn_haf))

# probability for n = (0,0), assuming that the Haf(B_s) = 1
P0_haf = 1/np.sqrt(det_s_Q)
print("\nP(0,0) = "+str(P0_haf))



#### CALCULATE PF PROB

#create data - EPR
m = np.zeros(4)
ndata = int(sys.argv[1])
print("\n========================\nCreating "+str(ndata)+" homodyne data...\n")

qp = []
for i in range(ndata):
    s = cov_randphase_epr(r)
    qp.append(np.random.multivariate_normal(m, s))
qp = np.array(qp)

#calculate prob
print("Calculating probability using PF...")
pfss, pfprods, pfsums, pfavgs, errs = pf_avg_prod_sum(qp, ns)

# save result
date = str(datetime.date.today())
hr = str(datetime.datetime.now().hour)
mint = str(datetime.datetime.now().minute)
time = date+"-"+hr+""+mint #timestamp file
reldir = "EPR_res/" #relative save directory
pattern = "("+str(ns[0])+","+str(ns[1])+")"

ress = [pfss, pfprods, pfsums, pfavgs, errs]
resn = ["pfss", "pfprods", "pfsums", "pfavgs", "errs"]
for i in range(len(ress)):
    np.save(reldir+resn[i]+pattern+"_"+str(ndata)+" - "+str(time), ress[i])

