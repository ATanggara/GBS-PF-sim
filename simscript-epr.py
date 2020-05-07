import numpy as np
from matplotlib import pyplot as plt
from itertools import product
from scipy import special
import hafnian
import datetime
import sys

from GBSPF.opers import *
from GBSPF.PF import *


#### Define GBS output pattern and interferometer params

ns = np.array([1,1])
rs = np.array([-0.25,0.25]) #squeezing param
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
print("P(0,0) = "+str(1/np.sqrt(det_s_Q)))



#### CALCULATE PF PROB

#create data - EPR
m = np.zeros(4)
r = 0.25
ndata = int(sys.argv[1])
print("Creating "+str(ndata)+" homodyne data...")

qp = []
for i in range(ndata):
    s = cov_randphase_epr(r)
    qp.append(np.random.multivariate_normal(m, s))
qp = np.array(qp)

#calculate prob
print("Calculating probability using PF...")
ns = np.array([1,1])
pfss, pfprods, pfsums, pfavgs, errs = pf_avg_prod_sum(qp, ns, prtf=int(ndata/10))

# save result
date = str(datetime.date.today())
hr = str(datetime.datetime.now().hour)
mint = str(datetime.datetime.now().minute)
time = date+"-"+hr+":"+mint #timestamp file
reldir = "EPR_res/" #relative save directory
pattern = "(1,1)"

ress = [pfss, pfprods, pfsums, pfavgs, errs]
resn = ["pfss", "pfprods", "pfsums", "pfavgs", "errs"]
for i in range(len(ress)):
    np.save(reldir+resn[i]+pattern+" - "+str(time), ress[i])

