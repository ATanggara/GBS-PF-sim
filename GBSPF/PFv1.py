"""
Pattern Function

@author: andrewtanggara
"""

import numpy as np
from scipy import special
import itertools
import hafnian
import datetime
import sys
import math

from GBSPF.opers import *
from GBSPF.opers_haf import prob_haf_gen
from GBSPF.PF import *

PI = np.pi

def wav(q,n):
    """
    Reg wave function
    """
    if n==0:
        return PI**(-1/4) * np.exp(-(q**2)/2)
    else:
        return (1/np.sqrt(2*(n-1)+2)) * (q*wav(q,n-1) - dwav(q,n-1))

def irwav(q,n):
    """
    Irreg wave function
    """
    if n==0:
        return PI**(3/4) * np.exp(-(q**2)/2) * special.erfi(q)
    else:
        return (1/np.sqrt(2*(n-1)+2)) * (q*irwav(q,n-1) - dirwav(q,n-1))
    
def dwav(q,n):
    """
    Derivative of reg wave functions
    """
    if n==0:
        return PI**(-1/4) * (-q) * np.exp(-(q**2)/2)
    else:
        wf = wav(q,n-1)
        return (1/np.sqrt(2*(n-1)+2)) * (
            (2*(n-1)+2)*wf + q*dwav(q,n-1) - (q**2)*wf
        ) 

def dirwav(q,n):
    """
    Derivative of irreg wave functions
    """
    if n==0:
        return PI**(3/4) * np.exp(-(q**2)/2) * ((2/np.sqrt(PI))*np.exp(q**2) - q*special.erfi(q))
    else:
        wf = irwav(q,n-1)
        return (1/np.sqrt(2*(n-1)+2)) * (
            (2*(n-1)+2)*wf + q*dirwav(q,n-1) - (q**2)*wf
        ) 

def pfwav(q,n):
    """
    Pattern Function
    """
    return dwav(q,n) * irwav(q,n) + wav(q,n) * dirwav(q,n)


##### PF calculation sim
    
def pf_stdev(N, n, Pfs, PfAvgs):
    """
    calculate s_N
    Ps is approximation
    prob is probability from GBS
    """
    assert (N == Pfs.shape[0]) and (N == PfAvgs.shape[0])
    summ = 0
    for i in range(N):
        summ = summ + (Pfs[i] - PfAvgs[i])**2
    return np.sqrt((1/(N**2-N)) * summ)


def pf_avg_prod_sum(qp, ns):
    """
    - qp is complete homodyne data for both quadratures of all modes
        qp dimension is (N,2m) for number of data N and m modes
    - ns is output pattern
    
    returned values:
    - pfss: size=(N,m), pattern function for each mode of each hom data
    - pfprods: size=(N,1), product of PFs for each mode for each hom data, j-th element is PF for the j-th data
    - pfsums: size=(N,1), j-th element sum until j-th data
    - pfavgs: size=(N,1), j-th element is average of PFs until j-th data
    - errs: size=(N,1), j-th element is stdev until j-th data
    - errsums: size(N,1), j-th element is sum for errs calculation
    """
    m = ns.shape[0]
    pfss = [] 
    pfprods = [] 
    pfsums = [] 
    pfavgs = []
    errs = []
    errsums = []
    for j in range(qp.shape[0]): #loop over homodyne data
        if (j%int((qp.shape[0]/10)) == 0):
            print(str(j)+" homodyne data")
        ## calculate pattern function
        pfprod = 1
        pfs = [] # size=m, pattern function for each mode of j-th hom data
        for k in range(m): #loop over modes
            pf = pfwav(qp[j,2*k], ns[k])
            pfs.append(pf)
            pfprod = pfprod * pfs[-1]
        pfss.append(pfs)
        pfprods.append(pfprod)
        ## calculate sum
        if len(pfsums) == 0:
            pfsums.append(pfprods[j])
        else:
            pfsums.append(pfsums[-1] + pfprods[j])
        ## calculate average
        pfavgs.append(pfsums[j]/(j+1))
        ## calculate error
        if j==0:
            errsums.append((pfprods[0] - pfavgs[0])**2)
            errs.append( np.sqrt(errsums[j]) )
        else:
            errsums.append(errsums[-1] + pfprods[-1]**2 - 
                           (j+1)*pfavgs[j]**2 + j*pfavgs[j-1]**2)
            errs.append( np.sqrt( errsums[j]/((j+1)**2-(j+1)) ) )
    return (np.array(pfss), np.array(pfprods), 
            np.array(pfsums), np.array(pfavgs), 
            np.array(errs), np.array(errsums))
    

def generate_hom_samples(N, rs, bs_arr, t, n_bar, t_noi):
    """
    N: number of hom data
    rs: squeezing of input
    bs_arr: description of interferometer
    t: transmissivity of all BS in interferometer
    n_bar: noise amount
    t_noi: transmissivity of noise BS
    
    outputs N many m-mode homodyne samples, an array of size (ns,m)
    """
    m = rs.shape[0]
    #compute first and second moments
    mu = np.zeros(m*2)
    s = msqueezer(rs).T@msqueezer(rs)
    for i in range(bs_arr.shape[0]):
        s = beamsplitter(s,t,bs_arr[i])
    
    #put noise into cov matrix
    for i in range(1,m+1):
        s,_ = therm_vac_noise(st=s, fm=i, n_bar=n_bar, t=t_noi)
    np.set_printoptions(precision=3, linewidth=200)
    
    #create data
#    print("\n========================\nCreating "+str(N)+" homodyne data...\n")
    qp = []
    for i in range(N):
        #rand phase shift
        ts = np.ones(m)*np.random.rand(1)
        pi = np.ones(m)*np.pi
        s_rand = mphasor(np.multiply(2*pi,ts)).T@s@mphasor(np.multiply(2*pi,ts))
        qp.append(np.random.multivariate_normal(mu, s_rand*0.5))
    qp = np.array(qp)
    return qp

def no_cols(m,k):
    """
    Generate all possible combinations of no-collision output with k photons
    in an m modes system.
    """
    nss_all = np.array(list(map(list, itertools.product([0, 1], repeat=m))))
    nss_tot = np.sum(nss_all,1)
    return nss_all[np.where(nss_tot==k)] #arrays of {0,1} of size m with k 1s

def prod_erf(vs):
    return np.prod(special.erf(vs/np.sqrt(2)))

def sd_coef(c_mk, gam):
    """
    calculate set of coefficients of the standard deviation sum
    """
    vs = np.zeros(c_mk)
    incr = c_mk/100
    while prod_erf(vs/np.sqrt(2)) <= (1-gam):
        vs = vs + (np.ones(c_mk) * incr)
    vs_erf = prod_erf(vs/np.sqrt(2))
    print("Had vs satisfied the condition with erf(vs/sqrt(2)) = "+str(vs_erf)+
          " for 1-gamma = "+str(1-gam)+".")
    print("Tightening stdev coefficients vs...")
    incr = incr/1000
    while prod_erf(vs/np.sqrt(2)) > (1-gam):
        vs = vs - (np.ones(c_mk) * incr)
    vs = vs + 2*incr
    vs_erf = prod_erf(vs/np.sqrt(2))
    print("erf(vs/sqrt(2)) = "+str(vs_erf))
    return vs

def bs_str_to_array(bs_arr_str):
    """
    Convert string of beamsplitter decription "1,2,3,2,1" to an 
    ndarray [1,2,3,2,1].
    """
    bs_arr_str = str.split(bs_arr_str, sep=',')
    bs_arr = []
    for i in range(len(bs_arr_str)):
        bs_arr.append(int(bs_arr_str[i]))
    return np.array(bs_arr)


def pf_avg_prod_sum_bound_error(m ,k, beta, vs, cir, stepN=200, maxN=50000, printmul=5):
    """
    - m: number of modes
    - k: number of input squeezed state
    - beta: PF-GBS TV distance error
    - vs: size = c_mk, array of coefficients for stdev
    - cir: 5-tuple describing the Homodyne GBS circuit 
            (input to generate_hom_samples)
    - stepN: increment of Homodyne data to reach the stdev sum bound
    - maxN: maximum Homodyne data
    - printmul: result printed every printmul*stepN Homodyne samples
    
    returned values:
    - pfsss: size=(c_mk,N,m), pattern function for each mode of each hom data
    - pfprodss: size=(c_mk,N,1), product of PFs for each mode for each hom data,
            j-th element is PF for the j-th data
    - pfsumss: size=(c_mk,N,1), j-th element sum until j-th data
    - pfavgss: size=(c_mk,N,1), j-th element is average of PFs until j-th data
    - errss: size=(c_mk,N,1), j-th element is stdev until j-th data
    - errsumss: size(c_mk,N,1), j-th element is sum for errs calculation
    
    for c_mk size of no-collision set of states, and
    N total number of Homodyne data
    """
    beta = 2*beta

    nss_mk = no_cols(m,k) #array of no-collision
    c_mk = nss_mk.shape[0] #size of no-collision set (i.e: nss_mk)
    
    pfsss = [] 
    pfprodss = [] 
    pfsumss = [] 
    pfavgss = []
    errss = []
    errsumss = []
    
    stdev_sum = np.sum(vs)
    N = 0
    incr = stepN #increment of Homodyne data per loop
    print("\n*** For beta = "+str(beta)+", start generating Homodyne samples...\n")
    while (stdev_sum > beta): #keep creating Homodyne samples while stdev sum exceeds bound
        if (N >= maxN) and (maxN>0):
            break
        if (N%(incr*printmul) == 0):
            print("\n** "+str(stdev_sum)+" > beta = "+str(beta)+"\n"+
                  "Creating "+str(incr)+" more Homodyne samples ...")
        qp = generate_hom_samples(incr, cir[0], cir[1], 
                                      cir[2], cir[3], cir[4])
        for l in range(c_mk): #loop over no-collision outcomes
            pfss = []
            pfprods = []
            pfsums = []
            pfavgs = []
            errs = []
            errsums = []
            try:
                pfss = pfsss[l]
                pfprods = pfprodss[l] 
                pfsums = pfsumss[l] 
                pfavgs = pfavgss[l]
                errs = errss[l]
                errsums = errsumss[l]
            except IndexError:
                pfsss.append([])
                pfprodss.append([]) 
                pfsumss.append([])
                pfavgss.append([])
                errss.append([])
                errsumss.append([])
                pfss = pfsss[l]
                pfprods = pfprodss[l] 
                pfsums = pfsumss[l] 
                pfavgs = pfavgss[l]
                errs = errss[l]
                errsums = errsumss[l]
            ns = nss_mk[l]
            
            for j in range(qp.shape[0]): #loop over homodyne data
                ## calculate pattern function
                pfprod = 1
                pfs = [] # size=m, pattern function for each mode of j-th hom data
                for k in range(m): #loop over modes
                    pf = pfwav(qp[j,2*k], ns[k])
                    pfs.append(pf)
                    pfprod = pfprod * pfs[-1]
                pfss.append(pfs)
                pfprods.append(pfprod)
                ## calculate sum
                if len(pfsums) == 0:
                    pfsums.append(pfprods[-1])
                else:
                    pfsums.append(pfsums[-1] + pfprods[-1])
                ## calculate average
                pfavgs.append(pfsums[-1]/(N+j+1))
                ## calculate error
                if (N+j)==0:
                    errsums.append((pfprods[0] - pfavgs[0])**2)
                    errs.append( np.sqrt(errsums[-1]) )
                elif (N+j)==1:
                    errsums.append(errsums[-1] + pfprods[-1]**2 - 
                                   (N+j+1)*pfavgs[N+j]**2 + (N+j)*pfavgs[N+j-1]**2)
                    errs.append( np.sqrt( errsums[N+j]/((N+j)*(N+j)) ) )
                else:
                    errsums.append(errsums[-1] + pfprods[-1]**2 - 
                                   (N+j+1)*pfavgs[N+j]**2 + (N+j)*pfavgs[N+j-1]**2)
#                    errs.append( np.sqrt( errsums[N+j]/((N+j+1)**2-(N+j+1)) ) )
                    errs.append( np.sqrt( errsums[N+j]/((N+j)*(N+j-1)) ) )
            ## end loop over Homodyne samples
        ## end loop over no-collision states
        N = N + incr
        stdevs = np.array(errss)
        pav = np.array(pfavgss)
        if (N%(incr*printmul) == 0):
            print("loaded prev "+str(len(errs))+" values, N="+str(N))
            print("last stdevs "+str(stdevs[:,-1]))
            print("last PF avgs "+str(pav[:,-1]))
        stdev_sum = np.sum(np.dot(stdevs[:,-1], vs))
    ## end while loop
    print("\n**** Achieved stdev sum "+str(stdev_sum)+" <= beta = "+str(beta))
    print("last stdevs "+str(stdevs[:,-1]))
    pav = np.array(pfavgss)
    print("last PF avgs "+str(pav[:,-1]))
    return (np.array(pfsss), np.array(pfprodss), 
            np.array(pfsumss), np.array(pfavgss), 
            np.array(errss), np.array(errsumss))


