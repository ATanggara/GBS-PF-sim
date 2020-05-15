#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 00:07:31 2020

@author: andrewtanggara
"""

import numpy as np
from scipy import special

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
    - prtf print everytime n hom data is processed
    
    returned values:
    - pfss: size=(N,m), pattern function for each mode of each hom data
    - pfprods: size=(N,1), product of PFs for each mode for each hom data, j-th element is PF for the $
    - pfsums: size=(N,1), j-th element sum until j-th data
    - pfavgs: size=(N,1), j-th element is average of PFs until j-th data
    - stdevs: size=(N,1), j-th element is stdev until j-th data
    """
    m = ns.shape[0]
    pfss = [] 
    pfprods = [] 
    pfsums = [] 
    pfavgs = []
    errs = []
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
        sqsum = 0
        for l in range(j):
            sqsum = sqsum + (pfprods[l] - pfavgs[j])**2
        if j==0:
            errs.append( np.sqrt(sqsum) )
        else:
            errs.append( np.sqrt( sqsum/((j+1)**2-(j+1)) ) )
    return np.array(pfss), np.array(pfprods), np.array(pfsums), np.array(pfavgs), np.array(errs)


