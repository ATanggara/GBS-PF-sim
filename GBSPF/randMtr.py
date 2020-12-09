#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 12:47:24 2020

@author: andrewtanggara

Functions to generate random unitary matrix and random symplectic matrix
"""

import numpy as np

def randomUnitary(m, re=False):
    """
    returns a random mxm unitary matrix
    m is number of modes
    re true for orthogonal matrix (real entries), false for otherwise
    
    Based on: http://www.qetlab.com/RandomUnitary
    """
    U = np.random.randn(m,m)
    if not re:
        U = U + 1j*np.random.randn(m,m)
    q,r = np.linalg.qr(U)
    diag = np.diag(np.sign(np.diag(r)))
    U = q@diag
    return np.matrix(U)

def unitaryToSymplectic(U):
    u = U
    x = np.real(u)
    y = np.imag(-u)
    S = np.block([[x,y],
                  [-y,x]])
    return S

def randomSymplectic(m):
    """
    returns a random 2mx2m symplectic matrix
    m is number of modes
    """
    u = randomUnitary(m)
    return unitaryToSymplectic(u)

def transf_mtr(m):
    T = np.zeros((2*m,2*m))
    for i in range(m):
        T[i*2, i] = 1
        T[i*2+1, i+m] = 1
    return T

