#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:23:45 2020

@author: andrewtanggara
"""

import numpy as np
from GBSPF.opers import dirsum
from GBSPF.randMtr import randomUnitary
from GBSPF.randMtr import randomSymplectic
from GBSPF.randMtr import transf_mtr


np.set_printoptions(precision=2)

m = 4

s = randomSymplectic(m)
print(s)

I = np.eye(m)
zero = np.zeros((m,m))
omega = np.block([[zero, I],
                 [-I, zero]])
b = s@omega@s.T

print()
print(b)
#print(np.isclose(b,omega))




T = transf_mtr(int(m))

sT = T@s@T.T #transformed random symplectic matrix into q,p,q,p form
print("sT")
print(sT)
omega2 = np.array([[0,1],
                   [-1,0]])
om = omega2
for i in range(m-1):
    om = dirsum(om,omega2)
print(om)
so = sT@om@sT.T

print()
print(np.isclose(so,om))


f = np.arange(1,200000,1000)
f = np.divide(np.ones(f.shape[0]), f)
f = f[1:]
from matplotlib import pyplot as plt
plt.figure(1)
plt.plot(f)
plt.show()
print(f)






