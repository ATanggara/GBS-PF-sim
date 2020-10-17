#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 16:59:11 2020

to run this:
    $ python trim_npy.py sim_details n

where
    - sim_details: file name after the data description without ".npy"
    - n: only take every n-th element in array

@author: andrewtanggara
"""

import numpy as np
import sys
import os

m = sys.argv[1]
k = sys.argv[2]
sim_details = sys.argv[3]
timestamp = sys.argv[4]
n = int(sys.argv[5])

reldir = "gen_res/"
resn = ["pfavgss", "errss"]
lress = []

for i in range(len(resn)):
    lress.append(np.load(
            reldir+resn[i]+"("+m+","+k+")"+
            "_"+sim_details+" - "+timestamp+".npy"))

for i in range(len(lress)):
    lress[i] = np.array(lress[i])[:,::n]
    np.save(reldir+str(n)+"th"+"_"+resn[i]+"("+m+","+k+")"+"_"+sim_details+
            " - "+timestamp, lress[i])

