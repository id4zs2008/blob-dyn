#!/usr/bin/python
import numpy as np
from sklearn.metrics import mutual_info_score
import MDAnalysis
import os
import re

"""
Calculate conditional activity as implemented in Milo Lin JACS 2016"""

"""transition time function T(X, i) to get the time of Nth transition 
of a time series X, unit: 100 ps"""
def trans_time(X):
    t_t = []
    """increase 1D neighbor search radius nbr to filter out thermal fluctuations,
        assuming symmetric peaks, for unsymmetric case need left nbr and right nbr"""
    nbr = 5
    for i in range(1, len(X) - nbr):
        peak = 1
        for j in range(1,nbr+1):
            if X[i] < X[i - j] or X[i] < X[i + j]:
                peak = 0
        if peak == 1:
            t_t.append(i+1)
    return t_t

X = np.loadtxt('dih-sample', usecols = (0,), unpack = True)
T = trans_time(X)
for i in range(0,len(T)):
    print T[i]
