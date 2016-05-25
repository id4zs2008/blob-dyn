#!/usr/bin/python
import numpy as np
import os
import re
import math
import sys


"""get a list of transition time points in sequential order"""
def trans_time(X):
    t_t = []
    """increase 1D neighbor search radius nbr to filter out thermal fluctuations,
        assuming symmetric peaks, for unsymmetric case need left nbr and right nbr"""
    nbr = 4
    for i in range(0 + nbr, len(X) - nbr):
        peak = 1
        for j in range(1,nbr+1):
            if X[i] < X[i - j] or X[i] < X[i + j]:
                peak = 0
        if peak == 1:
            t_t.append(i+1)
    return t_t


"""
transition time function T(X, n) to get the time of
N th transition of a time series X time unit is 100 picosecond
"""
def trans_time_n(X, n):
    T = trans_time(X)
    return trans_time(X)[n-1]


"""
waiting time function W(X, t) to get the time interval
from t until the next transition/peak of X
"""
def wait_time(X, t):
    if t < 0 :
        sys.exit("Error: time needs to be a positive number")
    wait = 0
    T = trans_time(X)
    for i in range(0,len(T)-1):
        if t >= T[i] and t <= T[i+1]:
            wait = T[i+1] - t
            break
        elif t <= T[i]:
            wait = T[i] - t
            break
    return wait



X = np.loadtxt('dih-sample', usecols = (0,), unpack = True)
T = trans_time(X)
for i in range(0,len(T)):
    print T[i]
#print "#####"
#wait = wait_time(X, 93)
#print wait
