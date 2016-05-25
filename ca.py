#!/usr/bin/python
import numpy as np
import os
import re
import math
import sys

"""
get a list of transition time points in sequential order
"""
def trans_time(X):
    t_t = []
    """increase 1D neighbor search radius nbr to filter out thermal fluctuations,
        assuming symmetric peaks, for unsymmetric case need left nbr and right nbr"""
    nbr = 10
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
N th transition of a time series X time unit is in 100 picosecond
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
        if t > T[i] and t < T[i+1]:
            wait = T[i+1] - t
            break
        elif t < T[i]:
            wait = T[i] - t
            break
        elif t == T[i]:
            wait = T[i+1] - t
            break
    return wait


"""
get mean persistence time of X
"""
def  tau_p(X):
    T = trans_time(X)
    observ_t = T[-1]
    sum = 0
    for i in range(0,len(T)):
        sum += wait_time(X, trans_time_n(X,i+1))
    taup = math.pow(float(sum), 2)/float(2*observ_t)
    return taup


"""
get mean exchange time of X following the (i+1)th transition in Y
need to deal with the cases whereas after the transition time of Y,
no transtion occurred in X, in which case we need to remove this outlier
"""
def  tau_ex(X, Y):
    TX = trans_time(X)
    TY = trans_time(Y)
    observ_t_Y = TY[-1]
    sum = 0
    for i in range(0,len(TY)-1):
        w1 = wait_time(X, trans_time_n(Y,i+1))
        w2 = wait_time(Y, trans_time_n(Y,i))
        sum += w1 * w2
    tauex = float(sum)/float(observ_t_Y)
    return tauex


def get_ca(X, Y):
    CA = (-1) * math.log(tau_ex(X, Y)/tau_p(X))
    return CA

X = np.loadtxt('dih-sample22', usecols = (0,), unpack = True)
Y = np.loadtxt('dih-sample29', usecols = (0,), unpack = True)
CA = get_ca(X, Y)
print 'CA is %.3f' % (CA)
