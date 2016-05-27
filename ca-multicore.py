#!/usr/bin/python
import numpy as np
from sklearn.metrics import mutual_info_score
import MDAnalysis
import os
import re
import math
import sys
from itertools import combinations_with_replacement,permutations
from concurrent.futures import ProcessPoolExecutor, Future, wait


usecpus = 10#how many cores to use
frms_num = 10000
u = MDAnalysis.Universe('ini.pdb','allpdb.trr')
f = open('CA-out.txt', 'w')
b = np.zeros((352,352))
#for i in range(0,352):
#    for j in range(0,352):
#        b[i][j] = 100


def new_dihedral(p):
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


"""get a list of transition time points in sequential order"""
def trans_time(X):
    t_t = []
    """increase 1D neighbor search radius nbr to filter out thermal fluctuations,
        assuming symmetric peaks, for unsymmetric case need left nbr and right nbr"""
    nbr = 10
    for i in range(0 + nbr, len(X) - nbr):
        peak = 1
        for j in range(1, nbr+1):
            if X[i] < X[i - j] or X[i] < X[i + j]:
                peak = 0
                break
        if peak == 1:
            t_t.append(i+1)
    find_basin = t_time(X, t_t)
    return find_basin
def t_time(X, t_t):
    rg = 1
    k_k = []
    for i in range(0, len(X)):
        peak1 = min(t_t, key=lambda x:abs(x-i))
        peak2 = min(t_t, key=lambda x:abs(x-(i+rg)))
        if peak1 != peak2:
            k_k.append(i)
    return k_k


"""
transition time function T(X, n) to get the time of
N th transition of a time series X time unit is 100 picosecond
"""
def trans_time_n(X, n, tx):
    T = tx
    return T[n-1]


"""
waiting time function W(X, t) to get the time interval
from t until the next transition/peak of X
"""
def wait_time(X, t, tx):
    #if t < 0 :
    #    sys.exit("Error: time needs to be a positive number")
    wait = 0
    T = tx
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
        #elif t > T[-1]:
        #    wait = -100
    return wait


"""
get mean persistence time of X
"""
def  tau_p(X, tx):
    T = tx
    observ_t = T[-1]
    sum = 0
    for i in range(0,len(T)):
        sum += wait_time(X, trans_time_n(X,i+1,T), T)
    taup = math.pow(float(sum), 2)/float(2*observ_t)
    return taup


"""
get mean exchange time of X following the (i+1)th transition in Y
in the cases whereas after the transition time of Y,
no transtion occurred in X, wait time is assigned to 0
"""
def  tau_ex(X, Y, tx, ty):
    TX = tx
    TY = ty
    observ_t_Y = TY[-1]
    sum = 0
    for i in range(0,len(TY)-1):
        w1 = wait_time(X, trans_time_n(Y,i+1,TY), TX)
        w2 = wait_time(Y, trans_time_n(Y,i,TY), TY)
        sum += w1 * w2
    tauex = float(sum)/float(observ_t_Y)
    return tauex
    
    
def get_ij_ca(res_i,res_j,fs):
    protein = u.select_atoms('backbone')
    phi_sel = protein.residues[res_i].phi_selection()
    phi_sel2 = protein.residues[res_j].phi_selection()
    resi_phi = []
    resj_phi = []
    for ts in u.trajectory:
        frame = ts.frame
        if frame >= fs:
            break
        k = new_dihedral(phi_sel.positions)
        resi_phi.append(k)
        p = new_dihedral(phi_sel2.positions)
        resj_phi.append(p)
    #get TX, TY, lists of transition time and pass them
    X = np.array(resi_phi)
    Y = np.array(resj_phi)
    TX = trans_time(X)
    TY = trans_time(Y)
    CA = (-1) * math.log(tau_ex(X, Y, TX, TY)/tau_p(X, TX))
    #CA = get_ca(np.array(resi_phi),np.array(resj_phi),TX,TY)
    pair = str(res_i) + '-' + str(res_j) + '.ax'
    all = str(res_i) + '\t' + str(res_j) + '\t' + str(CA) + '\n'
    f1 = open(pair, 'w')
    f1.write(all)
    f1.close()


def main():
    with ProcessPoolExecutor(max_workers=usecpus) as executer:
        a = []
        for i, j in permutations(range(2,8), 2):
            future = executer.submit(get_ij_ca, i, j, frms_num)
            a.append(future)
        wait(a)
    #join small files together
    os.system('cat *.ax > temp-all')
    f2 = open("temp-all")
    for line in f2.readlines():
        a = re.split('\t|\n',line)
        s0 = int(a[0])
        s1 = int(a[1])
        s2 = float(a[2])
        b[s0][s1] = s2
        #b[s1][s0] = s2
    f2.close()

    for i in range(0,352):
        for j in range(0,352):
            p = str(i) + '\t' + str(j) + '\t' + str(b[i][j]) + '\n'
            f.write(p)
    f.close()
    os.system('mv *.ax crap/')
    os.system('rm temp-all')


if __name__ == '__main__':
    main()
