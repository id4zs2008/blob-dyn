#!/usr/bin/python
import numpy as np
from sklearn.metrics import mutual_info_score
import MDAnalysis
from itertools import combinations_with_replacement,permutations
from concurrent.futures import ProcessPoolExecutor, Future, wait
import os
import re

usecpus = 16 #how many cores to use
frms_num = 100
u = MDAnalysis.Universe('ini.pdb','allpdb.trr')
bins = 30
r = (-180, 180)
range_dihe = np.vstack((list(r),list(r)))
f = open('MI-out.txt', 'w')
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


def auto_mi(nbins, X, Y, range=2*[[-180., 180.]]):
    joint_mtx =  np.histogram2d(X,Y,nbins,range=range)
    return mutual_info_score(None, None, contingency=joint_mtx[0])


def get_ij_mi(res_i,res_j,fs):
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
    mi = auto_mi(bins, resi_phi, resj_phi,range=range_dihe)
    pair = str(res_i) + '-' + str(res_j) + '.px'
    all = str(res_i) + '\t' + str(res_j) + '\t' + str(mi) + '\n'
    f1 = open(pair, 'w')
    f1.write(all)
    f1.close()


def main():
    pt = 0
    with ProcessPoolExecutor(max_workers=usecpus) as executer:
        a = []
        for i, j in combinations_with_replacement(range(2,351), 2):
        #for i, j in combinations_with_replacement(range(2,20), 2):
            future = executer.submit(get_ij_mi, i, j, frms_num)
            a.append(future)
            pt += 1
            if pt%1000 == 0:
                print pt
        wait(a)
    #join small files together
    os.system('cat *.px > temp-all')
    f2 = open("temp-all")
    for line in f2.readlines():
        a = re.split('\t|\n',line)
        s0 = int(a[0])
        s1 = int(a[1])
        s2 = float(a[2])
        b[s0][s1] = s2
        b[s1][s0] = s2
    f2.close()

    for i in range(0,352):
        for j in range(0,352):
            p = str(i) + '\t' + str(j) + '\t' + str(b[i][j]) + '\n'
            f.write(p)
    f.close()
    os.system('rm *.px')
    os.system('rm temp-all')


if __name__ == '__main__':
    main()

