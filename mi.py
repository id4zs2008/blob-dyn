#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pylab as pb
from matplotlib.collections import LineCollection
from mdentropy.entropy import mi


"""
adopted from the following source
http://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy
compared mi1(stackoverflow formula, using log2 but changed to natural log for comparison) and mi2 (mdentropy module, natural log)
FOUND that mi1 is exactly equal to mi2 using the same data, indicating that the algorithm implementation is good
"""
def calc_MI(X,Y,bins,range1,range2):
    c_XY = np.histogram2d(X,Y,bins,range=range2)[0]
    c_X = np.histogram(X,bins,range=range1)[0]
    c_Y = np.histogram(Y,bins,range=range1)[0]
    H_X = shan_entropy(c_X)
    H_Y = shan_entropy(c_Y)
    H_XY = shan_entropy(c_XY)
    MI = H_X + H_Y - H_XY
    return MI


def shan_entropy(c):
    c_normalized = c / float(np.sum(c))
    c_normalized = c_normalized[np.nonzero(c_normalized)]
    H = -sum(c_normalized* np.log(c_normalized))
    return H


def main():
    x = [-4,-2, 0, 2, 3,-3, 2, 3]
    y = [ 5, 4, 2,-3, 1, 3, 2, 4]
    bins = 5
    range1 = (-5, 5)
    range2 = ((-5, 5),(-5, 5))
    range3 = np.vstack((list(range1),list(range1)))
    mi1 = calc_MI(x, y, bins,range1,range2)
    mi2 = mi(bins, x, y,range=range3)
    print "mi1 = %.2f" % (mi1)
    print "mi2 = %.2f" % (mi2)


if __name__ == '__main__':
    main()

