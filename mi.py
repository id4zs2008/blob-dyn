import numpy as np
from sklearn.metrics import mutual_info_score

"""
A mini version to get mi
Sample output:
mi = 0.63
"""


def auto_mi(nbins, X, Y, range=2*[[-180., 180.]]):
    joint_mtx =  np.histogram2d(X,Y,nbins,range=range)
    return mutual_info_score(None, None, contingency=joint_mtx[0])


def main():
    x = [-4,-2, 0, 2, 3,-3, 2, 3]
    y = [ 5, 4, 2,-3, 1, 3, 2, 4]
    bins = 5
    range1 = (-5, 5)
    range2 = np.vstack((list(range1),list(range1)))
    mi = auto_mi(bins, x, y,range=range2)
    print "mi = %.2f" % (mi)


if __name__ == '__main__':
    main()
