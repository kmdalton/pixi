#!/usr/bin/env python


import numpy as np
from sys import argv
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_parm(parm):
    A,B,C = parm.A,parm.B,parm.C
    A = A/np.linalg.norm(A)
    B = B/np.linalg.norm(B)
    C = C/np.linalg.norm(C)

    f = plt.figure()
    a = Axes3D(f)

    a.plot3D([0,1], [0,0], [0,0], c='b')
    a.plot3D([0,0], [0,1], [0,0], c='g')
    a.plot3D([0,0], [0,0], [0,1], c='r')

    a.text3D(1, 0, 0, 'X', color='b')
    a.text3D(0, 1, 0, 'Y', color='g')
    a.text3D(0, 0, 1, 'Z', color='r')

    a.plot3D([0, A[0]], [0, A[1]], [0, A[2]], c='b', linestyle='--')
    a.plot3D([0, B[0]], [0, B[1]], [0, B[2]], c='g', linestyle='--')
    a.plot3D([0, C[0]], [0, C[1]], [0, C[2]], c='r', linestyle='--')
    a.text3D(A[0], A[1], A[2], "A", color='b')
    a.text3D(B[0], B[1], B[2], "B", color='g')
    a.text3D(C[0], C[1], C[2], "C", color='r')
