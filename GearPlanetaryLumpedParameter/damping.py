from __future__ import division
import numpy as np


def damping(N, Rot, mp):
    D = np.zeros([N, N])
    # PIANETA 1
    # etap1
    D[4, 5] = 2*Rot*mp
    # xip1
    D[5, 5] = -2*Rot*mp
    # PIANETA 2
    # etap2
    D[7, 8] = 2*Rot*mp
    # xip2
    D[8, 8] = -2*Rot*mp
    # PIANETA 3
    # etap3
    D[10, 11] = 2*Rot*mp
    # xip3
    D[11, 11] = -2*Rot*mp
    return D
