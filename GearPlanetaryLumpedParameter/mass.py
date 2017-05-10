from __future__ import division
import numpy as np


def mass(N, ms, Js, mp, Jp, mc, Jc, mr, Jr, rc):
    M = np.zeros([N, N])
    # SOLARE
    # tetas
    M[0, 0] = Js
    # xs
    M[1, 1] = ms
    #M[1, 1] = -ms
    # ys
    M[2, 2] = ms
    #M[2, 2] = -ms

    # PIANETA 1
    # tetap1
    M[3, 3] = Jp
    # etap1
    M[4, 4] = -mp
    # xip1
    M[5, 5] = -mp

    # PIANETA 2
    # tetap2
    M[6, 6] = Jp
    # etap2
    M[7, 7] = -mp
    # xip2
    M[8, 8] = -mp

    # PIANETA 3
    # tetap3
    M[9, 9] = Jp
    # etap3
    M[10, 10] = -mp
    # xip3
    M[11, 11] = -mp

    # PORTATRENO
    # thetac
    M[12, 12] = -Jc-3*mp*rc**2
    # xc
    #M[13, 13] = -mc -3*mp?
    M[13, 13] = -mc
    # yc
    #M[14, 14] = -mc -3*mp?
    M[14, 14] = -mc

    # RING
    # tetar
    M[15, 15] = Jr
    #M[15, 15] = -Jr
    # xr
    M[16, 16] = -mr
    # yr
    M[17, 17] = -mr
    return M
