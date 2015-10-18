# -*- coding: utf-8 -*-
"""
Title:    DixonPrice.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Description:

Dixon Price test function for design optimization

xOpt = []
fOpt = 0.0

------------------------------------------------------------------------------------------------------------------------
"""

from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    n = len(x)
    j = np.arange(2, n+1)
    x2 = 2 * x**2
    f = np.sum(j*(x2[1:]-x[:-1])**2)+(x[0]-1)**2
    g = []
    return(f, g)


x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                        Alg="SLSQP", StatusReport=True, OptNameAdd="DixonPrice",
                        DoE=False, SBDO=False, ResultReport=False,
                        deltax=1e-6)
