# -*- coding: utf-8 -*-
"""
Title:    Griewank.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Griewank test function for design optimization

xOpt = []
fOpt = 0.0
-------------------------------------------------------------------------------
"""
from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    fr = 4000
    n = len(x)
    j = np.arange(1., n+1)
    s = np.sum(x**2)
    p = np.prod(np.cos(x/np.sqrt(j)))
    f = s/fr - p + 1
    g = []
    return(f, g)

x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, Output = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                        Alg="SLSQP", StatusReport=True, OptNameAdd="Griewank",
                        ResultReport=False,
                        deltax=1e-6)
