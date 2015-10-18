# -*- coding: utf-8 -*-
"""
Title:    Rastrigin.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Description:

Rastrigin test function for design optimization

xOpt = []
fOpt = 0.0

------------------------------------------------------------------------------------------------------------------------
"""

from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    n = len(x)
    f = 10*n+np.sum(x**2-10*np.cos(2.*np.pi*x))
    g = []
    return(f, g)


x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                        Alg="SLSQP", StatusReport=True, OptNameAdd="Rastrigin",
                        DoE=False, SBDO=False, ResultReport=True,
                        deltax=1e-6)
