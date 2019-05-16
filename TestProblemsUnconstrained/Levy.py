# -*- coding: utf-8 -*-
"""
Title:    Levy.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Levy test function for design optimization

xOpt = []
fOpt = 0.0

-------------------------------------------------------------------------------
"""
from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    n = len(x)
    z = 1.+(x-1.)/4.
    f = (np.sin(np.pi*z[0])**2 +
         np.sum((z[:-1]-1.)**2*(1. + 10.*np.sin(np.pi*z[:-1] + 1.)**2)) +
         (z[-1]-1.)**2*(1. + np.sin(2.*np.pi*z[-1])**2))
    g = []
    return(f, g)


x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                        Alg="NLPQLP", StatusReport=False, OptNameAdd="Levy",
                        ResultReport=False, deltax=1e-6)
