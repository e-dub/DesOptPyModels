# -*- coding: utf-8 -*-
"""
Title:    Schwefel.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Schwefel test function for design optimization

xOpt = []
fOpt = 0.0
-------------------------------------------------------------------------------
"""
from __future__ import absolute_import, division, print_function
from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    n = len(x)
    f = 418.9829*n-np.sum(x*np.sin(np.sqrt(np.abs(x))))
    g = []
    return(f, g)

x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, Output = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                            Alg="SLSQP", StatusReport=False,
                            OptNameAdd="Schwefel", ResultReport=False,
                            deltax=1e-6)
