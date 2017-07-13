# -*- coding: utf-8 -*-
"""
Title:    Himmelblau.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Himmelblau test function for design optimization

xOpt = [3.0, 2.0]
xOpt = [-2.805118, 3.131312]
xOpt = [-3.779310, -3.283186]
xOpt = [3.584428, -1.848126]
fOpt = 0.0
-------------------------------------------------------------------------------
"""
from __future__ import absolute_import, division, print_function
from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    f = (x[0]**2+x[1]-11)**2+(x[0]+x[1]**2-7)**2
    g = []
    return(f, g)


x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, Output = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                            Alg="PyGMO_de", StatusReport=False,
                            OptNameAdd="Himmelblau", ResultReport=False,
                            deltax=1e-6)
