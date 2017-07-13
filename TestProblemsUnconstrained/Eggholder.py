# -*- coding: utf-8 -*-
"""
Title:    Eggholder.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Eggholder test function for design optimization

xOpt = []
fOpt = 0.0

-------------------------------------------------------------------------------
"""
from __future__ import absolute_import, division, print_function
from DesOptPy import DesOpt
import numpy as np
from numpy import sin, sqrt, abs, ones, zeros

def SysEq(z, gc):
    x = z[0]
    y = z[1]
    f = -(y+47.)*sin(np.sqrt(abs(y+x/2.+47.)))-x*sin(sqrt(abs(x-(y+47.))))
    return(f, [])


x0 = zeros([2, ])
xL = ones([2, ])*-5
xU = ones([2, ])*5
gc = []
xOpt, fOpt, Output = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                            deltax=1e-3, Alg="SLSQP", StatusReport=False,
                            ResultReport=False, OptNameAdd="Eggholder")
