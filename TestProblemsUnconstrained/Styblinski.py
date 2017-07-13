# -*- coding: utf-8 -*-
"""
Title:    Styblinski.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Styblinski--Tang test function for design optimization

xOpt = []
fOpt = 0.0
-------------------------------------------------------------------------------
"""
from __future__ import absolute_import, division, print_function
from DesOptPy import DesOpt
import numpy as np


def SysEq(z, gc):
    x = z[0]
    y = z[1]
    f = 1./2.*(x**4-16.*x**2+5.*x+y**4-16.*y**2+5.*y)
    g = []
    return(f, g)

x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, Output = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                            Alg="SLSQP", deltax=1e-3, StatusReport=False,
                            OptNameAdd="Styblinski", ResultReport=False)
