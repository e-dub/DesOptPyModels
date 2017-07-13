# -*- coding: utf-8 -*-
"""
Title:    Michalewicz.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Michalewicz test function for design optimization

xOpt = []
fOpt = 0.0

-------------------------------------------------------------------------------
"""
from __future__ import absolute_import, division, print_function
from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    michalewicz_m = 2
    n = len(x)
    j = np.arange(1., n+1)
    f = -np.sum(np.sin(x)*np.sin(j*x**2/np.pi)**(2.*michalewicz_m))
    g = []
    return(f, g)


x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                        Alg="SLSQP", deltax=1e-6, StatusReport=False,
                        OptNameAdd="Michalewicz", ResultReport=True)
