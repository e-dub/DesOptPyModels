# -*- coding: utf-8 -*-
"""
Title:    Ackley.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Ackley test function for design optimiazation

xOpt = [0.0, 0.0]
fOpt = 0.0

-------------------------------------------------------------------------------
"""
import numpy as np
from DesOptPy import DesOpt


def SysEq(x, gc):
    a = 20
    b = 0.2
    c = 2*np.pi
    n = len(x)
    s1 = sum(x**2)
    s2 = sum(np.cos(c * x))
    f = -a*np.exp(-b*np.sqrt(s1 / n)) - np.exp(s2 / n) + a + np.exp(1)
    g = []
    return(f, g)

x0 = np.ones([2, ])*2
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, Output = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                            DesVarNorm=True, Alg="NLPQLP", deltax=1e-6,
                            StatusReport=True, ResultReport=True,
                            OptNameAdd="Ackely")
