# -*- coding: utf-8 -*-
"""
Title:    RosenbrockSens.py
Units:    -
Author:   E.J. Wehrle
Date:     July 5, 2016
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Description:

Rosenbrock test function for design optimization

xOpt = []
fOpt = 0.0

------------------------------------------------------------------------------------------------------------------------
"""
import sys
import numpy as np
from DesOptPy import DesOpt
from DesOptPy import OptAlgOptions
from scipy.optimize import rosen, rosen_der


def SysEq(x, gc):
    f = rosen(x)
    return f, []

    
def SensEq(x, f, g, gc):
    dfdx = rosen_der(x)
    return dfdx, []


x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                        Alg="SLSQP", StatusReport=False,   DesVarNorm=True, OptNameAdd="FDUnnorm",
                        DoE=False, SBDO=False, ResultReport=False,
                        deltax=1e-6)
