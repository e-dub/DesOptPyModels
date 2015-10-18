# -*- coding: utf-8 -*-
"""
Title:    Eggholder.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Description:

Eggholder test function for design optimization

xOpt = []
fOpt = 0.0

------------------------------------------------------------------------------------------------------------------------
"""

from DesOptPy import DesOpt
import numpy as np


def SysEq(z, gc):
    x = z[0]
    y = z[1]
    return(-(y+47.)*np.sin(np.sqrt(np.abs(y+x/2.+47.)))-x*np.sin(np.sqrt(np.abs(x-(y+47.)))), [])


x0 = np.zeros([2, ])
xL = np.ones([2, ])*-5
xU = np.ones([2, ])*5
gc = []
xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                        Alg="SLSQP", StatusReport=True, OptNameAdd="Eggholder",
                        DoE=False, SBDO=False, ResultReport=False,
                        deltax=1e-3)
