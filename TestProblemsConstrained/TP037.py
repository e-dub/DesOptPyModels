# -*- coding: utf-8 -*-
"""
Title:    TP037.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Test function for design optimization

-------------------------------------------------------------------------------
"""
from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    f = -x[0]*x[1]*x[2]
    g = np.ones([2, ])*0.
    g[0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
    g[1] = -x[0] - 2.*x[1] - 2.*x[2]
    return(f, g)


x0 = np.ones([3, ])*4.
xL = np.ones([3, ])*0.
xU = np.ones([3, ])*42.
gc = np.ones([2, ])*0.
xOpt, fOpt, Output = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                            Alg="NLPQLP", DesVarNorm=True, deltax=1e-6,
                            Debug=False, ResultReport=False,
                            StatusReport=False)
