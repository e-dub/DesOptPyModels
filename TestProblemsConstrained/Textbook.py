# -*- coding: utf-8 -*-
"""
Title:    Textbook.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Testbook test function for design optimization

Source:  https://dakota.sandia.gov//sites/default/files/docs/6.1/html-ref/textbook.html

Unconstrained:
xOpt = [1.0, 1.0]
fOpt = 0.0

Constrained:
xOpt = [0.5, 0.5]
fOpt = 0.125

-------------------------------------------------------------------------------
"""

from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    f = (x[0]-1.)**4 + (x[1]-1.)**4
    g1 = x[0]**2 - x[1]/2.
    g2 = x[1]**2 - x[0]/2.
    return(f, [g1, g2])


xL = np.array([0.5, -2.9])
xU = np.array([5.8, 2.9])
gc = np.array([0.0, 0.0])
x0 = xU
xOpt, fOpt, Output = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                            Alg="SLSQP", StatusReport=True, DesVarNorm="xLxU",
                            ResultReport=False, deltax=1e-10,
                            OptNameAdd="Textbook")
