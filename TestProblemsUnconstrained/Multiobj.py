# -*- coding: utf-8 -*-
"""
Title:    Multiobj.py
Units:    -
Author:   E. J. Wehrle
Date:     July 25, 2015
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Description:
Simple multi-objective test function for design optimization with use of weighting factors resulting in a Pareto front

------------------------------------------------------------------------------------------------------------------------
"""

from DesOptPy import DesOpt
import numpy as np

xL = np.ones([1,])*0
xU = np.ones([1,])*10
x0 = np.average((xL, xU), 0)
gc = []
Alg = "SLSQP"
nPareto = 21
xOpt = [[]]*nPareto
fOpt = [[]]*nPareto
f1 = [[]]*nPareto
f2 = [[]]*nPareto
nEval = [[]]*nPareto
gamma = np.logspace(-1, 1, nPareto)
for ii in range(nPareto):
    def SysEq(x, gc):
        f1[ii] = x
        f2[ii] = 100.0-x**2
        f = -1*((gamma[ii]*f1[ii]) + f2[ii])
        g = []
        return f, g


    xOpt[ii], fOpt[ii], SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq, DesVarNorm=True,
                                    Alg=Alg, deltax=1e-6, DoE=None, SBDO=False, OptNameAdd=str(ii),
                                    Debug=True, ResultReport=False, StatusReport=False, Alarm=False, PrintOut=False)
    x0 = xOpt[ii]  # Better start value
    SysEq(xOpt[ii], fOpt[ii])  # Save the optimal values of f1 and f2 instead of last values of optimization
ParetoFront = np.array([f1, f2])
print ParetoFront
