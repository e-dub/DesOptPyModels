# -*- coding: utf-8 -*-
"""
Title:    CylinderHead.py
Units:    -
Author:   E.J. Wehrle
Date:     November 30, 2014
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Cylinder head test function for design optimization

Source:  Dakota User's Guide §20.4

xOpt = [2.122, 1.769]
fOpt = -2.461
-------------------------------------------------------------------------------
"""
from __future__ import absolute_import, division, print_function
from DesOptPy import DesOpt
import numpy as np


def SysEq(x, gc):
    exhaust_offset = 1.34
    exhaust_dia = 1.556
    intake_offset = 3.25
    warranty = 100000. + 15000. * (4. - x[1])
    cycle_time = 45. + 4.5*pow(4. - x[1], 1.5)
    wall_thickness = intake_offset-exhaust_offset-(x[0]+exhaust_dia)/2.
    horse_power = 250.+200.*(x[0]/1.833-1.)
    max_stress = 750. + pow(np.fabs(wall_thickness), -2.5)
    f = -1.*(horse_power/250.+warranty/100000)
    g1 = max_stress/1500.-1.
    g2 = 1.-warranty/100000.
    g3 = cycle_time/60. - 1.
    return(f, [g1, g2, g3])


xL = np.array([1.5, 0.0])
xU = np.array([2.164, 4.0])
gc = np.array([0.0, 0.0, 0.0])
x0 = np.array([1.8, 1.0])
xOpt, fOpt, Output = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq,
                            DesVarNorm=True, Alg="COBYLA", deltax=1e-6,
                            Debug=False, ResultReport=False,
                            StatusReport=False, OptNameAdd="CylinderHead")

