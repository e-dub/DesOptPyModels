# -*- coding: utf-8 -*-
"""
Title:              CalcGear.py
Units:              N, mm, t, MPa, s
Author:             E.J. Wehrle
Date:               May 4, 2017
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

-------------------------------------------------------------------------------
"""
from __future__ import absolute_import, division, print_function
import numpy as np
import Solver
from PlanetGear import CalcMatrices
    
GearSet = "TQ070i3"
GearSet = "MP080i6"
GearSet = "Benchmark"

if GearSet=="TQ070i3":
    print("TQ070i3")
    rS=32.4/2
    rC=24.5
    rR=64.8/2
    rP=15.6/2
    mS = 0.262e-3
    JS = 21e-3
    mC = 0.725e-3
    JC = 189e-3
    mR = 1.87e-3
    JR = 3131e-3
    mP = 0.025e-3
    JP = 9.11e-3
    kSx = 213e3
    kSy = 213e3
    kStheta = 148e3
    kCx = 270e3
    kCy = 270e3
    kCtheta = 125e3
    kRx = 650e8
    kRy = 650e8
    kRtheta = 1e10
    kSP = 213e3
    kCP = 1250e3
    kRP = 207e3
    nP = 3.
    alpha = 20.
    thetaCdot = 1000./60.*2*np.pi
    thetaCdot = 0
    thetaCdotdot = 0
elif GearSet=="MP080i6":
    print("MP080i6")
    rS = 32.4/2.
    rC = 24.5
    rR = 64.8/2
    rP = 15.6/2
    mS = 0.262e-3
    JS = 21e-3
    mC = 0.44e-3
    JC = 189e-3
    mR = 2e-3
    JR = 3270e-3
    mP = 0.025e-3
    JP = 9.11e-3
    kSx = 170e3
    kSy = 170e3
    kStheta = 4700e3
    kCx = 270e3
    kCy = 270e3
    kCtheta = 32000e3
    kRx = 650e8
    kRy = 650e8
    kRtheta = 1e10
    kSP = 170e3
    kCP = 1250e3
    kRP = 170e3
    nP=3.
    alpha=20.
    thetaCdot=1000./60.*2*np.pi
    thetaCdotdot=0
elif GearSet=="Benchmark":
    print("New benchmark")
    rS = 50.
    rC = 75.
    rR = 100.
    rP = 25.
    mS = 1.0e-3
    JS = 1.0e-3
    mC = 5.0e-3
    JC = 5.0e-3
    mR = 5.0e-3
    JR = 5.0e-3
    mP = 0.5e-3
    JP = 0.5e-3
    kSx = 200e3
    kSy = 200e3
    kStheta = 5000e3
    kCx = 250e3
    kCy = 250e3
    kCtheta = 5000e3
    kRx = 5000e3
    kRy = 5000e3
    kRtheta = 5000e3
    kSP = 200e3
    kCP = 1000e3
    kRP = 200e3
    nP=3.
    alpha=20.
    thetaCdot=20000./60.*2*np.pi
    thetaCdot=0
    thetaCdotdot=0
M, D, K = CalcMatrices(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                       kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                       kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                       kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                       rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                       thetaCdotdot=thetaCdotdot)
omega, Phi = Solver.ModalDamped(M=M, D=D, K=K, Hz=True, All=True, Sum=True)    
print("\nOmega damped")
for ii in omega:
    print(ii)
