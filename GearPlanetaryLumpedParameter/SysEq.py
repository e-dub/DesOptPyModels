# -*- coding: utf-8 -*-
"""
Title:              SysEq.py
Units:              N, mm, t, MPa, s
Author:             E.J. Wehrle
Date:               May 4, 2017
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:
    Modal solver for damped and undamped systems
    ...coming soon to an optimization near you: analytical sensitivities!

Important:
    Jmax = m*r**2/2
-------------------------------------------------------------------------------
"""
from __future__ import absolute_import, division, print_function
from DesOptPy import DesOpt, OptAlgOptions
from Solver import ModalSensUndampedHaftka1, rad2hz, ModalSensDamped1
from LumPy import Solver
import PlanetaryGearDatabase
from LumPy.PlanetaryGear import CalcMatrices, PlotGearSet 
import numpy as np


GearSet = "TQ070i3"
#GearSet = "MP080i6"
#GearSet = "Me"
#GearSet = "TurboProp"
#GearSet = "HeliOH58"
GearSet = "Benchmark"
#GearSet = "EricsonParker2013B"
#GearSet = "EricsonParker2013B"
iEval =0
global iEval

def SysEval(mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, 
            kCtheta, kRx, kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, 
            rP, thetaCdot, thetaCdotdot, Hz=True, Sum=True, All=False, 
            Damping=True, Plot=False):
    M, D, K, dof = CalcMatrices(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, 
                                mP=mP, JP=JP, kSx=kSx, kSy=kSy, 
                                kStheta=kStheta, kCx=kCx, kCy=kCy, 
                                kCtheta=kCtheta, kRx=kRx, kRy=kRy, 
                                kRtheta=kRtheta, kSP=kSP, kCP=kCP, kRP=kRP, 
                                nP=nP, alpha=alpha, rS=rS, rC=rC, rR=rR, rP=rP, 
                                thetaCdot=thetaCdot, thetaCdotdot=thetaCdotdot)
    if Damping:
        omega, Phi = Solver.ModalDamped(M=M, D=D, K=K, Hz=Hz, All=All, Sum=Sum)    
    else:
        omega, Phi = Solver.ModalUndamped(M=M, D=D, K=K,  Hz=Hz, All=All, 
                                          Sum=Sum) 
    if Plot:
        PlotGearSet(nP, rS, rC, rR, rP, Phi, factor=10.0, FileType="svg")
    global iEval
    iEval += 1
    return(omega, Phi)


def SysEqI(x, gc):
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = x[0]
    JS = x[1]
    mC = x[2]
    JC = x[3]
    mR = x[4]
    JR = x[5]
    mP = x[6]
    JP = x[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                      kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                      kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                      kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                      rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                      thetaCdotdot=thetaCdotdot, Hz=True, Sum=True, All=False, 
                      Damping=True)
    f = -f0[0]
    g = []
    return(f, g)

def SensEqI(x, f, g, gc):
    dfdx = np.zeros(np.shape(x))
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = x[0]
    JS = x[1]
    mC = x[2]
    JC = x[3]
    mR = x[4]
    JR = x[5]
    mP = x[6]
    JP = x[7]
    M, D, K, dof = CalcMatrices(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, 
                                JR=JR, mP=mP, JP=JP, kSx=kSx, kSy=kSy, 
                                kStheta=kStheta, kCx=kCx, kCy=kCy, 
                                kCtheta=kCtheta, kRx=kRx, kRy=kRy, 
                                kRtheta=kRtheta, kSP=kSP, kCP=kCP, 
                                kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                                rC=rC, rR=rR, rP=rP, 
                                thetaCdot=thetaCdot, thetaCdotdot=0)
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                      kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                      kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                      kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                      rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                      thetaCdotdot=0, Hz=False, Sum=True, All=False, 
                      Damping=True) 
    for ii in range(len(x)):
        dprocent = 1e-3
        import copy
        xN = copy.deepcopy(x)
        xN[ii] *= (1+dprocent)
        dx = x[ii]*dprocent
        mS = xN[0]
        JS = xN[1]
        mC = xN[2]
        JC = xN[3]
        mR = xN[4]
        JR = xN[5]
        mP = xN[6]
        JP = xN[7]
        M1, D1, K1, dof = CalcMatrices(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, 
                                       JR=JR, mP=mP, JP=JP, kSx=kSx, kSy=kSy, 
                                       kStheta=kStheta, kCx=kCx, kCy=kCy, 
                                       kCtheta=kCtheta, kRx=kRx, kRy=kRy, 
                                       kRtheta=kRtheta, kSP=kSP, kCP=kCP, 
                                       kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                                       rC=rC, rR=rR, rP=rP, 
                                       thetaCdot=thetaCdot, thetaCdotdot=0)
        dMdx = (M1-M)/dx
        dDdx = (D1-D)/dx
        dKdx = (K1-K)/dx
        domegadx = ModalSensUndampedHaftka1(f0, Phi, M, K, dMdx, dKdx, n=1)
        #domegadx = ModalSensDamped1(f0, Phi, M, D, K, dMdx, dDdx, dKdx, n=1)
        dfdx[ii] = rad2hz(domegadx[0])
    return dfdx, []


def SysEqII(x, gc):
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = x[0]
    JS = x[1]
    mC = x[2]
    JC = x[3]
    mR = x[4]
    JR = x[5]
    mP = x[6]
    JP = x[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                       kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                       kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                       kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                       rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                       thetaCdotdot=thetaCdotdot)     

#    fBand = np.array([[1200.0, 1800.0],
#                      [2400.0, 3600.0]])
#    fBand = np.array([[ 300.0,  450.0],
#                      [ 600.0,  900.0],
#                      [1200.0, 1800.0],
#                      [2400.0, 3600.0],
#                      [4800.0, 7200.0]])
    fBand = np.array([[  75.0,  100.0],
                      [ 150.0,  200.0],
                      [ 300.0,  400.0],
                      [ 600.0,  800.0],
                      [1200.0, 1600.0],
                      [2400.0, 3200.0]])
    g = np.zeros(np.shape(gc))
    ig = 0
    nf0 = len(f0)
    for ii in range(nf0):
        for jj in range(np.shape(fBand)[0]):
            g[ig] = np.prod((f0[ii]/fBand[jj,0]-1, 
                             fBand[jj,1]/f0[ii]-1))*(fBand[jj,1]+fBand[jj,0])/2.
            ig += 1
    f = x[0]+x[2]+x[4]+3*x[6]
    return f, g

def SensEqII(x, f, g, gc):
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = x[0]
    JS = x[1]
    mC = x[2]
    JC = x[3]
    mR = x[4]
    JR = x[5]
    mP = x[6]
    JP = x[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                       kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                       kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                       kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                       rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                       thetaCdotdot=thetaCdotdot)    
#    fBand = np.array([[1200.0, 1800.0],
#                      [2400.0, 3600.0]])
#    fBand = np.array([[ 300.0,  450.0],
#                      [ 600.0,  900.0],
#                      [1200.0, 1800.0],
#                      [2400.0, 3600.0],
#                      [4800.0, 7200.0]])
    fBand = np.array([[  75.0,  100.0],
                      [ 150.0,  200.0],
                      [ 300.0,  400.0],
                      [ 600.0,  800.0],
                      [1200.0, 1600.0],
                      [2400.0, 3200.0]])
    g = np.zeros(np.shape(gc))
    ig = 0
    nf0 = len(f0)
    for ii in range(nf0):
        for jj in range(np.shape(fBand)[0]):
            dgdx[ig] = -c*fBand[jj,1]*(-1 + f0[ii]/fBand[jj,0])*Derivative(omega(x), x)/omega(x)**2 + c*(omega_U/omega(x) - 1)*Derivative(omega(x), x)/omega_L
            #g[ig] = np.prod((f0[ii]/fBand[jj,0]-1, fBand[jj,1]/f0[ii]-1))*fBand[jj,1]+fBand[jj,0])/2.
            ig += 1
    dfdx = np.ones(np.shape(x))
    return dfdx, dgdx




def SysEqIII(x, gc):
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = x[0]
    JS = x[1]
    mC = x[2]
    JC = x[3]
    mR = x[4]
    JR = x[5]
    mP = x[6]
    JP = x[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                       kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                       kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                       kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                       rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                       thetaCdotdot=thetaCdotdot)     

#    fBand = np.array([[1200.0, 1800.0],
#                      [2400.0, 3600.0]])
#    fBand = np.array([[ 300.0,  450.0],
#                      [ 600.0,  900.0],
#                      [1200.0, 1800.0],
#                      [2400.0, 3600.0],
#                      [4800.0, 7200.0]])
    fBand = np.array([[  75.0,  100.0],
                      [ 150.0,  200.0],
                      [ 300.0,  400.0],
                      [ 600.0,  800.0],
                      [1200.0, 1600.0],
                      [2400.0, 3200.0]])
    g = np.zeros(np.shape(gc))
    ig = 0
    nf0 = len(f0)
    for ii in range(nf0):
        for jj in range(np.shape(fBand)[0]):
            g[ig] = np.prod((f0[ii]/fBand[jj,0]-1, 
                             fBand[jj,1]/f0[ii]-1))*(fBand[jj,1]+fBand[jj,0])/2.
            ig += 1
    f = -f0[0]
    return f, g

#ModalSensUndampedHaftka

def FuzzySysEq(p, x, ir):
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = x[0]
    JS = x[1]
    mC = x[2]
    JC = x[3]
    mR = x[4]
    JR = x[5]
    mP = x[6]
    JP = x[7]
    kSP = p[0]
    kRP = p[1]
    kSx = p[2] 
    kSy = p[2]
    kStheta = p[3]
    kCx = p[4]
    kCy = p[4]
    kCtheta = p[5]
    kRx = p[6] 
    kRy = p[6]
    kRtheta = p[7]
    kCP = p[8]
    thetaCdot=0
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                      kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                      kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                      kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                      rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                      thetaCdotdot=thetaCdot,  Hz=True, Sum=True, All=0, 
                      Damping=1) 
    return(f0[ir]*10000)


def FuzzySensEq(p, r, g, x, ir):
    dfdx = np.zeros(np.shape(p))
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    kSP = p[0]
    kRP = p[1]
    kSx = p[2] 
    kSy = p[2]
    kStheta = p[3]
    kCx = p[4]
    kCy = p[4]
    kCtheta = p[5]
    kRx = p[6] 
    kRy = p[6]
    kRtheta = p[7]
    kCP = p[8]
    M, D, K, dof = CalcMatrices(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, 
                                JR=JR, mP=mP, JP=JP, kSx=kSx, kSy=kSy, 
                                kStheta=kStheta, kCx=kCx, kCy=kCy, 
                                kCtheta=kCtheta, kRx=kRx, kRy=kRy, 
                                kRtheta=kRtheta, kSP=kSP, kCP=kCP, 
                                kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                                rC=rC, rR=rR, rP=rP, 
                                thetaCdot=thetaCdot, thetaCdotdot=0)
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                      kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                      kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                      kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                      rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                      thetaCdotdot=thetaCdot, Hz=False, Sum=True, All=0, 
                      Damping=1) 
    for ii in range(len(p)):
        dprocent = 0.0001
        import copy
        pN = copy.deepcopy(p)
        pN[ii] *= (1+dprocent)
        dx = p[ii]*dprocent
        kSP = pN[0]
        kRP = pN[1]
        kSx = pN[2] 
        kSy = pN[2]
        kStheta = pN[3]
        kRx = pN[4] 
        kRy = pN[4]
        kRtheta = pN[5]
        kCx = pN[6]
        kCy = pN[6]
        kCtheta = pN[7]
        kCP = pN[8]
        M1, D1, K1, dof = CalcMatrices(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, 
                                       JR=JR, mP=mP, JP=JP, kSx=kSx, kSy=kSy, 
                                       kStheta=kStheta, kCx=kCx, kCy=kCy, 
                                       kCtheta=kCtheta, kRx=kRx, kRy=kRy, 
                                       kRtheta=kRtheta, kSP=kSP, kCP=kCP, 
                                       kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                                       rC=rC, rR=rR, rP=rP, 
                                       thetaCdot=thetaCdot, thetaCdotdot=0)
        dMdx = (M1-M)/dx
        dDdx = (D1-D)/dx
        dKdx = (K1-K)/dx
        #domegadx = ModalSensUndampedHaftka1(f0, Phi, M, K, dMdx, dKdx, n=1, ir=ir)
        domegadx = ModalSensDamped1(f0, Phi, M, D, K, dMdx, dDdx, dKdx, n=1, 
                                    ir=ir)
        #print(domegadx)
        dfdx[ii] = rad2hz(domegadx)*10000
    #print(dfdx)
    return dfdx, []

def SysEqIunc(x, gc):
    from FuzzAnPy import FuzzAn
    from FuzzAnPy import FuzzyNumber
    nAlpha = 1
    nr = 18
    pFuzz = np.zeros([9, nAlpha, 2])
    pFuzz[0, :, :] = FuzzyNumber.Int(150e3, 250e3, nAlpha)
    pFuzz[1, :, :] = FuzzyNumber.Int(150e3, 250e3, nAlpha)
    pFuzz[2, :, :] = FuzzyNumber.Int(190e3, 210e3, nAlpha)
    pFuzz[3, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[4, :, :] = FuzzyNumber.Int(240e3, 260e3, nAlpha)
    pFuzz[5, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[6, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[7, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[8, :, :] = FuzzyNumber.Int(900e3, 1100e3, nAlpha)
    rFuzz, OutputData = FuzzAn(FuzzySysEq=FuzzySysEq, pFuzz=pFuzz, nr=nr,
                               Alg="NLPQLP", nAlpha=nAlpha, deltax=1e-2,
                               paraNorm=True, SBFA=False, Surr=False,
                               epsStop=1e-10, para=x, SensCalc="FD")
    f = rFuzz[0,0,0]/10000
    return -f, []

def SysEqIIunc(x, gc):
    from FuzzAnPy import FuzzAn
    from FuzzAnPy import FuzzyNumber
    nAlpha = 1
    nr = 18
    pFuzz = np.zeros([9, nAlpha, 2])
    pFuzz[0, :, :] = FuzzyNumber.Int(150e3, 250e3, nAlpha)
    pFuzz[1, :, :] = FuzzyNumber.Int(150e3, 250e3, nAlpha)
    pFuzz[2, :, :] = FuzzyNumber.Int(190e3, 210e3, nAlpha)
    pFuzz[3, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[4, :, :] = FuzzyNumber.Int(240e3, 260e3, nAlpha)
    pFuzz[5, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[6, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[7, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[8, :, :] = FuzzyNumber.Int(900e3, 1100e3, nAlpha)
    rFuzz, OutputData = FuzzAn(FuzzySysEq=FuzzySysEq, pFuzz=pFuzz, nr=nr,
                               Alg="MMA", nAlpha=nAlpha, deltax=1e-2,
                               paraNorm=True, SBFA=False, Surr=False,
                               epsStop=1e-10, para=x, SensCalc="FD")
    rFuzz *= 1/10000

#    fBand = np.array([[1200.0, 1800.0],
#                      [2400.0, 3600.0]])
#    fBand = np.array([[ 300.0,  450.0],
#                      [ 600.0,  900.0],
#                      [1200.0, 1800.0],
#                      [2400.0, 3600.0],
#                      [4800.0, 7200.0]])
    fBand = np.array([[  75.0,  100.0],
                      [ 150.0,  200.0],
                      [ 300.0,  400.0],
                      [ 600.0,  800.0],
                      [1200.0, 1600.0],
                      [2400.0, 3200.0]])
    g = np.zeros(np.shape(gc))
    ig = 0
    f0 = rFuzz
    nf0 = 18#len(f0)
    for ii in range(nf0):
        for jj in range(np.shape(fBand)[0]):
            g[ig] = np.prod((f0[ii,0,0]/fBand[jj,0]-1, 
                             fBand[jj,1]/f0[ii,0,0]-1))*(fBand[jj,1]+fBand[jj,0])/2.
            g[ig+1] = np.prod((f0[ii,0,1]/fBand[jj,0]-1,
                               fBand[jj,1]/f0[ii,0,1]-1))*(fBand[jj,1]+fBand[jj,0])/2.
            ig += 2
    f = x[0]+x[2]+x[4]+3*x[6]
    return f, g

def SysEqIIIunc(x, gc):
    from FuzzAnPy import FuzzAn
    from FuzzAnPy import FuzzyNumber
    nAlpha = 1
    nr = 18
    pFuzz = np.zeros([9, nAlpha, 2])
    pFuzz[0, :, :] = FuzzyNumber.Int(150e3, 250e3, nAlpha)
    pFuzz[1, :, :] = FuzzyNumber.Int(150e3, 250e3, nAlpha)
    pFuzz[2, :, :] = FuzzyNumber.Int(190e3, 210e3, nAlpha)
    pFuzz[3, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[4, :, :] = FuzzyNumber.Int(240e3, 260e3, nAlpha)
    pFuzz[5, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[6, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[7, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
    pFuzz[8, :, :] = FuzzyNumber.Int(900e3, 1100e3, nAlpha)
    rFuzz, OutputData = FuzzAn(FuzzySysEq=FuzzySysEq, pFuzz=pFuzz, nr=nr,
                               Alg="MMA", nAlpha=nAlpha, deltax=1e-2,
                               paraNorm=True, SBFA=False, Surr=False,
                               epsStop=1e-10, para=x, SensCalc="FD")
    rFuzz *= 1/10000

#    fBand = np.array([[1200.0, 1800.0],
#                      [2400.0, 3600.0]])
#    fBand = np.array([[ 300.0,  450.0],
#                      [ 600.0,  900.0],
#                      [1200.0, 1800.0],
#                      [2400.0, 3600.0],
#                      [4800.0, 7200.0]])
    fBand = np.array([[  75.0,  100.0],
                      [ 150.0,  200.0],
                      [ 300.0,  400.0],
                      [ 600.0,  800.0],
                      [1200.0, 1600.0],
                      [2400.0, 3200.0]])
    g = np.zeros(np.shape(gc))
    ig = 0
    f0 = rFuzz
    nf0 = 18#len(f0)
    for ii in range(nf0):
        for jj in range(np.shape(fBand)[0]):
            g[ig] = np.prod((f0[ii,0,0]/fBand[jj,0]-1, 
                             fBand[jj,1]/f0[ii,0,0]-1))*10*(fBand[jj,1]+fBand[jj,0])/2.
            g[ig+1] = np.prod((f0[ii,0,1]/fBand[jj,0]-1,
                               fBand[jj,1]/f0[ii,0,1]-1))*10*(fBand[jj,1]+fBand[jj,0])/2.
            ig += 2
    f = rFuzz[0,0,0]/10000
    return -f, g


OptForm = 5
Alg = "COBYLA"
Alg = "NLPQLP"
#Alg = "SLSQP"
#Alg = "NSGA2"
Alg = "MMA"
#Alg = "NSGA2"
#Alg = "GCMMA"
#Alg = "KSOPT"
#Alg = "ALPSO"
AlgOptions = OptAlgOptions.setDefault(Alg)
AlgOptions.setSimple(stopTol=1e-6)  
x0 = np.array([1.0e-3, 1e0,  5.0e-3, 5e0,  10.0e-3,  10e0, 0.5e-3, 5e-2])
xL = np.array([1.0e-4, 1e-1, 5.0e-4, 5e-1, 10.0e-4, 10e-1, 0.5e-4, 5e-3])
xU = np.array([1.0e-2, 1e1,  5.0e-2, 5e1,  10.0e-2,  10e1, 0.5e-2, 5e-1])
x0 = x0
from FuzzAnPy import FuzzAn
from FuzzAnPy import FuzzyNumber
nAlpha = 1
nr = 18
pFuzz = np.zeros([9, nAlpha, 2])
pFuzz[0, :, :] = FuzzyNumber.Int(150e3, 250e3, nAlpha)
pFuzz[1, :, :] = FuzzyNumber.Int(150e3, 250e3, nAlpha)
pFuzz[2, :, :] = FuzzyNumber.Int(190e3, 210e3, nAlpha)
pFuzz[3, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
pFuzz[4, :, :] = FuzzyNumber.Int(240e3, 260e3, nAlpha)
pFuzz[5, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
pFuzz[6, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
pFuzz[7, :, :] = FuzzyNumber.Int(4900e3, 5100e3, nAlpha)
pFuzz[8, :, :] = FuzzyNumber.Int(900e3, 1100e3, nAlpha)
if OptForm == 0:
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, \
        kRx, kRy, kRtheta, kSP, kCP, kRP, nP, alpha, rS, rC, rR, rP, \
        thetaCdot, thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    x=x0
    mS = x[0]
    JS = x[1]
    mC = x[2]
    JC = x[3]
    mR = x[4]
    JR = x[5]
    mP = x[6]
    JP = x[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                      kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                      kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                      kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                      rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                      thetaCdotdot=thetaCdotdot, Damping=True, All=0, Hz=1, 
                      Plot=1)  
    for ii in f0:
        print(ii)
    print(mS+ mC+ mR+3*mP)
    print(iEval)
elif OptForm == -1:
    rFuzzFD, OutputData = FuzzAn(FuzzySysEq=FuzzySysEq, pFuzz=pFuzz, nr=nr,
                                   Alg="NLPQLP", nAlpha=nAlpha, deltax=1e-2,
                                   paraNorm=True, SBFA=False, Surr=False,
                                   epsStop=1e-10, para=x0, SensCalc="FD")
    for ii in range(nr):
        print(rFuzzFD[ii]/10000)
    print(OutputData['nEval'])
    print(iEval)
    
elif OptForm == -10:
    from FuzzAnPy import FuzzAn
    from FuzzAnPy import FuzzyNumber
    rFuzz, OutputData = FuzzAn(FuzzySysEq=FuzzySysEq, pFuzz=pFuzz, nr=nr,
                                   Alg="MMA", nAlpha=nAlpha, deltax=1e-3,
                                   paraNorm=True, SBFA=False, Surr=False,
                                   epsStop=1e-6, para=x0, 
                                   FuzzySensEq=FuzzySensEq)
    for ii in range(nr):
        print(rFuzz[ii]/10000)
    print(OutputData['nEval'])
    print(iEval)
elif OptForm == 1:
  
    gc = []
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEqI, Alg=Alg,
                            StatusReport=True, DesVarNorm=True, DoE=False, 
                            SBDO=False, Debug=False, ResultReport=False,
                            deltax=1e-3, AlgOptions=AlgOptions)
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = xOpt[0]
    JS = xOpt[1]
    mC = xOpt[2]
    JC = xOpt[3]
    mR = xOpt[4]
    JR = xOpt[5]
    mP = xOpt[6]
    JP = xOpt[7]
    EE, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                      kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                      kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                      kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                      rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                      thetaCdotdot=0, Hz=True, Sum=True, All=True, 
                      Damping=True)
    for ii in EE:
        print(ii)
    print(SysEqI(xOpt, gc)[0])
    print(mS+ mC+ mR+3*mP)
    print(iEval)

elif OptForm == 10:
    gc = []
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEqI, Alg=Alg,
                            SensCalc="SensEq", SensEq=SensEqI,
                            StatusReport=True, DesVarNorm=True, DoE=False, 
                            SBDO=False, Debug=False, ResultReport=False,
                            deltax=1e-3, AlgOptions=AlgOptions)
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = xOpt[0]
    JS = xOpt[1]
    mC = xOpt[2]
    JC = xOpt[3]
    mR = xOpt[4]
    JR = xOpt[5]
    mP = xOpt[6]
    JP = xOpt[7]
    EE, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                      kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                      kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                      kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                      rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                      thetaCdotdot=0, Hz=True, Sum=True, All=False, 
                      Damping=True)
    for ii in EE:
        print(ii)
    print(SysEqI(xOpt, gc)[0])
    print(mS+ mC+ mR+3*mP)
    print(iEval)
elif OptForm == 2:
    gc = np.ones([18*6,])
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEqII, 
                            Alg=Alg, AlgOptions=AlgOptions, StatusReport=True, 
                            DesVarNorm=True,
                            DoE=False, SBDO=False, Debug=False, 
                            ResultReport=False, deltax=1e-3)
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = xOpt[0]
    JS = xOpt[1]
    mC = xOpt[2]
    JC = xOpt[3]
    mR = xOpt[4]
    JR = xOpt[5]
    mP = xOpt[6]
    JP = xOpt[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                        kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                        kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                        kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                        rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                        thetaCdotdot=thetaCdotdot)
    for ii in f0:
        print(ii)
    print(mS+ mC+ mR+3*mP)
    print(iEval)

elif OptForm == 20:
    gc = np.ones([18*6,])
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEqII, Alg=Alg,
                            SensCalc="SensEq", SensEq=SensEqII,
                            StatusReport=True, DesVarNorm=True, DoE=False, 
                            SBDO=False, Debug=False, ResultReport=False,
                            deltax=1e-3, AlgOptions=AlgOptions)
    mS = xOpt[0]
    JS = xOpt[1]
    mC = xOpt[2]
    JC = xOpt[3]
    mR = xOpt[4]
    JR = xOpt[5]
    mP = xOpt[6]
    JP = xOpt[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                        kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                        kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                        kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                        rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                        thetaCdotdot=thetaCdotdot)
    for ii in f0:
        print(ii)
    print(mS+ mC+ mR+3*mP)
    print(iEval)
        
elif OptForm==3:
    print("Formulation III Deterministic")
    gc = np.ones([18*6,])
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEqIII, 
                            Alg=Alg, StatusReport=True, DesVarNorm=True,
                            DoE=False, SBDO=False, Debug=False, 
                            ResultReport=False, deltax=1e-6,
                            AlgOptions=AlgOptions)
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = xOpt[0]
    JS = xOpt[1]
    mC = xOpt[2]
    JC = xOpt[3]
    mR = xOpt[4]
    JR = xOpt[5]
    mP = xOpt[6]
    JP = xOpt[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                        kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                        kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                        kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                        rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                        thetaCdotdot=thetaCdotdot)
    for ii in f0:
        print(ii)
    print(mS+ mC+ mR+3*mP)
    print(iEval)
elif OptForm == 4:
  
    gc = []
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEqIunc, Alg=Alg,
                            StatusReport=True, DesVarNorm=True, DoE=False, 
                            SBDO=False, Debug=False, ResultReport=False,
                            deltax=1e-3, AlgOptions=AlgOptions)
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = xOpt[0]
    JS = xOpt[1]
    mC = xOpt[2]
    JC = xOpt[3]
    mR = xOpt[4]
    JR = xOpt[5]
    mP = xOpt[6]
    JP = xOpt[7]
    EE, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                      kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                      kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                      kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                      rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                      thetaCdotdot=0, Hz=True, Sum=True, All=True, 
                      Damping=True)
    for ii in EE:
        print(ii)
    print(SysEqI(xOpt, gc)[0])
    print(mS+ mC+ mR+3*mP)
    print(iEval)
elif OptForm==5:
    gc = np.ones([18*6*2,])
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEqIIunc, 
                            Alg=Alg, StatusReport=True, DesVarNorm=True,
                            DoE=False, SBDO=False, Debug=False, 
                            ResultReport=False, deltax=1e-6, 
                            AlgOptions=AlgOptions)
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = xOpt[0]
    JS = xOpt[1]
    mC = xOpt[2]
    JC = xOpt[3]
    mR = xOpt[4]
    JR = xOpt[5]
    mP = xOpt[6]
    JP = xOpt[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                        kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                        kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                        kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                        rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                        thetaCdotdot=thetaCdotdot)
    for ii in f0:
        print(ii)
    print(mS+ mC+ mR+3*mP)
    print(iEval)
    
elif OptForm==6:
    gc = np.ones([18*6*2,])
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEqIIIunc, 
                            Alg=Alg, StatusReport=True, DesVarNorm=True,
                            DoE=False, SBDO=False, Debug=False, 
                            ResultReport=False, deltax=1e-6, 
                            AlgOptions=AlgOptions)
    mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, kRx,\
        kRy, kRtheta, kSP, kCP, kRP, nP, alpha,  rS, rC, rR, rP, thetaCdot,\
        thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)
    mS = xOpt[0]
    JS = xOpt[1]
    mC = xOpt[2]
    JC = xOpt[3]
    mR = xOpt[4]
    JR = xOpt[5]
    mP = xOpt[6]
    JP = xOpt[7]
    f0, Phi = SysEval(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                        kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                        kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                        kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                        rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                        thetaCdotdot=thetaCdotdot)
    for ii in f0:
        print(ii)
    print(mS+ mC+ mR+3*mP)
    print(iEval)
"""
if symbolic:
    print M
    print D
    print K
else:
    print((K.transpose() == K).all())
    print(np.allclose(K, K.T, atol=1e-4))
    #import LumPy
    import Solver as OldSolver
    from LumPy import Solver
    omega, Phi = Solver.ModalDamped(M, D, K, Hz=True, All=False, Sum=True)
    for ii in omega:
        print ii
    print("\n")
"""
#    omega, Phi = Solver.ModalUndamped(M, K, Hz=True, All=True, Sum=True)
#    for ii in omega:
#        print ii
#    print("yo")
#    omega, Phi = OldSolver.ModalDamped(M, D, K, Hz=True, All=False)
#    for ii in omega:
#        print ii
    
#    omega, Phi = Solver.ModalDamped(M, D, K, Hz=True, All=True, Sum=True)
#    for ii in omega:
#        print ii
"""
from GearAnalysis import SetGearSet, GearEval
from mass import mass
from stiffness import stiffness
from damping import damping
x0 = np.array([0.262, 21.0, 0.025,  9.11, 0.44, 189.0, 2.0, 3.27e3])
nP, Rot, nf, phi, alpha, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry, rs, rp, rpt, rr, Damping = SetGearSet()
ms = x0[0]
Js = x0[1]
mp = x0[2]
Jp = x0[3]
mc = x0[4]
Jc = x0[5]
mr = x0[6]
Jr = x0[7]
Rot = 200
M1 = mass(nf, ms, Js, mp, Jp, mc, Jc, mr, Jr,rpt)
K1 = stiffness(nf, Rot, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry,
                  mp, rs, rp, rpt, rr, alpha, phi)
D1 = damping(nf,Rot, mp)
"""