# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 09:29:58 2016

@author: wehrle
"""
from __future__ import division
from numpy import(append, pi, max, zeros, array, ceil, sort, concatenate, ones,
                  abs, shape, min, prod, eye, set_printoptions)
import scipy.linalg as linalg
from Solver import ModalDamped
#from Solver import ModalDampedI
from Solver import ModalUndamped
#from Solver import ModalUndampedI
from Solver import ModalSensDamped
from Solver import ModalSensUndamped
#import numpy as np
#def GearPlantetaryAnalysis():
from mass import mass
from damping import damping
from stiffness import stiffness
import sys
import copy
import sys
#sys.path.append("/home/wehrle/opt/")
from FuzzAnPy import FuzzAn
from FuzzAnPy import FuzzyNumber
from termcolor import colored

#sys.path.append('/home/wehrle/opt/StAnPy')
#from StructureAnalyzer import StructureModalResponse
delta = 0.001
def SetGearSet():
    Damping=1
    # Concept definition
    nP = 3          # Number of planets
    # Velocity
    Rot = 3000/3    #velocita media del portatreno [rpm]
    # Geometry definition
    rs = 32.4/2     # radius of sun [mm]
    rp = 15.6/2     # radius of planet gears
    rpt = 24.3      # radius of planet carrier
    rr = 64.8/2     # radius of ring gear
    # Stiffness definition (later with finite element analysis)
    ks = 170e3      #[N/mm]
    kp = 1250e3     # prova
    kpt = 2*135e3   #170000
    kr = 170e3      #
    kTs = 4700e3    # prova
    kTr = 10e9      # prova
    kTpt = 32000e3  # prova
    ksx = 170e3     # 165000#170000 # prova
    ksy = ksx       #
    krx = 65e9      # prova
    kry = krx       #
    Rot = Rot/60*2*pi   # Angular velocity in rad/s
    nf = 3*(nP+3)       #Degrees of freedom
    phi = zeros(nP,)
    for ii in range(nP-1):
        phi[ii+1] = 360.0/nP*(1+ii)/180.0*pi
    alpha = 20/180*pi
    return(nP, Rot, nf, phi, alpha, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry, rs,
           rp, rpt, rr, Damping)

def CompSens(SysEq, SensEq, x, gc, delta):
    fIt, gIt = SysEq(x, gc)
    dfdxFD = zeros(shape(x))
    dgdxFD = zeros((len(x),len(gIt)))
    for ii in range(len(x)):
        xFD = copy.deepcopy(x)
        xDelta = delta #delta*x[ii]
        xFD += eye(len(x))[:,ii]*xDelta
        f, g = SysEq(xFD, gc)
        dfdxFD[ii] = (f-fIt)/xDelta
        if len(g)>0:
            dgdxFD[ii] = (g-gIt)/xDelta
        else:
            dgdxFD = []
    dfdxAn, dgdxAn = SensEq(x, fIt, gIt, gc)
    return dfdxFD, dgdxFD, dfdxAn, dgdxAn


def GearEval(nP, Rot, nf, phi, alpha, rs, rp, rpt, rr, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx,
             ksy, krx, kry, ms, Js, mp, Jp, mc, Jc, mr, Jr, Damping=True):
    M = mass(nf, ms, Js, mp, Jp, mc, Jc, mr, Jr,rpt)
    K = stiffness(nf, Rot, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry,
                  mp, rs, rp, rpt, rr, alpha, phi)
    D = damping(nf,Rot, mp)
    if Damping:
        f0, Phi = ModalDamped(M, D, K, Hz=False, All=False)
    else:
        f0, Phi = ModalUndamped(M, K, Hz=False, All=True)
    return f0, Phi

def SensEq(x, f, g, gc):
    ms = x[0]
    Js = x[1]
    mp = x[2]
    Jp = x[3]
    mc = x[4]
    Jc = x[5]
    mr = x[6]
    Jr = x[7]
    nP, Rot, nf, phi, alpha, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry, rs, rp, rpt, rr, Damping = SetGearSet()
    M0 = mass(nf, ms, Js, mp, Jp, mc, Jc, mr, Jr, rpt)
    K0 = stiffness(nf, Rot, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy,       krx, kry,
                mp, rs, rp, rpt, rr, alpha, phi)
    D0 = damping(nf,Rot, mp)

    if Damping:
        f0, Phi = ModalDamped(M0, D0, K0, Hz=False, All=False)
    else:
        f0, Phi = ModalUndamped(M0, K0, Hz=False, All=True)
    domegadx = zeros(shape(x))
    for ii in range(len(x)):
        xFD = copy.deepcopy(x)
        xDelta = delta #xFD[ii]*delta
        xFD += eye(len(x))[:,ii]*xDelta
        set_printoptions(precision=4)
        #print(xFD)
        msFD = xFD[0]
        JsFD = xFD[1]
        mpFD = xFD[2]
        JpFD = xFD[3]
        mcFD = xFD[4]
        JcFD = xFD[5]
        mrFD = xFD[6]
        JrFD = xFD[7]
        M = mass(nf, msFD, JsFD, mpFD, JpFD, mcFD, JcFD, mrFD, JrFD, rpt)
        K = stiffness(nf, Rot, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx,
                      kry, mp, rs, rp, rpt, rr, alpha, phi)
        D = damping(nf,Rot, mp)
        f0_ = f0[0]
        Phi_ = Phi[:,0]
        dMdx = (M-M0)/xDelta
        dKdx = (K-K0)/xDelta
        dDdx = (D-D0)/xDelta
        if Damping:
            domegadx[ii] = ModalSensDamped(f0_, Phi_, M, D, K, dMdx, dDdx, dKdx)
        else:
            domegadx[ii] = ModalSensUndamped(f0_, Phi_, M, K, dMdx, dKdx)
    dgdx = []
    return -domegadx, dgdx


def SysEq1(x, gc):
    nP, Rot, nf, phi, alpha, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry, rs, rp, rpt, rr, Damping = SetGearSet()
    ms = x[0]
    Js = x[1]
    mp = x[2]
    Jp = x[3]
    mc = x[4]
    Jc = x[5]
    mr = x[6]
    Jr = x[7]
    EE, Phi = GearEval(nP, Rot, nf, phi, alpha, rs, rp, rpt, rr, ks, kp, kpt, kr, kTs, kTr,
                       kTpt, ksx, ksy, krx, kry, ms, Js, mp, Jp, mc, Jc, mr,
                       Jr, Damping)
    g = []
    f = EE[0]
    return -f, g

def SysEq2(x, gc):
    nP, Rot, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry, rs, rp, rpt, rr, Damping = SetGearSet()
    ms = x[0]
    Js = x[1]
    mp = x[2]
    Jp = x[3]
    mc = x[4]
    Jc = x[5]
    mr = x[6]
    Jr = x[7]
    EE, Phi = GearEval(nP, Rot, rs, rp, rpt, rr, ks, kp, kpt, kr, kTs, kTr,
                       kTpt, ksx, ksy, krx, kry, ms, Js, mp, Jp, mc, Jc, mr, Jr)
    #1st harmonic
    omegal = 2000.0/60.0*108.0/3.0
    omegaU = 2900.0/60.0*108.0/3.0
    fBand = array([[1200.0, 1800.0],
                   [2400.0, 3600.0]])
    f0 = array(sorted(EE))
    g = zeros(shape(gc))
    ig = 0
    nf0 = 36
    for ii in range(nf0):
        for jj in range(shape(fBand)[0]):
            g[ig] = prod((f0[ii]/fBand[jj,0]-1, fBand[jj,1]/f0[ii]-1))*10000.
            ig += 1
    f = x[0]+3*x[2]+x[4]+x[6]
    return f, g

def FuzzySysEq(p, x, ir):
    nP, Rot, nf, phi, alpha, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry, rs, rp, rpt, rr, Damping = SetGearSet()
    ms = x[0]
    Js = x[1]
    mp = x[2]
    Jp = x[3]
    mc = x[4]
    Jc = x[5]
    mr = x[6]
    Jr = x[7]

    ks = p[0]
    kp = p[1]
    kpt = p[2]
    kr = p[3]
    kTs = p[4]
    kTr = p[5]
    kTpt = p[6]
    ksx = ksy = p[7]
    krx = kry = p[8]
    Rot = p[9]/3
    EE, Phi = GearEval(nP, Rot, nf, phi, alpha, rs, rp, rpt, rr, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx,
                           ksy, krx, kry, ms, Js, mp, Jp, mc, Jc, mr, Jr, Damping)
    return(EE[ir])




from DesOptPy import DesOpt, OptAlgOptions
Alg="COBYLA"
Alg="SLSQP"
#Alg = "SLSQP"
AlgOptions = OptAlgOptions.setDefault(Alg)
AlgOptions.setSimple(stopTol=1e-6)
#Alg = "SLSQP"
#Alg = "NSGA2"
#Alg = "MMA"
#Alg = "NSGA2"
#Alg = "GCMMA"
#Alg = "KSOPT"
#Alg = "ALPSO"
x0 = array([0.262, 21.0, 0.025,  9.11, 0.44, 189.0, 2.0, 3.27e3])
xL = array([0.150, 10.0, 0.015,  6.00, 0.30, 170.0, 1.0, 2.50e3])
xU = array([0.400, 30.0, 0.035, 12.00, 0.60, 210.0, 3.0, 4.00e3])
OptForm = 0
if OptForm == 0:
    nP, Rot, nf, phi, alpha, ks, kp, kpt, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry, rs, rp, rpt, rr, Damping = SetGearSet()
    ms = x0[0]
    Js = x0[1]
    mp = x0[2]
    Jp = x0[3]
    mc = x0[4]
    Jc = x0[5]
    mr = x0[6]
    Jr = x0[7]
    EE, Phi = GearEval(nP, Rot, nf, phi, alpha, rs, rp, rpt, rr, ks, kp, kpt, kr, kTs, kTr,
                       kTpt, ksx, ksy, krx, kry, ms, Js, mp, Jp, mc, Jc, mr,
                       Jr, Damping)
    for ii in EE:
        print ii
elif OptForm == 1:
    gc = []
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq1, Alg=Alg,
                        StatusReport=True, DesVarNorm=True, DoE=False, SBDO=False,
                        Debug=False, ResultReport=True, deltax=1e-3)

elif OptForm == 11:
    gc = []
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq1, SensEq=SensEq, Alg=Alg,
                        StatusReport=True, DesVarNorm=True, DoE=False, SBDO=False,
                        Debug=False, ResultReport=True, deltax=1e-3)
elif OptForm == 2:
    gc = ones([36*2,])
    xOpt, fOpt, SP = DesOpt(x0=x0, xL=xL, xU=xU, gc=gc, SysEq=SysEq2, Alg=Alg,
                        StatusReport=False, DesVarNorm=True, DoE=False, SBDO=False,
                        Debug=False, ResultReport=False, deltax=1e-3)

elif OptForm == 3:
    dfdxFD, dgdxFD, dfdxAn, dgdxAn = CompSens(SysEq1, SensEq, x0, gc=[], delta=0.001)
    print dfdxFD
    print dfdxAn
elif OptForm == 4:
    nAlpha = 1
    #    ks = 170e3      #[N/mm]
    #    kp = 1250e3     # prova
    #    kpt = 2*135e3   #170000
    #    kr = 170e3      #
    #    kTs = 4700e3    # prova
    #    kTr = 10e9      # prova
    #    kTpt = 32000e3  # prova
    #    ksx = 170e3     # 165000#170000 # prova
    #    ksy = ksx       #
    #    krx = 65e9      # prova
    #    kry = krx       #
    import numpy as np
    pFuzz = np.zeros([10, nAlpha, 2])
    pFuzz[0, :, :] = FuzzyNumber.Int(169e3, 171e3, nAlpha)
    pFuzz[1, :, :] = FuzzyNumber.Int(1240e3, 1255e3, nAlpha)
    pFuzz[2, :, :] = FuzzyNumber.Int(2*134e3, 2*136e3, nAlpha)
    pFuzz[3, :, :] = FuzzyNumber.Int(169e3, 171e3, nAlpha)
    pFuzz[4, :, :] = FuzzyNumber.Int(4699e3, 4701e3, nAlpha)
    pFuzz[5, :, :] = FuzzyNumber.Int(9.5e9, 10.5e9, nAlpha)
    pFuzz[6, :, :] = FuzzyNumber.Int(31900e3, 32100e3, nAlpha)
    pFuzz[7, :, :] = FuzzyNumber.Int(169e3, 171e3, nAlpha)
    pFuzz[8, :, :] = FuzzyNumber.Int(68e9, 72e9, nAlpha)
    pFuzz[9, :, :] = FuzzyNumber.Int(1500, 3500, nAlpha)
    rFuzz, OutputData = FuzzAn(FuzzySysEq=FuzzySysEq,
                               pFuzz=pFuzz, nr=18, para=x0,
                               Alg="NLPQLP", nAlpha=nAlpha)
    print rFuzz
    print OutputData["nEval"]
#ms = xOpt[0]
#Js = xOpt[1]
#mp = xOpt[2]
#Jp = xOpt[3]
#mc = xOpt[4]
#Jc = xOpt[5]
#mr = xOpt[6]
#Jr = xOpt[7]
##f, g = SysEq2(xOpt, gc)
#f0Opt, Phi = GearEval(ms, Js, mp, Jp, mc, Jc, mr, Jr)
#print f0Opt



#if __name__ == "__main__":
#
#    ms = 0.262
#    Js = 21
#    mp = 0.025
#    Jp = 9.11
#    mc = 0.44
#    Jc = 189
#    mr = 2
#    Jr = 3.27e3
#    EE = GearEval(ms, Js, mp, Jp, mc, Jc, mr, Jr)
#    print EE
#
