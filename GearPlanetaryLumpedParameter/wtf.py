#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 09:00:42 2017

@author: wehrle
"""
from Solver import ModalSensUndampedHaftka1, rad2hz, ModalSensDamped1
from LumPy import Solver, PlanetaryGearDatabase
#import 
from LumPy.PlanetaryGear import CalcMatrices, PlotGearSet
import numpy as np


GearSet = "Benchmark"
Damping = True
All = True
Sum = True
Hz=True
mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, kCtheta, \
    kRx, kRy, kRtheta, kSP, kCP, kRP, nP, alpha, rS, rC, rR, rP, \
    thetaCdot, thetaCdotdot = PlanetaryGearDatabase.GetParam(GearSet)

M, D, K, dof = CalcMatrices(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, 
                            JP=JP, kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, 
                            kCy=kCy, kCtheta=kCtheta, kRx=kRx, kRy=kRy, 
                            kRtheta=kRtheta, kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, 
                            alpha=alpha, rS=rS, rC=rC, rR=rR, rP=rP, 
                            thetaCdot=thetaCdot, thetaCdotdot=thetaCdotdot)
if Damping:
    omega, Phi = Solver.ModalDamped(M=M, D=D, K=K, Hz=Hz, All=All, Sum=Sum)    
else:
    omega, Phi = Solver.ModalUndamped(M=M, D=D, K=K,  Hz=Hz, All=All, 
                                      Sum=Sum) 
print omega