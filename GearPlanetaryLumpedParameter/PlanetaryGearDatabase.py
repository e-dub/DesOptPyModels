# -*- coding: utf-8 -*-
"""
Title:              PlanetaryGearDatabase.py
Units:              N, mm, t, MPa, s
Author:             E.J. Wehrle
Date:               May 4, 2017
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Planetary gear database
-------------------------------------------------------------------------------
"""
from __future__ import absolute_import, division, print_function
import numpy as np

def GetParam(GearSet="Benchmark"):
    if GearSet=="TQ070i3":
        rS=32.4/2
        rC=24.5
        rR=64.8/2
        rP=15.6/2
        mS = 0.262e-3
        JS = 21e-5
        mC = 0.725e-3
        JC = 189e-6
        mR = 1.87e-3
        JR = 3131e-7
        mP = 0.025e-3
        JP = 9.11e-5
        kSx = 213e3
        kSy = 213e3
        kStheta = 148e3
        kCx = 270e3
        kCy = 270e3
        kCtheta = 125e3
        kRx = 6e10
        kRy = 6e10
        kRtheta = 1e10
        kSP = 213e3
        kCP = 1250e3
        kRP = 207e3
        nP = 3.
        alpha = 20.
        thetaCdot = 1000./60.*2*np.pi
        thetaCdotdot = 0
    elif GearSet=="MP080i6":
        rS = 32.4/2.
        rC = 24.5
        rR = 64.8/2
        rP = 15.6/2
        mS = 0.262e-3
        JS = 21e-3
        mC = 0.44e-3
        JC = 189e-3
        mR = 2e-3
        JR = 3.27
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
        #thetaCdot = 0
        thetaCdotdot=0
    elif GearSet=="Me":
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
    elif GearSet=="TurboProp":
        rS = 55.75/2.
        rC = 88.6/2
        rR = 109.7/2.
        rP = 27./2.
        mS = 3e-3
        JS = 1.75e-3
        mC = 12e-3
        JC = 6.8e-3
        mR = 7.64e-3
        JR = 8.09e-3
        mP = 1.86e-3
        JP = 1.25e-3
        kSx = 20e3
        kSy = 20e3
        kStheta = 20e3
        kCx = 50e3
        kCy = 50e3
        kCtheta = 100e3
        kRx = 100e3
        kRy = 100e3
        kRtheta = 100e3
        kSP = 100e3
        kCP = 10e3
        kRP = 100e3
        nP=3.
        alpha=20.
        thetaCdot=20000./60.*2*np.pi
        thetaCdot=0
        thetaCdotdot=0
    elif GearSet=="HeliOH58":    
        rS = 77.42/2.
        rC = 177.8/2.
        rR = 275.03/2.
        rP = 100.35/2.
        mS = 0.4e-3
        JS = 0.39e-3*rS**2
        mC = 5.43e-3
        JC = 6.29e-3*rC**2
        mR = 2.35e-3
        JR = 3.00e-3*rR**2
        mP = 0.66e-3
        JP = 0.61e-3*rP**2
        kSx = 100e3
        kSy = 100e3
        kStheta = 1000e3
        kStheta = 0
        kCx = 100e3
        kCy = 100e3
        kCtheta = 1000e3
        kCtheta = 0
        kRx = 100e3
        kRy = 100e3
        kRtheta = 1000e3    
        kSP = 500e3
        kCP = 100e3
        kRP = 500e3
        nP=3.
        alpha=24.6
        thetaCdot=1000./60.*2*np.pi
        thetaCdotdot=0.
    elif GearSet=="Benchmark":
        rS = 50.
        rC = 75.
        rR = 100.
        rP = 25.
        mS = 1.0e-3
        JS = 1.0
        mC = 5.0e-3
        JC = 5.0
        mR = 10.0e-3
        JR = 10.0
        mP = 0.5e-3
        JP = 0.05
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
        thetaCdot = 523.5987
        #thetaCdot = 0
        thetaCdotdot=0
    elif GearSet=="EricsonParker2013A":
        mS=7.657
        mP=1.724
        mR=110.1
        mC=28.87
        JS = 29.89e-3
        JP = 2.378e-3
        JR = 0.0
        JC = 346.3e-3
        
        rS=71.25/2*1e-3
        rR=262.7/2*1e-3
        rP=91.29/2*1e-3
        rC=rS+rP

        kSx = 98e6
        kSy = 98e6
        kStheta = 316e3
        kCx = 10.1e6
        kCy = 10.1e6
        kCtheta = 231e3
        kRx = 1.72e9
        kRy = 10.2e9
        kRtheta = 1e10
        kSP = 89.1e6
        kCP = 315e6
        kRP = 120e6
        nP = 5.
        alpha = 25.
        thetaCdot = 1000./60.*2*np.pi
        thetaCdot = 0
        thetaCdotdot = 0
    elif GearSet=="EricsonParker2013B":
        mS = 12.44
        mP = 1.433
        mR = 110.2
        mC = 28.38
        JS = 64.04e-3
        JP = 1.433e-3
        JR = 0.0
        JC = 398.7e-3
        
        rS = 71.25/2*1e-3
        rR = 262.7/2*1e-3
        rP = 91.29/2*1e-3
        rC = rS+rP

        kSx = 98e6
        kSy = 98e6
        kStheta = 316e3
        kCx = 10.1e6
        kCy = 10.1e6
        kCtheta = 231e3
        kRx = 1.72e9
        kRy = 10.2e9
        kRtheta = 0
        kSP = 89.1e6
        kCP = (146e6+427e6)/2.
        kRP = 120e6
        nP = 5.
        alpha = 25.
        thetaCdot = 1000./60.*2*np.pi
        thetaCdot = 0
        thetaCdotdot = 0
    return(mS, JS, mC, JC, mR, JR, mP, JP, kSx, kSy, kStheta, kCx, kCy, 
           kCtheta, kRx, kRy, kRtheta, kSP, kCP, kRP, nP, alpha, rS, rC, rR, 
           rP, thetaCdot, thetaCdotdot)   
if __name__ == "__main__":
    mS, JS, mC, JC, mR, JR, mP, JP, \
               kSx, kSy, kStheta,   \
               kCx, kCy, kCtheta,   \
               kRx, kRy, kRtheta,   \
               kSP, kCP, kRP,       \
               nP, alpha,           \
               rS, rC, rR, rP,      \
               thetaCdot, thetaCdotdot = GetParam()
    print(mS*rS**2/2.)
    print(JS)
    print(mR*rR**2/2.)
    print(JR)
    print(mC*rC**2/2.)
    print(JC)
    print(mP*rP**2/2.)
    print(JP)
    print(3.14*rS**2*7.8e-9*100)
    print(mS)