# -*- coding: utf-8 -*-
"""
Title:              PlanetGear.py
Units:              N, mm, t, MPa, s
Author:             E.J. Wehrle
Date:               May 4, 2017
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Description:

Kratos model "taken apart"
 Questions etc. for understanding, further developement.

-------------------------------------------------------------------------------
"""
from __future__ import division
import numpy as np
import sys

def CalcMatrices(mS=0., JS=0., mC=0., JC=0., mR=0., JR=0., mP=0., JP=0., 
                 kSx=0., kSy=0., kStheta=0., kCx=0., kCy=0., kCtheta=0., 
                 kRx=0., kRy=0., kRtheta=0., kSP=0., kCP=0., kRP=0, 
                 nP=3., alpha=0., rS=0., rC=0., rR=0., rP=0., 
                 thetaCdot=0., thetaCdotdot=0., symbolic=False, debug=False):
    if symbolic:
        import sympy as sp
        mS, JS, mC, JC, mR, JR, mP, JP = sp.symbols('m_S J_S m_C J_C m_R J_R m_P J_P')
        kSx, kSy, kStheta = sp.symbols('k_S_x, k_S_y, k_S_theta')
        kCx, kCy, kCtheta = sp.symbols('k_C_x k_C_y, k_C_theta')
        kRx, kRy, kRtheta = sp.symbols('k_R_x, k_R_y, k_R_theta')
        kSP, kCP, kRP = sp.symbols(' k_SP, k_CP, k_RP')
        alpha, rS, rC, rR, rP = sp.symbols('alpha r_S r_C r_R r_P')
        thetaCdot, thetaCdotdot= sp.symbols('theta_C_dot theta_C_dotdot')
        from sympy import cos, sin 
    else:
        from numpy import cos, sin, deg2rad
        alpha = deg2rad(alpha)
def GetTransMat():
    QS = np.zeros([3, 9+int(nP)*3])         # Transformation matrix sun gear
    QS[0, 0] = 1.
    QS[1, 1] = 1.
    QS[2, 2] = 1.
    QC = np.zeros([3, 9+int(nP)*3])         # Transformation matrix carrier
    QC[0, 3] = 1.
    QC[1, 4] = 1.
    QC[2, 5] = 1.
    QR = np.zeros([3, 9+int(nP)*3])         # Transformation matrix ring gear
    QR[0, 6] = 1.
    QR[1, 7] = 1.
    QR[2, 8] = 1.
    MS = np.array([[mS, 0., 0.],            # Mass matrix sun gear
                   [0., mS, 0.],
                   [0., 0., JS]])
    MC = np.array([[mC, 0., 0.],            # Mass matrix carrier
                   [0., mC, 0.],
                   [0., 0., JC]])
    MR = np.array([[mR, 0., 0.],            # Mass matrix ring gear
                   [0., mR, 0.],
                   [0., 0., JR]])
    MP = {}
    QP = {}
    for ii in range(int(nP)):
        MP[ii] = np.array([[mP, 0., 0.],
                           [0., mP, 0.],
                           [0., 0., JP]])
        QP[ii] = np.zeros([3, 9+int(nP)*3])
        QP[ii][0, 9+3*ii+0] = 1
        QP[ii][1, 9+3*ii+1] = 1
        QP[ii][2, 9+3*ii+2] = 1
    M = QS.T.dot(MS.dot(QS)) + \
        QC.T.dot(MC.dot(QC)) + \
        QR.T.dot(MR.dot(QR))
    for ii in range(int(nP)):
        M += QP[ii].T.dot(MP[ii].dot(QP[ii]))
    # Gyroscopic damping matrices
    DGyroP = {}
    for ii in range(int(nP)):
        DGyroP[ii] = 2*mP*thetaCdot*np.array([[0., -1., 0.],
                                              [1.,  0., 0.],
                                              [0.,  0., 0.]])
    D = np.zeros([9+int(nP)*3, 9+int(nP)*3])
    for ii in range(int(nP)):
        D = D + QP[ii].T.dot(DGyroP[ii].dot(QP[ii]))
    # Stiffness from fixed central bearing
    KS = np.array([[kSx,  0.,      0.],
                   [ 0., kSy,      0.],
                   [ 0.,  0., kStheta]])
    KC = np.array([[kCx,  0.,      0.],
                   [ 0., kCy,      0.],
                   [ 0.,  .0, kCtheta]])
    KR = np.array([[kRx,  0.,      0.],
                   [ 0., kRy,      0.],
                   [ 0.,  0., kRtheta]])
    KP = {}
    for ii in range(int(nP)):
        KP[ii] = np.array([[kCP, 0., 0.],
                           [0., kCP, 0.],
                           [0., 0.,  0.]])
    KSP = {}
    QSP = {}
    KCP = {}
    QCP = {}
    KRP = {}
    QRP = {}
    KGyroP = {}
    KDragP = {}
    for ii in range(int(nP)):
        phi = 2*np.pi/nP*ii 
        # Sun-planet stiffnesses  
        KSP[ii] = kSP*np.array([[ sin(phi-alpha)*sin(phi-alpha), -sin(phi-alpha)*cos(phi-alpha), -sin(phi-alpha)*rS,  sin(phi-alpha)*sin(alpha),  sin(phi-alpha)*cos(alpha), -sin(phi-alpha)*rP],
                                [-cos(phi-alpha)*sin(phi-alpha),  cos(phi-alpha)*cos(phi-alpha),  cos(phi-alpha)*rS, -cos(phi-alpha)*sin(alpha), -cos(phi-alpha)*cos(alpha),  cos(phi-alpha)*rP],
                                [    -rS*sin(phi-alpha),                      rS*cos(phi-alpha),              rS*rS,             -rS*sin(alpha),             -rS*cos(alpha),              rS*rP],
                                [ sin(alpha)*sin(phi-alpha),         -sin(alpha)*cos(phi-alpha),     -sin(alpha)*rS,      sin(alpha)*sin(alpha),      sin(alpha)*cos(alpha),     -sin(alpha)*rP],
                                [ cos(alpha)*sin(phi-alpha),         -cos(alpha)*cos(phi-alpha),     -cos(alpha)*rS,      cos(alpha)*sin(alpha),      cos(alpha)*cos(alpha),     -cos(alpha)*rP],
                                [    -rP*sin(phi-alpha),                      rP*cos(phi-alpha),              rP*rS,             -rP*sin(alpha),             -rP*cos(alpha),              rP*rP]])   
        QSP[ii] = np.zeros([6, 9+int(nP)*3])
        QSP[ii][0, 0] = 1.
        QSP[ii][1, 1] = 1.
        QSP[ii][2, 2] = 1.
        QSP[ii][3, 9+3*ii+0] = 1.
        QSP[ii][4, 9+3*ii+1] = 1.
        QSP[ii][5, 9+3*ii+2] = 1.
        print(KSP[ii][0,0])
        # Carrier-planet stiffnesses
        KCP[ii] = kCP*np.array([[         1,         0, -rC*sin(phi), -cos(phi),  sin(phi), 0.],
                                [         0,         1,  rC*cos(phi), -sin(phi), -cos(phi), 0.],
                                [-rC*sin(phi), rC*cos(phi),      rC*rC,      0.,     -rC, 0.],
                                [   -cos(phi),   -sin(phi),         0.,      1.,      0., 0.],
                                [    sin(phi),   -cos(phi),        -rC,      0.,      1., 0.],
                                [        0.,        0.,         0.,      0.,      0., 0.]])
        QCP[ii] = np.zeros([6, 9+int(nP)*3])
        QCP[ii][0, 3] = 1.
        QCP[ii][1, 4] = 1.
        QCP[ii][2, 5] = 1.
        QCP[ii][3, 9+3*ii+0] = 1.
        QCP[ii][4, 9+3*ii+1] = 1.
        QCP[ii][5, 9+3*ii+2] = 1.
        KRP[ii] = kRP*np.array([[ sin(phi+alpha)*sin(phi+alpha), -sin(phi+alpha)*cos(phi+alpha), -sin(phi+alpha)*rR, -sin(phi+alpha)*sin(alpha),  sin(phi+alpha)*cos(alpha),  sin(phi+alpha)*rP],
                                [-cos(phi+alpha)*sin(phi+alpha),  cos(phi+alpha)*cos(phi+alpha),  cos(phi+alpha)*rR,  cos(phi+alpha)*sin(alpha), -cos(phi+alpha)*cos(alpha), -cos(phi+alpha)*rP],
                                [    -rR*sin(phi+alpha),      rR*cos(phi+alpha),      rR*rR,      rR*sin(alpha),     -rR*cos(alpha),     -rR*rP],
                                [-sin(alpha)*sin(phi+alpha),  sin(alpha)*cos(phi+alpha),  sin(alpha)*rR,  sin(alpha)*sin(alpha), -sin(alpha)*cos(alpha), -sin(alpha)*rP],
                                [ cos(alpha)*sin(phi+alpha), -cos(alpha)*cos(phi+alpha), -cos(alpha)*rR, -cos(alpha)*sin(alpha),  cos(alpha)*cos(alpha),  cos(alpha)*rP],
                                [     rP*sin(phi+alpha),     -rP*cos(phi+alpha),     -rP*rR,     -rP*sin(alpha),      rP*cos(alpha),      rP*rP]])
        QRP[ii] = np.zeros([6, 9+int(nP)*3])
        QRP[ii][0, 6] = 1.
        QRP[ii][1, 7] = 1.
        QRP[ii][2, 8] = 1.
        QRP[ii][3, 9+3*ii+0] = 1.
        QRP[ii][4, 9+3*ii+1] = 1.
        QRP[ii][5, 9+3*ii+2] = 1.
        # Gyroscopic stiffness matrices
        KGyroP[ii] = -mP*thetaCdot**2*np.array([[1., 0., 0.],
                                                [0., 1., 0.],
                                                [0., 0., 0.]])
        # Dynamic drag stiffness matrices
        KDragP[ii] = mP*thetaCdotdot*np.array([[0., 1., 0.],
                                               [1., 0., 0.],
                                               [0., 0., 0.]])
    All={}
    IncludeGyro = 1
    IncludeDrag = 0
    K = QS.T.dot(KS.dot(QS)) + \
        QC.T.dot(KC.dot(QC)) + \
        QR.T.dot(KR.dot(QR))
    All["KS"] = QS.T.dot(KS.dot(QS))
    All["KC"] = QC.T.dot(KC.dot(QC))
    All["KR"] = QR.T.dot(KR.dot(QR))
    for ii in range(int(nP)):
        K += QSP[ii].T.dot(KSP[ii].dot(QSP[ii]))
        All["KSP"+str(ii)] = QSP[ii].T.dot(KSP[ii].dot(QSP[ii]))
    for ii in range(int(nP)):
        K += QCP[ii].T.dot(KCP[ii].dot(QCP[ii]))
        All["KCP"+str(ii)] = QCP[ii].T.dot(KCP[ii].dot(QCP[ii]))
    for ii in range(int(nP)):
        K += QRP[ii].T.dot(KRP[ii].dot(QRP[ii]))
        All["KRP"+str(ii)] = QRP[ii].T.dot(KRP[ii].dot(QRP[ii]))
    for ii in range(int(nP)):
        if IncludeDrag:
            K += QP[ii].T.dot(KDragP[ii].dot(QP[ii]))
        All["KDragP"+str(ii)] = QP[ii].T.dot(KDragP[ii].dot(QP[ii]))
        if IncludeGyro:
            K += QP[ii].T.dot(KGyroP[ii].dot(QP[ii]))
        All["KGyroP"+str(ii)] = QP[ii].T.dot(KGyroP[ii].dot(QP[ii]))
    return(M, D, K)


if __name__ == "__main__":
    print("Testing of the calculation and assembly of the system matricies\n")
    try: 
        GearSet = sys.argv[1]    
    except: 
        GearSet = "TQ070i3"
    if GearSet=="TQ070i3":
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
        kRx = 6e10
        kRy = 6e10
        kRtheta = 1e10
        kSP = 213e3
        kCP = 1250e3
        kRP = 207e3
        nP = 3.
        alpha = 20.
        thetaCdot = 0.
        thetaCdotdot = 0.
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
        thetaCdot=0.
        thetaCdotdot=0.
    M, D, K = CalcMatrices(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                           kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                           kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                           kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                           rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                           thetaCdotdot=thetaCdotdot, symbolic=False)
    print("K is symmetric: " + str((K.transpose() == K).all()))
    print("M is symmetric:"  + str((M.transpose() == M).all()))
    print("D is skewed symmetric: "+ str((D.transpose() == -D).all()))