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
    
    # Transformation
    QP = {}
    QSP = {}
    QCP = {}
    QRP = {}
    QS = np.zeros([3, 9+int(nP)*3])         # Transformation matrix sun gear
    QS[:, 0:3] = np.eye(3)
    QC = np.zeros([3, 9+int(nP)*3])         # Transformation matrix carrier
    QC[:, 3:6] = np.eye(3)
    QR = np.zeros([3, 9+int(nP)*3])         # Transformation matrix ring gear
    QR[:, 6:9] = np.eye(3)
    for ii in range(int(nP)):
        QP[ii] = np.zeros([3, 9+int(nP)*3])
        
        QP[ii][0, 9+3*ii+0] = 1
        QP[ii][1, 9+3*ii+1] = 1
        QP[ii][2, 9+3*ii+2] = 1        
        QSP[ii] = np.concatenate((QS, QP[ii]),axis=0)
        QCP[ii] = np.concatenate((QC, QP[ii]),axis=0)
        QRP[ii] = np.concatenate((QR, QP[ii]),axis=0)
        
    # Mass
    MS = [mS, mS, JS]*np.eye(3)
    MC = [mC, mC, JC]*np.eye(3)
    MR = [mR, mR, JR]*np.eye(3)
    MP = {}
    for ii in range(int(nP)):
        MP[ii] = [mP, mP, JP]*np.eye(3)
        
    # Stiffness
    KS = [kSx, kSy, kStheta]*np.eye(3)
    KC = [kCx, kCy, kCtheta]*np.eye(3)
    KR = [kRx, kRy, kRtheta]*np.eye(3)
    KP = {}
    KSP = {}
    KCP = {}
    KRP = {}
    KGyroP = {}
    KDragP = {}
    for ii in range(int(nP)):
        phi = 2*np.pi/nP*ii+0.25  
        KP[ii] = [kCP, kCP, 0]*np.eye(3)
        # Sun-planet stiffnesses  
        KSP[ii] = kSP*np.array([[ sin(phi-alpha)*sin(phi-alpha), -sin(phi-alpha)*cos(phi-alpha), -sin(phi-alpha)*rS,  sin(phi-alpha)*sin(alpha),  sin(phi-alpha)*cos(alpha), -sin(phi-alpha)*rP],
                                [-cos(phi-alpha)*sin(phi-alpha),  cos(phi-alpha)*cos(phi-alpha),  cos(phi-alpha)*rS, -cos(phi-alpha)*sin(alpha), -cos(phi-alpha)*cos(alpha),  cos(phi-alpha)*rP],
                                [    -rS*sin(phi-alpha),                      rS*cos(phi-alpha),              rS*rS,             -rS*sin(alpha),             -rS*cos(alpha),              rS*rP],
                                [ sin(alpha)*sin(phi-alpha),         -sin(alpha)*cos(phi-alpha),     -sin(alpha)*rS,      sin(alpha)*sin(alpha),      sin(alpha)*cos(alpha),     -sin(alpha)*rP],
                                [ cos(alpha)*sin(phi-alpha),         -cos(alpha)*cos(phi-alpha),     -cos(alpha)*rS,      cos(alpha)*sin(alpha),      cos(alpha)*cos(alpha),     -cos(alpha)*rP],
                                [    -rP*sin(phi-alpha),                      rP*cos(phi-alpha),              rP*rS,             -rP*sin(alpha),             -rP*cos(alpha),              rP*rP]])   

        # Carrier-planet stiffnesses
        KCP[ii] = kCP*np.array([[         1,         0, -rC*sin(phi), -cos(phi),  sin(phi), 0.],
                                [         0,         1,  rC*cos(phi), -sin(phi), -cos(phi), 0.],
                                [-rC*sin(phi), rC*cos(phi),      rC*rC,      0.,     -rC, 0.],
                                [   -cos(phi),   -sin(phi),         0.,      1.,      0., 0.],
                                [    sin(phi),   -cos(phi),        -rC,      0.,      1., 0.],
                                [        0.,        0.,         0.,      0.,      0., 0.]])

        KRP[ii] = kRP*np.array([[ sin(phi+alpha)*sin(phi+alpha), -sin(phi+alpha)*cos(phi+alpha), -sin(phi+alpha)*rR, -sin(phi+alpha)*sin(alpha),  sin(phi+alpha)*cos(alpha),  sin(phi+alpha)*rP],
                                [-cos(phi+alpha)*sin(phi+alpha),  cos(phi+alpha)*cos(phi+alpha),  cos(phi+alpha)*rR,  cos(phi+alpha)*sin(alpha), -cos(phi+alpha)*cos(alpha), -cos(phi+alpha)*rP],
                                [    -rR*sin(phi+alpha),      rR*cos(phi+alpha),      rR*rR,      rR*sin(alpha),     -rR*cos(alpha),     -rR*rP],
                                [-sin(alpha)*sin(phi+alpha),  sin(alpha)*cos(phi+alpha),  sin(alpha)*rR,  sin(alpha)*sin(alpha), -sin(alpha)*cos(alpha), -sin(alpha)*rP],
                                [ cos(alpha)*sin(phi+alpha), -cos(alpha)*cos(phi+alpha), -cos(alpha)*rR, -cos(alpha)*sin(alpha),  cos(alpha)*cos(alpha),  cos(alpha)*rP],
                                [     rP*sin(phi+alpha),     -rP*cos(phi+alpha),     -rP*rR,     -rP*sin(alpha),      rP*cos(alpha),      rP*rP]])
        # Gyroscopic stiffness matrices
        KGyroP[ii] = -mP*thetaCdot**2*np.array([[1., 0., 0.],
                                                [0., 1., 0.],
                                                [0., 0., 0.]])
        # Dynamic drag stiffness matrices
        KDragP[ii] = mP*thetaCdotdot*np.array([[0., 1., 0.],
                                               [1., 0., 0.],
                                               [0., 0., 0.]])
    # Dymamic matrix
    DGyroP = {}
    for ii in range(int(nP)):
        DGyroP[ii] = 2*mP*thetaCdot*np.array([[0., -1., 0.],
                                              [1.,  0., 0.],
                                              [0.,  0., 0.]])


    M = QS.T.dot(MS.dot(QS)) + \
        QC.T.dot(MC.dot(QC)) + \
        QR.T.dot(MR.dot(QR))
    for ii in range(int(nP)):
        M += QP[ii].T.dot(MP[ii].dot(QP[ii]))
    # Gyroscopic damping matrices
    D = QP[ii].T.dot(DGyroP[ii].dot(QP[ii]))
    All={}
    IncludeGyro = 1
    IncludeDrag = 0
    K = QS.T.dot(KS.dot(QS)) + \
        QC.T.dot(KC.dot(QC)) + \
        QR.T.dot(KR.dot(QR))
    All["KS"] = QS.T.dot(KS.dot(QS))
    All["KC"] = QC.T.dot(KC.dot(QC))
    All["KR"] = QR.T.dot(KR.dot(QR))
    print((K.transpose() == K).all())
    for ii in range(int(nP)):
        K += QSP[ii].T.dot(KSP[ii].dot(QSP[ii]))
        All["KSP"+str(ii)] = QSP[ii].T.dot(KSP[ii].dot(QSP[ii]))
    print((K.transpose() == K).all())
    for ii in range(int(nP)):
        K += QCP[ii].T.dot(KCP[ii].dot(QCP[ii]))
        All["KCP"+str(ii)] = QCP[ii].T.dot(KCP[ii].dot(QCP[ii]))
    print((K.transpose() == K).all())
    for ii in range(int(nP)):
        K += QRP[ii].T.dot(KRP[ii].dot(QRP[ii]))
        All["KRP"+str(ii)] = QRP[ii].T.dot(KRP[ii].dot(QRP[ii]))
    print((K.transpose() == K).all())
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
    print K