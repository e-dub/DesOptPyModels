from __future__ import division
import numpy as np
import sympy as sp


symbolic = False
GearSet = "TQ070i3"
GearSet = "MP080i6"
GearSet = "Me"

if GearSet=="TQ070i3":
    rS=32.4/2
    rC=24.3
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
    thetaCdot = 1000./60.*2*np.pi
    thetaCdotdot = 0
elif GearSet=="MP080i6":
    rS = 32.4/2.
    rC = 24.3
    rR = 64.8/2
    rP = 15.6/2
    mS = 0.262e-3
    JS = 21e-3
    mC = 0.44e-3
    JC = 190e-3
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
elif GearSet=="Me":
    rS = 25.
    rC = 37.5
    rR = 100.
    rP = 12.5
    mS = 0.5e-3
    JS = 25.0e-3
    mC = 0.6e-3
    JC = 200.0e-3
    mR = 2.0e-3
    JR = 3000.0e-3
    mP = 0.05e-3
    JP = 10.0e-3
    kSx = 200e3
    kSy = 200e3
    kStheta = 5000e3
    kCx = 250e3
    kCy = 250e3
    kCtheta = 30000e3
    kRx = 750e8
    kRy = 750e8
    kRtheta = 1e10
    kSP = 200e3
    kCP = 1000e3
    kRP = 200e3
    nP=4.
    alpha=20.
    thetaCdot=20000./60.*2*np.pi
    thetaCdotdot=0

def CalcMatrices(mS=0, JS=0, mC=0, JC=0, mR=0, JR=0, mP=0, JP=0, 
                 kSx=0, kSy=0, kStheta=0, kCx=0, kCy=0, kCtheta=0, 
                 kRx=0, kRy=0, kRtheta=0, kSP=0, kCP=0, kRP=0, 
                 nP=3, alpha=20, rS=0, rC=0, rR=0, rP=0, 
                 thetaCdot=0, thetaCdotdot=0):
    if symbolic:
        mS, JS, mC, JC, mR, JR, mP, JP = sp.symbols('mS JS mC JC mR JR mP JP')
        kSx, kSy, kStheta = sp.symbols('kSx, kSy, kStheta')
        kCx, kCy, kCtheta = sp.symbols('kCx kCy, kCtheta')
        kRx, kRy, kRtheta = sp.symbols('kRx, kRy, kRtheta')
        kSP, kCP, kRP = sp.symbols(' kSP, kCP, kRP')
        alpha, rS, rC, rR, rP = sp.symbols('alpha rS rC rR rP')
        thetaCdot, thetaCdotdot= sp.symbols('thetaCdot thetaCdotdot')
        from sympy import cos, sin 
    else:
        from numpy import cos, sin, deg2rad
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
    # Sun-planet stiffnesses
    KSP = {}
    QSP = {}
    for ii in range(int(nP)):
        phi = 360./nP*ii
        if symbolic:
            A = phi + sp.symbols('90-alpha')
            B = sp.symbols('90-alpha')
        else:     
            A = deg2rad(phi+90-alpha)
            B = deg2rad(90-alpha)
        KSP[ii] = kSP*np.array([[-cos(A)*cos(A), -cos(A)*sin(A), -rS*cos(A),  cos(A)*sin(B),  cos(A)*cos(B), -rP*cos(A)],
                                [-sin(A)*cos(A), -sin(A)*sin(A), -rS*sin(A),  sin(A)*sin(B),  sin(A)*cos(B), -rP*sin(A)],
                                [    -rS*cos(A),     -rS*sin(A),     -rS**2,      rS*sin(B),      rS*cos(B),     -rP*rS],
                                [ sin(B)*cos(A),  sin(B)*sin(A),  rS*sin(B), -sin(B)*sin(B), -sin(B)*cos(B),  rP*sin(B)],
                                [ cos(B)*cos(A),  cos(B)*sin(A),  rS*cos(B), -cos(B)*sin(B), -cos(B)*cos(B),  rP*cos(B)],
                                [    -rP*cos(A),     -rP*sin(A),     -rS*rP,      rP*sin(B),      rP*cos(B),     -rP**2]])    
        QSP[ii] = np.zeros([6, 9+int(nP)*3])
        QSP[ii][0, 0] = 1.
        QSP[ii][1, 1] = 1.
        QSP[ii][2, 2] = 1.
        QSP[ii][3, 9+3*ii+0] = 1.
        QSP[ii][4, 9+3*ii+1] = 1.
        QSP[ii][5, 9+3*ii+2] = 1.
    # Carrier-planet stiffnesses
    KCP = {}
    QCP = {}
    for ii in range(int(nP)):
        phi = 360./nP*ii
        if symbolic:
            A = phi
        else:
            A = deg2rad(phi)
        KCP[ii] = kCP*np.array([[-cos(A)*cos(A)-sin(A)*sin(A),  cos(A)*sin(A)+sin(A)*cos(A),  rC*sin(A),  cos(A), -sin(A), 0.],
                                [ sin(A)*cos(A)+cos(A)*sin(A), -sin(A)*sin(A)-cos(A)*cos(A), -rC*cos(A), -sin(A),  cos(A), 0.],
                                [                   rC*sin(A),                   -rC*cos(A),     -rC**2,      0.,      rC, 0.],
                                [                      cos(A),                      -sin(A),         0.,     -1.,      0., 0.],
                                [                     -sin(A),                       cos(A),         rC,      0.,     -1., 0.],
                                [                          0.,                           0.,         0.,      0.,      0., 0.]])
        QCP[ii] = np.zeros([6, 9+int(nP)*3])
        QCP[ii][0, 3] = 1.
        QCP[ii][1, 4] = 1.
        QCP[ii][2, 5] = 1.
        QCP[ii][3, 9+3*ii+0] = 1.
        QCP[ii][4, 9+3*ii+1] = 1.
        QCP[ii][5, 9+3*ii+2] = 1.
    # Ring-planet stiffnesses
    KRP = {}
    QRP = {}
    for ii in range(int(nP)):
        phi = 360/nP*ii
        if symbolic:
            A = phi + sp.symbols('270-alpha')
            B = sp.symbols('90-alpha')
        else:   
            A = deg2rad(phi+270-alpha)
            B = deg2rad(90-alpha)
        KRP[ii] = kRP*np.array([[ cos(A)*cos(A),  cos(A)*sin(A),  rR*cos(A), -cos(A)*cos(B), -cos(A)*sin(B), -rP*cos(A)],
                                [ sin(A)*cos(A),  sin(A)*sin(A),  rR*sin(A), -sin(A)*cos(B), -sin(A)*sin(B), -rP*sin(A)],
                                [     rR*cos(A),      rR*sin(A),      rR**2,     -rR*cos(B),     -rR*sin(B),     -rR*rP],
                                [-cos(B)*cos(A), -cos(B)*sin(A), -rR*cos(B),  cos(B)*cos(B),  cos(B)*sin(B),  rP*cos(B)],
                                [-sin(B)*cos(A), -sin(B)*sin(A), -rR*sin(B),  sin(B)*cos(B),  sin(B)*sin(B),  rP*sin(B)],
                                [    -rP*cos(A),     -rP*sin(A),    -rR*rP,       rP*cos(B),      rP*sin(B),     -rP**2]])
        QRP[ii] = np.zeros([6, 9+int(nP)*3])
        QRP[ii][0, 6] = 1.
        QRP[ii][1, 7] = 1.
        QRP[ii][2, 8] = 1.
        QRP[ii][3, 9+3*ii+0] = 1.
        QRP[ii][4, 9+3*ii+1] = 1.
        QRP[ii][5, 9+3*ii+2] = 1.
    # Gyroscopic stiffness matrices
    KGyroP = {}
    for ii in range(int(nP)):
        KGyroP[ii] = -mP*thetaCdot**2*np.array([[1., 0., 0.],
                                                [0., 1., 0.],
                                                [0., 0., 0.]])
    # Dynamic drag stiffness matrices
    KDragP = {}
    for ii in range(int(nP)):
        KDragP[ii] = mP*thetaCdotdot*np.array([[0., 1., 0.],
                                               [1., 0., 0.],
                                               [0., 0., 0.]])
    IncludeGyro = 1
    IncludeDrag = 1
    K = QS.T.dot(KS.dot(QS)) + \
        QC.T.dot(KC.dot(QC)) + \
        QR.T.dot(KR.dot(QR))
    for ii in range(int(nP)):
        K += QSP[ii].T.dot(KSP[ii].dot(QSP[ii]))
        K += QCP[ii].T.dot(KCP[ii].dot(QCP[ii]))
        K += QRP[ii].T.dot(KRP[ii].dot(QRP[ii]))
        if IncludeDrag:
            K += QP[ii].T.dot(KDragP[ii].dot(QP[ii]))
        if IncludeGyro:
            K += QP[ii].T.dot(KGyroP[ii].dot(QP[ii]))
    return(M, D, K)


def SysEval():
    M, D, K = CalcMatrices(mS=mS, JS=JS, mC=mC, JC=JC, mR=mR, JR=JR, mP=mP, JP=JP, 
                           kSx=kSx, kSy=kSy, kStheta=kStheta, kCx=kCx, kCy=kCy, 
                           kCtheta=kCtheta, kRx=kRx, kRy=kRy, kRtheta=kRtheta, 
                           kSP=kSP, kCP=kCP, kRP=kRP, nP=nP, alpha=alpha, rS=rS, 
                           rC=rC, rR=rR, rP=rP, thetaCdot=thetaCdot, 
                           thetaCdotdot=0)
    import Solver as OldSolver
    from LumPy import Solver
    omega, Phi = Solver.ModalDamped(M, D, K, Hz=True, All=False, Sum=True)
    return(omega, Phi)

omega, Phi = SysEval()
print omega
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