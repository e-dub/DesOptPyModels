from __future__ import division
import numpy as np
from numpy import cos, sin


def stiffness(N, Rot, ks, kp, kc, kr, kTs, kTr, kTpt, ksx, ksy, krx, kry, mp,
              rs, rp, rc, rr, alpha, phi):
    K = np.zeros([N, N])
    # Sun gear
    # theta_s
    K[0, 0] = kTs+ks*rs**2                                                      # theta_s
    K[0, 1] = ks*rs*(sin(alpha-phi[0])+sin(alpha-phi[1])+sin(alpha-phi[2]))     # x_s
    K[0, 2] = ks*rs*(cos(alpha-phi[0])+cos(alpha-phi[1])+cos(alpha-phi[2]))     # y_s
    K[0, 3] = ks*rs*rp                                                          # theta_p1
    K[0, 4] = -ks*rs*sin(alpha)                                                 # eta_p1
    K[0, 5] = -ks*rs*cos(alpha)                                                 # xi_p1
    K[0, 6] = ks*rs*rp                                                          # theta_p2
    K[0, 7] = -ks*rs*sin(alpha)                                                 # eta_p2
    K[0, 8] = -ks*rs*cos(alpha)                                                 # xi_p2
    K[0, 9] = ks*rs*rp                                                          # theta_p3
    K[0, 10] = -ks*rs*sin(alpha)                                                # eta_p3
    K[0, 11] = -ks*rs*cos(alpha)                                                # xi_p3
    # x_s
    K[1, 0] = ks*rs*(sin(alpha-phi[0])+sin(alpha-phi[1])+sin(alpha-phi[2]))      # theta_s
    K[1, 1] = ksx+ks*sin(alpha-phi[0])*sin(alpha-phi[0])+ks*sin(alpha-phi[1])*sin(alpha-phi[1])+ks*sin(alpha-phi[2])*sin(alpha-phi[2])  # x_s
    K[1, 2] = ks*sin(alpha-phi[0])*cos(alpha-phi[0])+ks*sin(alpha-phi[1])*cos(alpha-phi[1])+ks*sin(alpha-phi[2])*cos(alpha-phi[2])  # y_s
    K[1, 3] = ks*rp*sin(alpha-phi[0])                                           # theta_p1
    K[1, 4] = -ks*sin(alpha)*sin(alpha-phi[0])                                  # eta_p1
    K[1, 5] = -ks*cos(alpha)*sin(alpha-phi[0])                                  # xi_p1
    K[1, 6] = ks*rp*sin(alpha-phi[1])                                           # theta_p2
    K[1, 7] = -ks*sin(alpha)*sin(alpha-phi[1])                                  # eta_p2
    K[1, 8] = -ks*cos(alpha)*sin(alpha-phi[1])                                  # xi_p2
    K[1, 9] = ks*rp*sin(alpha-phi[2])                                           # theta_p3
    K[1, 10] = -ks*sin(alpha)*sin(alpha-phi[2])                                 # eta_p3
    K[1, 11] = -ks*cos(alpha)*sin(alpha-phi[2])                                 # xi_p3
    # y_s
    K[2, 0] = ks*rs*cos(alpha-phi[0])+ks*rs*cos(alpha-phi[1])+ks*rs*cos(alpha-phi[2])
    K[2, 1] = ks*sin(alpha-phi[0])*cos(alpha-phi[0])+ks*sin(alpha-phi[1])*cos(alpha-phi[1])+ks*sin(alpha-phi[2])*cos(alpha-phi[2])
    K[2, 2] = ksy+ks*cos(alpha-phi[0])*cos(alpha-phi[0])+ks*cos(alpha-phi[1])*cos(alpha-phi[1])+ks*cos(alpha-phi[2])*cos(alpha-phi[2])
    K[2, 3] = ks*rp*cos(alpha-phi[0])
    K[2, 4] = -ks*sin(alpha)*cos(alpha-phi[0])
    K[2, 5] = -ks*cos(alpha)*cos(alpha-phi[0])
    K[2, 6] = ks*rp*cos(alpha-phi[1])
    K[2, 7] = -ks*sin(alpha)*cos(alpha-phi[1])
    K[2, 8] = -ks*cos(alpha)*cos(alpha-phi[1])
    K[2, 9] = ks*rp*cos(alpha-phi[2])
    K[2, 10] = -ks*sin(alpha)*cos(alpha-phi[2])
    K[2, 11] = -ks*cos(alpha)*cos(alpha-phi[2])

    # Planet 1
    # theta_p1
    K[3, 0] = ks*rp*rs
    K[3, 1] = ks*rp*sin(alpha-phi[0])
    K[3, 2] = ks*rp*cos(alpha-phi[0])
    K[3, 3] = ks*rp**2+kr*rp**2
    K[3, 4] = -ks*rp*sin(alpha)-kr*rp*sin(alpha)
    K[3, 5] = -ks*rp*cos(alpha)+kr*rp*cos(alpha)
    K[3, 15] = -kr*rp*rr
    K[3, 16] = kr*rp*cos(alpha+phi[0])
    K[3, 17] = -kr*rp*sin(alpha+phi[0])
    # eta_p1
    K[4, 0] = ks*rs*sin(alpha)
    K[4, 1] = ks*sin(alpha)*sin(alpha-phi[0])
    K[4, 2] = ks*sin(alpha)*cos(alpha-phi[0])
    K[4, 3] = ks*rp*sin(alpha)+kr*rp*sin(alpha)
    K[4, 4] = -ks*sin(alpha)*sin(alpha)-kr*sin(alpha)*sin(alpha)-kp+Rot**2*mp
    K[4, 5] = -ks*sin(alpha)*cos(alpha)+kr*sin(alpha)*cos(alpha)
    K[4, 13] = -kp*cos(phi[0])
    K[4, 14] = -kp*sin(phi[0])
    K[4, 15] = -kr*rr*sin(alpha)
    K[4, 16] = -kr*sin(alpha)*cos(alpha+phi[0])
    K[4, 17] = kr*sin(alpha)*sin(alpha+phi[0])
    # xi_p1
    K[5, 0] = ks*rs*cos(alpha)
    K[5, 1] = ks*cos(alpha)*sin(alpha-phi[0])
    K[5, 2] = ks*cos(alpha)*cos(alpha-phi[0])
    K[5, 3] = ks*rp*cos(alpha)-kr*rp*cos(alpha)
    K[5, 4] = -ks*cos(alpha)*sin(alpha)+kr*sin(alpha)*cos(alpha)
    K[5, 5] = -kp-ks*cos(alpha)*cos(alpha)-kr*cos(alpha)*cos(alpha)-Rot**2*mp
    K[5, 13] = -kp*sin(phi[0])
    K[5, 14] = kp*cos(phi[0])
    K[5, 15] = -kr*rr*cos(alpha)
    K[5, 16] = kr*cos(alpha)*cos(alpha+phi[0])
    K[5, 17] = -kr*cos(alpha)*sin(alpha+phi[0])

    # Planet 2
    # theta_p2
    K[6, 0] = ks*rp*rs
    K[6, 1] = ks*rp*sin(alpha-phi[1])
    K[6, 2] = ks*rp*cos(alpha-phi[1])
    K[6, 6] = ks*rp**2+kr*rp**2
    K[6, 7] = -ks*rp*sin(alpha)-kr*rp*sin(alpha)
    K[6, 8] = -ks*rp*cos(alpha)+kr*rp*cos(alpha)
    K[6, 15] = -kr*rp*rr
    K[6, 16] = kr*rp*cos(alpha+phi[1])
    K[6, 17] = -kr*rp*sin(alpha+phi[1])
    # eta_p2
    K[7, 0] = ks*rs*sin(alpha)
    K[7, 1] = ks*sin(alpha)*sin(alpha-phi[1])
    K[7, 2] = ks*sin(alpha)*cos(alpha-phi[1])
    K[7, 6] = ks*rp*sin(alpha)+kr*rp*sin(alpha)
    K[7, 7] = -ks*sin(alpha)*sin(alpha)-kr*sin(alpha)*sin(alpha)-kp+Rot**2*mp
    K[7, 8] = -ks*sin(alpha)*cos(alpha)+kr*sin(alpha)*cos(alpha)
    K[7, 13] = -kp*cos(phi[1])
    K[7, 14] = -kp*sin(phi[1])
    K[7, 15] = -kr*rr*sin(alpha)
    K[7, 16] = -kr*sin(alpha)*cos(alpha+phi[1])
    K[7, 17] = kr*sin(alpha)*sin(alpha+phi[1])
    # xi_p2
    K[8, 0] = ks*rs*cos(alpha)
    K[8, 1] = ks*cos(alpha)*sin(alpha-phi[1])
    K[8, 2] = ks*cos(alpha)*cos(alpha-phi[1])
    K[8, 6] = ks*rp*cos(alpha)-kr*rp*cos(alpha)
    K[8, 7] = -ks*cos(alpha)*sin(alpha)+kr*sin(alpha)*cos(alpha)
    K[8, 8] = -kp-ks*cos(alpha)*cos(alpha)-kr*cos(alpha)*cos(alpha)-Rot**2*mp
    K[8, 13] = -kp*sin(phi[1])
    K[8, 14] = kp*cos(phi[1])
    K[8, 15] = -kr*rr*cos(alpha)
    K[8, 16] = kr*cos(alpha)*cos(alpha+phi[1])
    K[8, 17] = -kr*cos(alpha)*sin(alpha+phi[1])

    # Planet 3
    # theta_p3
    K[9, 0] = ks*rp*rs
    K[9, 1] = ks*rp*sin(alpha-phi[2])
    K[9, 2] = ks*rp*cos(alpha-phi[2])
    K[9, 9] = ks*rp**2+kr*rp**2
    K[9, 10] = -ks*rp*sin(alpha)-kr*rp*sin(alpha)
    K[9, 11] = -ks*rp*cos(alpha)+kr*rp*cos(alpha)
    K[9, 15] = -kr*rp*rr
    K[9, 16] = kr*rp*cos(alpha+phi[2])
    K[9, 17] = -kr*rp*sin(alpha+phi[2])
    # eta_p3
    K[10, 0] = ks*rs*sin(alpha)
    K[10, 1] = ks*sin(alpha)*sin(alpha-phi[2])
    K[10, 2] = ks*sin(alpha)*cos(alpha-phi[2])
    K[10, 9] = ks*rp*sin(alpha)+kr*rp*sin(alpha)
    K[10, 10] = -ks*sin(alpha)*sin(alpha)-kr*sin(alpha)*sin(alpha)-kp+Rot**2*mp
    K[10, 11] = -ks*sin(alpha)*cos(alpha)+kr*sin(alpha)*cos(alpha)
    K[10, 13] = -kp*cos(phi[2])
    K[10, 14] = -kp*sin(phi[2])
    K[10, 15] = -kr*rr*sin(alpha)
    K[10, 16] = -kr*sin(alpha)*cos(alpha+phi[2])
    K[10, 17] = kr*sin(alpha)*sin(alpha+phi[2])
    # xip3
    K[11, 0] = ks*rs*cos(alpha)
    K[11, 1] = ks*cos(alpha)*sin(alpha-phi[2])
    K[11, 2] = ks*cos(alpha)*cos(alpha-phi[2])
    K[11, 9] = ks*rp*cos(alpha)-kr*rp*cos(alpha)
    K[11, 10] = -ks*cos(alpha)*sin(alpha)+kr*sin(alpha)*cos(alpha)
    K[11, 11] = -kp-ks*cos(alpha)*cos(alpha)-kr*cos(alpha)*cos(alpha)-Rot**2*mp
    K[11, 13] = -kp*sin(phi[2])
    K[11, 14] = kp*cos(phi[2])
    K[11, 15] = -kr*rr*cos(alpha)
    K[11, 16] = kr*cos(alpha)*cos(alpha+phi[2])
    K[11, 17] = -kr*cos(alpha)*sin(alpha+phi[2])

    # Carrier
    # theta_c (was pt!)
    K[12, 5] = kp*rc            # eta_p1
    K[12, 8] = kp*rc            # eta_p2
    K[12, 11] = kp*rc           # eta_p3
    K[12, 12] = kTpt           #  theta_c
    K[12, 13] = kp*rc*sin(phi[0])+kp*rc*sin(phi[1])+kp*rc*sin(phi[2]) #x_c
    K[12, 14] = -kp*rc*cos(phi[0])-kp*rc*cos(phi[1])-kp*rc*cos(phi[2])  #y_c
    K[12, 12] = -kp*rc**2                      #WRONGWRONGWRONGWRONGWRONGWRONG FUCK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
    #missing k's?
    #K[12, 13] = -kTpt
    #K[12, 14] = -kTpt
    #x_c
    K[13, 4] = kp*cos(phi[0])
    K[13, 5] = -kp*sin(phi[0])
    K[13, 7] = kp*cos(phi[1])
    K[13, 8] = -kp*sin(phi[1])
    K[13, 10] = kp*cos(phi[2])
    K[13, 11] = -kp*sin(phi[2])
    K[13, 12] = kp*rc*sin(phi[0])+kp*rc*sin(phi[1])+kp*rc*sin(phi[2])
    K[13, 13] = -kc-kp*cos(phi[0])*cos(phi[0])-kp*cos(phi[1])*cos(phi[1])-kp*cos(phi[2])*cos(phi[2])-kp*sin(phi[0])*sin(phi[0])-kp*sin(phi[1])*sin(phi[1])-kp*sin(phi[2])*sin(phi[2])
    K[13, 14] = -kp*cos(phi[0])*sin(phi[0])-kp*cos(phi[1])*sin(phi[1])-kp*cos(phi[2])*sin(phi[2])+kp*sin(phi[0])*cos(phi[0])+kp*sin(phi[1])*cos(phi[1])+kp*sin(phi[2])*cos(phi[2])
    #y_c
    K[14, 4] = kp*sin(phi[0])
    K[14, 5] = kp*cos(phi[0])
    K[14, 7] = kp*sin(phi[1])
    K[14, 8] = kp*cos(phi[1])
    K[14, 10] = kp*sin(phi[2])
    K[14, 11] = kp*cos(phi[2])
    K[14, 12] = -kp*rc*cos(phi[0])-kp*rc*cos(phi[1])-kp*rc*cos(phi[2])
    K[14, 13] = -kp*sin(phi[0])*cos(phi[0])-kp*sin(phi[1])*cos(phi[1])-kp*sin(phi[2])*cos(phi[2])+kp*sin(phi[0])*cos(phi[0])+kp*sin(phi[1])*cos(phi[1])+kp*sin(phi[2])*cos(phi[2])
    K[14, 14] = -kc-kp*sin(phi[0])*sin(phi[0])-kp*sin(phi[1])*sin(phi[1])-kp*sin(phi[2])*sin(phi[2])-kp*cos(phi[0])*cos(phi[0])-kp*cos(phi[1])*cos(phi[1])-kp*cos(phi[2])*cos(phi[2])

    # RING
    # theta_r
    K[15, 3] = kr*rr*rp*cos(alpha+phi[0])
    K[15, 4] = -kr*rr*sin(alpha)*cos(alpha+phi[0])
    K[15, 5] = kr*rr*cos(alpha)*sin(alpha+phi[0])
    K[15, 6] = kr*rr*rp*cos(alpha+phi[1])
    K[15, 7] = -kr*rr*sin(alpha)*cos(alpha+phi[1])
    K[15, 8] = kr*rr*cos(alpha)*sin(alpha+phi[1])
    K[15, 9] = kr*rr*rp*cos(alpha+phi[2])
    K[15, 10] = -kr*rr*sin(alpha)*cos(alpha+phi[2])
    K[15, 11] = kr*rr*cos(alpha)*sin(alpha+phi[2])
    K[15, 15] = -kTr-kr*rr**2*cos(alpha+phi[0])-kr*rr**2*cos(alpha+phi[1])-kr*rr**2*cos(alpha+phi[2])
    K[15, 16] = -rr*kr*cos(alpha+phi[0])*cos(alpha+phi[0])-rr*kr*cos(alpha+phi[1])*cos(alpha+phi[1])-rr*kr*cos(alpha+phi[2])*cos(alpha+phi[2])
    K[15, 17] = +rr*kr*cos(alpha+phi[0])*sin(alpha+phi[0])+rr*kr*cos(alpha+phi[1])*sin(alpha+phi[1])+rr*kr*cos(alpha+phi[2])*sin(alpha+phi[2])
    # x_r
    K[16, 3] = -kr*rp*sin(alpha+phi[0])
    K[16, 4] = kr*sin(alpha)*sin(alpha+phi[0])
    K[16, 5] = -kr*cos(alpha)*sin(alpha+phi[0])
    K[16, 6] = -kr*rp*sin(alpha+phi[1])
    K[16, 7] = +kr*sin(alpha)*sin(alpha+phi[1])
    K[16, 8] = -kr*cos(alpha)*sin(alpha+phi[1])
    K[16, 9] = -kr*rp*sin(alpha+phi[2])
    K[16, 10] = +kr*sin(alpha)*sin(alpha+phi[2])
    K[16, 11] = -kr*cos(alpha)*sin(alpha+phi[2])
    K[16, 15] = +kr*rr*sin(alpha+phi[0])+kr*rr*sin(alpha+phi[1])+kr*rr*sin(alpha+phi[2])
    K[16, 16] = -krx+kr*cos(alpha+phi[0])*sin(alpha+phi[0])+kr*cos(alpha+phi[1])*sin(alpha+phi[1])+kr*cos(alpha+phi[2])*sin(alpha+phi[2])
    K[16, 17] = -kr*sin(alpha+phi[0])*sin(alpha+phi[0])-kr*sin(alpha+phi[1])*sin(alpha+phi[1])-kr*sin(alpha+phi[2])*sin(alpha+phi[2])
    # y_r
    K[17, 3] = kr*rp*cos(alpha+phi[0])
    K[17, 4] = -kr*sin(alpha)*cos(alpha+phi[0])
    K[17, 5] = kr*cos(alpha)*cos(alpha+phi[0])
    K[17, 6] = kr*rp*cos(alpha+phi[1])
    K[17, 7] = -kr*sin(alpha)*cos(alpha+phi[1])
    K[17, 8] = kr*cos(alpha)*cos(alpha+phi[1])
    K[17, 9] = kr*rp*cos(alpha+phi[2])
    K[17, 10] = kr*sin(alpha)*cos(alpha+phi[2])
    K[17, 11] = kr*cos(alpha)*cos(alpha+phi[2])
    K[17, 15] = -kr*rr*cos(alpha+phi[0])-kr*rr*cos(alpha+phi[1])-kr*rr*cos(alpha+phi[2])
    K[17, 16] = -kr*cos(alpha+phi[0])*cos(alpha+phi[0])-kr*cos(alpha+phi[1])*cos(alpha+phi[1])-kr*cos(alpha+phi[2])*cos(alpha+phi[2])
    K[17, 17] = -kry+kr*cos(alpha+phi[0])*sin(alpha+phi[0])+kr*cos(alpha+phi[1])*sin(alpha+phi[1])+kr*cos(alpha+phi[2])*sin(alpha+phi[2])
    return K
