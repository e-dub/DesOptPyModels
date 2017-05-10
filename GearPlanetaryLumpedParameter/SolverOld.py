# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 15:16:44 2016

@author: wehrle
"""
import numpy as np
import scipy.linalg as linalg

def ModalUndamped(M=[], D=[], K=[], Hz=False, All=True, Sum=True):
    L, Phi = linalg.eig(K, M)
    if Sum:
        L = np.abs(L.real+L.imag)
        Phi = (Phi.real+Phi.imag)
    iSort = L.argsort()
    if All:
        omega = np.sqrt(L)[iSort]
        Phi = (Phi.real+Phi.imag)[iSort,:]
        #Phi = (Phi.real+Phi.imag)[:,iSort]
    else:
        omega = np.sqrt(L)[iSort][::2]
        Phi = (Phi.real+Phi.imag)[iSort,:][:,::2]
    if Hz:
        omega *= 1/(2*np.pi)
    return omega, Phi
#
#def ModalUndampedI(M, K, Hz=False, All=True):
#    L, Phi = linalg.eig(K, M)
#    #Lreal = np.abs(L.real)
#    iSort = L.argsort()
#    if All:
#        omega0 = np.sqrt(L)[iSort]
#        Phi = (Phi)[:,iSort]
#    else:
#        omega0 = np.sqrt(L)[iSort][::2]
#        Phi = (Phi)[:,iSort][:,::2]
#    if Hz:
#        f0 = omega0/(2*np.pi)
#    else:
#        f0 = omega0
#    return f0, Phi

def ModalDamped(M=[], D=[], K=[], Hz=False, All=False, Sum=True):
    A = np.concatenate((np.concatenate((M, np.zeros(np.shape(M))), 1),
                        np.concatenate((np.zeros(np.shape(M)), -K), 1)), 0)
    B = np.concatenate((np.concatenate((np.zeros(np.shape(M)), M), 1),
                        np.concatenate((M, D), 1)), 0)
    omega, Phi = linalg.eig(A, B)
    if Sum:
        omega = np.abs(omega.real+omega.imag)
        Phi = (Phi.real+Phi.imag)
    iSort = omega.argsort()
    if All:
        omega = omega[iSort]
        Phi = Phi[iSort,:]
    else:
        omega = omega[iSort][::2]
        Phi = Phi[iSort,:][:,iSort][::2,::2]
    if Hz:
        omega *= 1/(2*np.pi)
    else:
        pass
    return omega, Phi

#def ModalDampedI(M, D, K, Hz=False, All=False):
#    A = np.concatenate((np.concatenate((M, np.zeros(np.shape(M))), 1),
#                        np.concatenate((np.zeros(np.shape(M)), -K), 1)), 0)
#    B = np.concatenate((np.concatenate((np.zeros(np.shape(M)), M), 1),
#                        np.concatenate((M, D), 1)), 0)
#    L, Phi = linalg.eig(A, B)
#    iSort = L.argsort()
#    if All:
#        omega0 = L[iSort]
#        Phi = (Phi)[0:18,iSort]
#    else:
#        omega0 = L[iSort][::2]
#        Phi = (Phi)[:,iSort][0:18,::2]
#    if Hz:
#        f0 = omega0/(2*np.pi)
#    else:
#        f0 = omega0
#    return f0, Phi

def ModalSensUndamped(omega, Phi, M, K, dMdx, dKdx, Hz=False, All=True):
    dlmbdadx = np.zeros(np.shape(omega))
    domegadx = np.zeros(np.shape(omega))
    if np.shape(omega) is ():
        omega = omega.reshape(1,1)
        dlmbdadx = dlmbdadx.reshape(1,1)
        domegadx = dlmbdadx.reshape(1,1)
        Phi = Phi.reshape(1,np.size(Phi))
    lmbda = omega**2
    for ii in range(np.size(omega)):
        dlmbdadx[ii] = Phi[ii,:].dot(dKdx-lmbda[ii]*dMdx).dot(Phi[ii,:])/Phi[ii,:].dot(M).dot(Phi[ii,:])
        domegadx[ii] = (Phi[ii,:].dot(dKdx-lmbda[ii]*dMdx).dot(Phi[ii,:])/Phi[ii,:].dot(M).dot(Phi[ii,:]))**0.5
        #domegadx[ii] = (Phi[ii,:].dot(dKdx+lmbda[ii]*dMdx).dot(Phi[ii,:])/Phi[ii,:].dot(M).dot(Phi[ii,:]))**0.5
    return domegadx

def ModalSensDamped(omega, Phi, M, D, K, dMdx, dDdx, dKdx, Hz=False, All=True):
    domegadx = np.zeros(np.shape(omega))
    if np.shape(omega) is ():
        omega = omega.reshape(1,1)
        domegadx = domegadx.reshape(1,1)
        Phi = Phi.reshape(1,np.size(Phi))
    for ii in range(np.size(omega)):
        domegadx[ii] = (-omega[ii]**2*Phi[ii,:].dot(dMdx).dot(Phi[ii,:])-
                        omega[ii]*Phi[ii,:]. dot(dDdx).dot(Phi[ii,:])-
                        Phi[ii,:].dot(dKdx).dot(Phi[ii,:]))/(2*omega[ii]*Phi[ii,:].dot(M).dot(Phi[ii,:])+
                                                             Phi[ii,:].dot(D).dot(Phi[ii,:]))
    return domegadx


def ModalSensDampedNew(omega, Phi, M, D, K, dMdx, dDdx, dKdx, Hz=False, All=True):
    domegadx = np.zeros(np.shape(omega))
    if np.shape(omega) is ():
        omega = omega.reshape(1,1)
        domegadx = domegadx.reshape(1,1)
        Phi = Phi.reshape(1,np.size(Phi))
    for ii in range(np.size(omega)):
        domegadx[ii] = Phi[ii,:].dot(omega[ii]**2*dMdx+omega[ii]*dDdx+dKdx).dot(Phi[ii,:])
    return domegadx