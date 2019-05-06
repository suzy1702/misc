#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:18:34 2019

@author: ty
"""
import numpy as np
pi = np.pi

###########################################
def comp_vol(vec):
    """
    """
    vol = vec[0,:].dot(np.cross(vec[1,:],vec[2,:]))
    
    return vol

###########################################
def make_kvec(rvec,rvol):
    """
    """
    kvec = np.zeros((3,3))
    kvec[0,:] = 2*pi*np.cross(rvec[1,:],rvec[2,:])/rvol
    kvec[1,:] = 2*pi*np.cross(rvec[2,:],rvec[0,:])/rvol
    kvec[2,:] = 2*pi*np.cross(rvec[0,:],rvec[1,:])/rvol
    kvol = comp_vol(kvec)
    
    return [kvec, kvol]

###########################################
def make_lattice(vec,n1,n2,n3):
    """
    """
    if n1%2 == 0:
        print('\tUSAGE ERROR: n1 should be an odd integer\n'
              '\t\tIncrementing n1 by 1\n')
        n1 = n1+1
    if n2%2 == 0:
        print('\tUSAGE ERROR: n2 should be an odd integer\n'
              '\t\tIncrementing n2 by 1\n')
        n2 = n2+1
    if n3%2 == 0:
        print('\tUSAGE ERROR: n3 should be an odd integer\n'
              '\t\tIncrementing n3 by 1\n')
        n3 = n3+1
    nc = n1*n2*n3
    pos = np.zeros((nc,3))
    count = 0
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                pos[count,:] = vec[0,:]*i+vec[1,:]*j+vec[2,:]*k
                count = count+1
                
    return pos