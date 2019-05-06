#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A code that takes real space lattice vectors and creates the k space
lattice and the Wigner-Seitz cells in both spaces. 

rvec
"""
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import latticeFX as fx
import matplotlib.pyplot as plt
import copy as cp

plot = 'yes'

rvec = np.array([[1,0,0],
                 [0,1,0],
                 [0,0,1]]) #simple cubic

rvol = fx.comp_vol(rvec)
kvec, kvol = fx.make_kvec(rvec,rvol)

n1, n2, n3 = [5 ,5 ,5] #should all be odd integer numbers
rpos = fx.make_lattice(rvec,n1,n2,n3)
kpos = fx.make_lattice(kvec,n1,n2,n3)

##############################################################
# try --- scipy.spatial.Voronoi
##############################################################
#pos = cp.deepcopy(rpos)
#num = len(pos[:,0])
#xmid = np.median(pos[:,0])
#ymid = np.median(pos[:,1])
#zmid = np.median(pos[:,2])
#
#mid = np.intersect1d(np.argwhere(pos[:,0] == xmid),
#                     np.intersect1d(np.argwhere(pos[:,1] == ymid),
#                                    np.argwhere(pos[:,2] == zmid)))
#
#nl = np.zeros(num) #list of neighbor ids
#nd = np.zeros(num) #list of distance to neighbors
#rij = np.zeros((num,3)) #list of vector to neighbors
#
#for i in range(num): #loop over particles
#    rx = xmid-pos[i,0] # vector from partice i to j
#    ry = ymid-pos[i,1] # vector from partice i to j
#    rz = zmid-pos[i,2] # vector from partice i to j
#    dist = np.sqrt(rx**2+ry**2+rz**2)
#    
#    rij[i,0] = -rx
#    rij[i,1] = -ry
#    rij[i,2] = -rz
#    
#    nl[i] = i
#    nd[i] = dist
#    
#nl = nl[np.argsort(nd)]
#rij = rij[np.argsort(nd),:]
#nd = nd[np.argsort(nd)]
#                                    
### chop off at 3rd NN
#third = np.argwhere(nd > np.unique(nd)[3]).min()
#nl = nl[0:third]
#nd = nd[0:third]
#rij = rij[0:third,:]
#
#rh = rij/2

###################################
if plot == 'yes':
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    rgb = np.array([0.25,0.25,1]).reshape(1,3) 
    #has to be like this or matplotlib bitches
    ax1.scatter(rpos[:,0],rpos[:,1],rpos[:,2],s=10,c=rgb,edgecolors='k')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_aspect('equal')
    ax1.set_title('Direct Lattice')
    plt.show()
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    rgb = np.array([1,0.25,0.25]).reshape(1,3) 
    #has to be like this or matplotlib bitches
    ax2.scatter(kpos[:,0],kpos[:,1],kpos[:,2],s=10,c=rgb,edgecolors='k')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_aspect('equal')
    ax2.set_title('Reciprocal Lattice')
    plt.show()
    





