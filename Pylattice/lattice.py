#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A code that takes real space lattice vectors and creates the k space
lattice and the Wigner-Seitz cells in both spaces. 

rvec
"""
import numpy as np
import latticeFX as fx
import matplotlib.pyplot as plt

a = 2.510


rvec = np.array([[   a,             0,    0], #primitive triang. lattice
                 [-a/2,     3**.5*a/2,    0],
                 [   0,             0,    1]]) 
rvol = fx.comp_vol(rvec)
kvec, kvol = fx.make_kvec(rvec,rvol)
n1, n2, n3 = [5 ,5 ,0] #should all be odd integer numbers
rpos = fx.make_lattice(rvec,n1,n2,n3)
kpos = fx.make_lattice(kvec,n1,n2,n3)



rvec2 = np.array([[a,            0,     0],     #ortho conve. cell
                 [ 0,      a*3**.5,     0],
                 [ 0,            0,     1]]) 
rvol2 = fx.comp_vol(rvec2)
kvec2, kvol2 = fx.make_kvec(rvec2,rvol2)
rpos2 = fx.make_lattice(rvec2,n1,n2,n3)
kpos2 = fx.make_lattice(kvec2,n1,n2,n3)



fig1, ax1 = plt.subplots()
rgb = np.array([0,0,1]).reshape(1,3) 
ax1.scatter(rpos[:,0],rpos[:,1],s=20,c=rgb,edgecolors='k')
rgb = np.array([1,0,0]).reshape(1,3) 
ax1.scatter(rpos2[:,0],rpos2[:,1],s=5,c=rgb,edgecolors='k')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
#ax1.set_aspect('equal')
ax1.set_title('Direct Lattice')
plt.show()
plt.savefig('directlattice.png',format='png',dpi=300)



fig2, ax2 = plt.subplots()
rgb = np.array([0,0,1]).reshape(1,3) 
ax2.scatter(kpos[:,0],kpos[:,1],s=20,c=rgb,edgecolors='k')
rgb = np.array([1,0,0]).reshape(1,3) 
ax2.scatter(kpos2[:,0],kpos2[:,1],s=5,c=rgb,edgecolors='k')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_aspect('equal')
ax2.set_title('Reciprocal Lattice')
plt.show()
plt.savefig('reciplattice.png',format='png',dpi=300)




rgb = np.array([1,0,0]).reshape(1,3) 
#has to be like this or matplotlib bitches
ax1.scatter(rpos[:,0],rpos[:,1],s=10,c=rgb,edgecolors='k')
plt.savefig('directlattice.png',format='png',dpi=300)

rgb = np.array([1,0,0]).reshape(1,3) 
#has to be like this or matplotlib bitches
ax2.scatter(kpos[:,0],kpos[:,1],s=10,c=rgb,edgecolors='k')
#    







