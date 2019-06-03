#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 12:41:27 2019

@author: ty
"""

#def makeGaN(nx,ny,nz,lammps='no'):
#""" 
#Same as make FCC, returns all the same stuff.
#"""
import numpy as np
import copy as cp

nx, ny, nz = [2,2,2]
lammps = 'data.relax'

masses = np.array([69.723,14.007]) 
a = 3.189 #a = b = diagonal length
c = np.round(np.sqrt(8/3.0)*a,decimals=3) #5.185/2.0 #ideal c = sqrt(8/3)*a

basis = np.array(([1,0,0,0], #Ga
                  [2,0,0,3], #N
                  [1,2,0,4], #Ga
                  [2,2,0,7], #N
                  [1,3,1,0], #Ga
                  [2,3,1,3], #N
                  [1,5,1,4], #Ga
                  [2,5,1,7])).astype(float) #N
uc = np.array([0,0,0,0,0,0,0,0])
ids = np.arange(0,8)

pos = cp.deepcopy(basis)
tmp = cp.deepcopy(basis)
tuc = cp.deepcopy(uc)
tmpids = cp.deepcopy(ids)
for i in range(nx-1):
    tmp[:,1] = tmp[:,1]+(6)
    pos = np.append(pos,tmp,axis=0)
    tuc = tuc+1
    uc = np.append(uc,tuc)
    ids = np.append(ids,tmpids,axis=0)
    
tmp = cp.deepcopy(pos)
tuc = cp.deepcopy(uc)
tmpids = cp.deepcopy(ids)
for i in range(ny-1):
    tmp[:,2] = tmp[:,2]+2
    tuc = tuc+nx
    uc = np.append(uc,tuc)
    pos = np.append(pos,tmp,axis=0)
    ids = np.append(ids,tmpids,axis=0)
    
tmp = cp.deepcopy(pos)
tuc = cp.deepcopy(uc)
tmpids = cp.deepcopy(ids)
for i in range(nz-1):
    tmp[:,3] = tmp[:,3]+8
    pos = np.append(pos,tmp,axis=0)
    tuc = tuc+nx*ny
    uc = np.append(uc,tuc)
    ids = np.append(ids,tmpids,axis=0)
    
pos[:,1] = pos[:,1]*np.sqrt(3)*a/6.0
pos[:,1] = np.round(pos[:,1],decimals=4)
pos[:,2] = pos[:,2]*a/2.0
pos[:,2] = np.round(pos[:,2],decimals=4)
pos[:,3] = pos[:,3]*c/8.0
pos[:,3] = np.round(pos[:,3],decimals=4)

num = len(pos)
pos = np.append(np.arange(1,num+1).reshape(num,1),pos,axis=1)

if lammps != 'no':    
    with open(lammps,'w') as fid:
        xbuff = 0.4603
        ybuff = 0.79725
        zbuff = 0.3255
        fid.write('LAMMPS DATA FILE FOR SED\n')
        fid.write('\n'+str(num)+' atoms\n')
        fid.write('\n'+str(len(masses))+' atom types\n')
        fid.write('\n'+str(pos[:,2].min()-xbuff)+' '+
                  str(pos[:,2].max()+xbuff)+' xlo xhi\n')
        fid.write(str(pos[:,3].min()-ybuff)+' '+str(pos[:,3].max()+ybuff)+
                  ' ylo yhi\n')
        fid.write(str(pos[:,4].min()-zbuff)+' '+str(pos[:,4].max()+zbuff)+
                  ' zlo zhi\n')
        fid.write('\nMasses\n\n')
        for i in range(len(masses)):
            fid.write(str(i+1)+' '+str(masses[i])+'\n')
        fid.write('\nAtoms\n\n')
        for i in range(num-1):
            fid.write(str(int(pos[i,0]))+' '+str(int(pos[i,1]))+' '
                      +str(pos[i,2])+' '+str(pos[i,3])+' '+
                      str(pos[i,4])+'\n')
        fid.write(str(int(pos[-1,0]))+' '+str(int(pos[-1,1]))+' '
                  +str(pos[-1,2])+' '+str(pos[-1,3])+' '+str(pos[-1,4]))

#return [num, pos, masses, uc, ids, a, c]  
