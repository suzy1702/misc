#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dynamic Structure Factor.

@author: ty
"""
import numpy as np
import dsfFX as fx
import matplotlib.pyplot as plt

fx.tic()

posfile = 'pos.dat'
datafile = 'data.vels'

split = 16 #times to split data for averaging
steps = 2**20 #total run time
dt = 1e-15 #lammps time step
dn = 2**3 #print frequency
prints = steps//dn #total times data is printed
tn = prints//split #times data is printed per chunk
thz = np.arange(0,tn)/(tn*dt*dn)*1e-12 #frequency in THz

natoms, ntypes, masses, pos, boxvec = fx.readData('data.vels') 
pos = np.delete(pos,1,axis=1)
pos = np.delete(pos,2,axis=1)

lx = boxvec[0,1]-boxvec[0,0]
ly = boxvec[1,1]-boxvec[1,0]
lz = boxvec[2,1]-boxvec[2,0]

nx, ny, nz = [24,0,0] #number of points to scan in each direction
kx = np.zeros((nx,3))
ky = np.zeros((ny,3))
kz = np.zeros((nz,3))

kx[:,0] = np.linspace(0,2*np.pi/lx,nx)

with open(posfile, 'r') as fid: 
    for i in range(1): #split): #loop over chunks to block average
        print('\n\tNow on chunk: '+str(i+1)+
              ' out of '+str(split)+'\n')

#        dat = np.zeros((tn,natoms,3))
#        print('\t\tNow reading positions...\n')
#        for j in range(tn): #loop over timesteps in block
#            if j!= 0 and j%(tn//10) == 0:
#                print('\t\t'+str(j/(tn//10)*10)+'% done reading '
#                        'velocites')
#            for k in range(9): #skip comments
#                fid.readline()
#            for k in range(natoms): #get atoms
#                tmp = fid.readline().strip().split()
#                dat[j,k,0] = float(tmp[2]) #vx
#                dat[j,k,1] = float(tmp[3]) #vy
#                dat[j,k,2] = float(tmp[4]) #vz
        
        isf = np.zeros((tn,nx)).astype(complex) #intermediate scat. fx
#        vh = np.zeros((tn,natoms,nx)).astype(complex) #Van Hove fx
#        for k in range(nx):
#            kvec = np.tile(kx[k,:],(tn,natoms,1))
#            vh[:,:,k] = np.exp(-1j*(np.multiply(dat[:,:,0],kvec[:,:,0])+
#                                    np.multiply(dat[:,:,1],kvec[:,:,1])+
#                                    np.multiply(dat[:,:,2],kvec[:,:,2])))
#        fx.toc()
        
#        for k in range(nx):
#            isf[:,k] = np.multiply(vh[:,:,k],np.conj(vh[:,:,k])).sum(axis=1)
        
                               
#        for k in range(nx):
#            for j in range(tn):
##                isf[j,k] = np.dot(vh[j,:,k],np.conj(vh[0,:,k]))
#                isf[j,k] = np.dot(np.conj(vh[j,:,k]),(vh[0,:,k]).T)
    
                
                
                
        
        
      
            
        


