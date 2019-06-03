#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 8 10:21:12 2019

@author: ty
"""
import numpy as np
import copy as cp

#nx = int(input('Enter number of unit cells in x:\n\t>'))
#ny = int(input('Enter numer of unit cells in y:\n\t>'))
#xyz = str(input('Enter the name of the xyz data file or "none" if '
#                       'you dont want one\n\t>'))
#lammps = str(input('Enter the name of the lammps data file or "none" if '
#                       'you dont want one\n\t>'))

nx = 100
ny = 10

xyz = 'example.xyz'
lammps = 'data.hbn'

ux = 1
uy = 1

## ORTHO UNIT CELL LATTICE VECTORS
ax = 2.510 #2*1.255
ay = 4.348 # 6*0.724574

## 4 ATOM CONVENTIONAL CELL -- BASIS
basis = np.array([[1,   1,    0.907,   0,   0,   0],   #B
                  [1,   2,   -0.907,   1,   1,   0],   #N
                  [1,   1,    0.907,   1,   3,   0],   #B
                  [1,   2,   -0.907,   0,   4,   0]])  #N
              #mol.id  #type  #charge   #x   #y   #z
            
masses = np.array([10.811,14.0067]) #masses in AMU

## REPLICATE TO FORM SUPERCELL
pos = cp.deepcopy(basis)
tmp = cp.deepcopy(basis)
for i in range(ux-1):
    tmp[:,3] = tmp[:,3]+2
    pos = np.append(pos,tmp,axis=0)
tmp = cp.deepcopy(pos)
for i in range(uy-1):
    tmp[:,4] = tmp[:,4]+6
    pos = np.append(pos,tmp,axis=0)

nb = len(pos[:,0])
uc = np.zeros(nb)
ids = np.arange(0,nb)
                
## REPLICATE TO FORM STRUCTURE
tmp = cp.deepcopy(pos)
tuc = cp.deepcopy(uc)
tids = cp.deepcopy(ids)
dx = pos[:,3].max()+1
for i in range(nx-1):
    tmp[:,3] = tmp[:,3]+dx
    pos = np.append(pos,tmp,axis=0)
    tuc = tuc+1
    uc = np.append(uc,tuc)
    ids = np.append(ids,tids)
    
tmp = cp.deepcopy(pos)
tuc = cp.deepcopy(uc)
duc = uc.max()+1
tids = cp.deepcopy(ids)
dy = pos[:,4].max()+2
for i in range(ny-1):
    tmp[:,4] = tmp[:,4]+dy
    pos = np.append(pos,tmp,axis=0)
    tuc = tuc+duc
    uc = np.append(uc,tuc)
    ids = np.append(ids,tids)
    
## ORGANIZE
num = len(pos[:,0])
pos = np.append(np.arange(1,num+1,1).reshape(num,1),pos,axis=1) #ids
pos = pos.astype(float)
pos[:,4] = pos[:,4]*ax/2 
pos[:,5] = pos[:,5]*ay/6

## DETERMINE BOND INFO
bx = ax/4 #buffer between periodic images
by = ay/6 # '' ''
xmin = pos[:,4].min()-bx
xmax = pos[:,4].max()+bx
ymin = pos[:,5].min()-by
ymax = pos[:,5].max()+by
zmin = pos[:,6].min()-1000
zmax = pos[:,6].max()+1000
    
lx = (xmax-xmin)/2 #distance to shift for pbc
ly = (ymax-ymin)/2

nl = np.zeros((num,4))
nd = np.zeros((num,4))
print('\n\tFinding nearest neighbors, this will take a while ...\n')
for i in range(num): #loop over particles
    if i != 0 and i%(num//10) == 0:
        print('\t\tNow '+str(np.round(10*i/(num//10),0))+
              ' % done finding neighbors.')
        
    rvec = pos[i,4:6]
    nnids = np.intersect1d(np.intersect1d(np.argwhere(pos[:,4] <= rvec[0]+1.5),
                           np.argwhere(pos[:,4] >= rvec[0]-1.5)),
                           np.intersect1d(np.argwhere(pos[:,5] <= rvec[1]+1.5),
                           np.argwhere(pos[:,5] >= rvec[1]-1.5)))
    nnids = np.append(nnids,np.argwhere(pos[:,4] <= xmin+1.5))
    nnids = np.append(nnids,np.argwhere(pos[:,4] >= xmax-1.5))
    nnids = np.append(nnids,np.argwhere(pos[:,5] <= ymin+1.5))
    nnids = np.append(nnids,np.argwhere(pos[:,5] >= ymax-1.5))
    nnids = np.unique(nnids)
    
    nn = len(nnids)
    tmp = np.zeros((2,nn)) # ids, distance
    for j in range(nn): #loop over neighbors
        rx = pos[nnids[j],4]-rvec[0] # vector from partice i to j
        ry = pos[nnids[j],5]-rvec[1] # vector from partice i to j

        ## Minimum image convention
        if rx >= lx:
            rx = rx-lx*2
        elif rx <= -lx:
            rx = rx+lx*2
        if ry >= ly:
            ry = ry-ly*2
        elif ry <= -ly:
            ry = ry+ly*2
        
        dist = np.sqrt(rx**2+ry**2)
        tmp[0,j] = pos[nnids[j],0]-1 #index
        tmp[1,j] = dist #distance
    
    nl[i,:] = tmp[0,np.argsort(tmp[1,:])[0:4]] #only keep the first few. 
    nd[i,:] = tmp[1,np.argsort(tmp[1,:])[0:4]] #saves memory.

print('\n\tDone finding neighbors!')
del ax, ay, basis, bx, by, dist, tmp, rx, ry, lx, ly

## BONDS
nbonds = num*3
bonds = np.zeros((nbonds,4))
d = 0
for i in range(num):
    for j in range(3):
        bonds[d,2] = i+1
        bonds[d,3] = nl[i,j+1]+1
        if (len(np.argwhere(bonds[:,2] == bonds[d,3])) != 0 and
            len(np.argwhere(bonds[:,3] == bonds[d,2])) != 0):
                bonds[d,0] = 1
        d = d+1

bonds = np.delete(bonds,np.argwhere(bonds[:,0] == 1),axis=0)    
nbonds = len(bonds[:,0])
bonds[:,0] = np.arange(1,nbonds+1)
bonds[:,1] = 1
        
## ANGLES
nangles = num*3
angles = np.zeros((nangles,5))
angles[:,0] = np.arange(1,nangles+1)
d = 0
for i in range(num):
    if pos[int(nl[i,0]),2] == 1:
        atype = 1
    else:
        atype = 2
    for j in range(3):
        angles[d,1] = atype
        angles[d,3] = i+1
        if j == 0:
            angles[d,2] = nl[i,j+1]+1
            angles[d,4] = nl[i,j+2]+1
        if j == 1:
            angles[d,2] = nl[i,j+1]+1
            angles[d,4] = nl[i,j+2]+1
        if j == 2:
            angles[d,2] = nl[i,j+1]+1
            angles[d,4] = nl[i,j-1]+1
        d = d+1
        
## IMPROPERS
nimprop = num
improp = np.zeros((num,6))
improp[:,0] = np.arange(1,nimprop+1)
for i in range(num):
    if pos[int(nl[i,0]),2] == 1:
        atype = 1
    else:
        atype = 2
    improp[i,1] = atype
    improp[i,2:6] = nl[i,0:4]+1

## WRITE XYZ  
if xyz != 'none':
    with open(xyz,'w') as fid:
        fid.write(str(num)+'\n')
        fid.write(str(xmin)+' '+str(xmax)+' '+str(ymin)+' '+str(ymax)+' '+
                  str(zmin)+' '+str(zmax)+'\n')
        for i in range(num-1):
            if pos[i,2] == 1:
                fid.write('B ')
            else: 
                fid.write('N ')
            fid.write(str(pos[i,4])+' '+str(pos[i,5])+' '+str(pos[i,6])+'\n')
        if pos[-1,2] == 1:
            fid.write('B ')
        else: 
            fid.write('N ')
        fid.write(str(pos[-1,4])+' '+str(pos[-1,5])+' '+str(pos[-1,6]))
        
## WRITE LAMMPS
if lammps != 'none':    
    with open(lammps,'w') as fid:
        fid.write('LAMMPS hBN\n')
        fid.write('\n'+str(num)+' atoms')
    
        fid.write('\n'+str(nbonds)+' bonds')
        fid.write('\n'+str(nangles)+' angles')
        fid.write('\n'+str(nimprop)+' impropers\n')
        
        fid.write('\n'+str(len(masses))+' atom types')
        fid.write('\n'+str(1)+' bond types')
        fid.write('\n'+str(2)+' angle types')
        fid.write('\n'+str(2)+' improper types\n')
        
        fid.write('\n'+str(xmin)+' '+str(xmax)+' xlo xhi\n')
        fid.write(str(ymin)+' '+str(ymax)+' ylo yhi\n')
        fid.write(str(zmin)+' '+str(zmax)+' zlo zhi\n')
        fid.write('\nMasses\n\n')
        
        for i in range(len(masses)):
            fid.write(str(i+1)+' '+str(masses[i])+'\n')
            
        fid.write('\nAtoms\n\n')
        for i in range(num):
            fid.write(str(int(pos[i,0]))+' '+str(int(pos[i,1]))+' '
                      +str(int(pos[i,2]))+' '+str(pos[i,3])+' '+
                      str(pos[i,4])+' '+str(pos[i,5])+' '+str(pos[i,6])+'\n')
        
        bonds = bonds.astype(int)
        fid.write('\nBonds\n\n')
        for i in range(nbonds):
            fid.write(str(bonds[i,0])+' '+str(bonds[i,1])+' '+str(bonds[i,2])+
                      ' '+str(bonds[i,3])+'\n')
            
        angles = angles.astype(int)
        fid.write('\nAngles\n\n')
        for i in range(nangles):
            fid.write(str(angles[i,0])+' '+str(angles[i,1])+' '+str(angles[i,2])+
                      ' '+str(angles[i,3])+' '+str(angles[i,4])+'\n')
            
        improp = improp.astype(int)
        fid.write('\nImpropers\n\n')
        for i in range(nimprop-1):
            fid.write(str(improp[i,0])+' '+str(improp[i,1])+' '+str(improp[i,2])+
                      ' '+str(improp[i,3])+' '+str(improp[i,4])+' '
                      +str(improp[i,5])+'\n')
        fid.write(str(improp[-1,0])+' '+str(improp[-1,1])+' '+str(improp[-1,2])+
                  ' '+str(improp[-1,3])+' '+str(improp[-1,4])+' '
                  +str(improp[-1,5]))