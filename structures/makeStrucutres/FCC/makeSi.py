#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
An example script to generate silicon orthorhombic supercells for LAMMPS 
calculations. 

"""
import numpy as np
import copy as cp

lammpsfile = 'data.lammps' # name of LAMMPS data file that will be created, 
# can be anything you want

nx = 2 # size of supercell in x
ny = 2 # size of supercell in y
nz = 2 # size of supercell in z

a = 5.431 # lattice constant
mass = 28.0855 # mass of silicon in AMU
 
#########################################################################
#######  You don't have to change anything below here ###################
#########################################################################
basis = np.array([[0,0,0], 
                  [0,2,2],
                  [2,0,2],
                  [2,2,0],
                  [1,1,1],
                  [3,3,1],
                  [1,3,3],
                  [3,1,3]]).astype(float) #8 atom conventional 
                                           #FCC-diamond cell

# replicate in x, y, z
# use cp.deepcopy() because python will create a reference otherwise, 
# i.e. changing the copy will change the original. cp.deepcopy() prevents this!
pos = cp.deepcopy(basis) 
tmp = cp.deepcopy(pos)
for i in range(nx-1): # replicate unit cell along x 
   tmp[:,0] = tmp[:,0]+4
   pos = np.append(pos,tmp,axis=0)

tmp = cp.deepcopy(pos)
for i in range(ny-1): # replicate x-direction structure along y to form plane
   tmp[:,1] = tmp[:,1]+4
   pos = np.append(pos,tmp,axis=0)

tmp = cp.deepcopy(pos)
for i in range(nz-1): # replicate x-y plane along z to form 3-d structure
   tmp[:,2] = tmp[:,2]+4
   pos = np.append(pos,tmp,axis=0)

pos = pos[np.lexsort((pos[:,2],pos[:,1],pos[:,0]))] # sort by coord just because
num = len(pos[:,0]) # number of atoms
pos = np.append(np.ones((num,1)),pos,axis=1) # only 1 type so all 1's
pos = np.append(np.arange(1,num+1).reshape(num,1),pos,axis=1) # add ids
# LAMMPS requires IDS in data file
pos[:,2] = pos[:,2]*a/4 # rescale to the real lattice constant
pos[:,3] = pos[:,3]*a/4
pos[:,4] = pos[:,4]*a/4

buffer = a/8 # you have to add this to include proper spacing between periodic
# images. If you don't add buffer, the atoms will be on top of each other at 
# the box boundary. a/8 is the right amount, i.e. half the bond distance along
# the box edge directions
xmax = pos[:,2].max()+buffer
xmin = 0-buffer
ymax = pos[:,3].max()+buffer
ymin = 0-buffer
zmax = pos[:,4].max()+buffer
zmin = 0-buffer

# this section writes the LAMMPS file. It is simple but LAMMPS has particular
# requirements for the structure.
with open(lammpsfile, 'w') as fid:
   fid.write(str('LAMMPS position FILE\n'))
   fid.write('\n' + str(num) + ' atoms\n')
   fid.write('\n1 atom types\n')
   fid.write('\n{} {} xlo xhi'.format(xmin,xmax))
   fid.write('\n{} {} ylo yhi'.format(ymin,ymax))
   fid.write('\n{} {} zlo zhi'.format(zmin,zmax))
   fid.write('\n\nMasses\n')
   fid.write('\n1 ' + str(mass))
   fid.write('\n\nAtoms\n\n')
   for i in range(num-1):
       fid.write('{} {} {} {} {}\n'.format(
             int(i+1),int(pos[i,1]),pos[i,2],pos[i,3],pos[i,4]))
   fid.write('{} {} {} {} {}'.format(
         int(num),int(pos[i,-1]),pos[i,-2],pos[i,-3],pos[i,-4]))