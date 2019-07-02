#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:----------(type name of code here)

Created on Fri Jun 28 15:49:01 2019

Author:---------Ty Sterling
Contact:--------ty.sterling@colorado.edu
Instituion:-----University of Colorado Boulder
--------------------Raman Spectroscopy and Neutron Scattering Laboratory
--------------------Professor Dmitry Reznik

Description: You have to go down to where 'masses' is introduced and enter the
masses of the atoms corresponding to the atoms in 'types'
"""

import numpy as np
from diffpy.Structure import loadStructure
from diffpy.Structure.expansion import supercell

infile = 'C_mp-667273_conventional_standard.cif'
outfile = 'c60'

c60 = loadStructure(infile)
c60_2x2x2 = supercell(c60,[12,2,2])
c60.write(outfile+'.xyz','xyz')
c60_2x2x2.write(outfile+'_2x2x2.xyz','xyz')

struct = c60_2x2x2
with open(outfile+'_2x2x2.xyz','r') as fid:
      nat = int(fid.readline()) # number of atoms
      fid.readline() # skip comment line
      
      pos = np.zeros((nat,4))
      types = []
      for i in range(nat):
            tmp = fid.readline().strip().split()
            if tmp[0] not in types:
                  types.append(tmp[0])
                  pos[i,0] = types.index(tmp[0])+1
            else:
                  pos[i,0] = types.index(tmp[0])+1
                  
            pos[i,1:4] = tmp[1:4]

nat = len(pos[:,0])
masses = ['']*len(types)
masses[0] = 12.0107
cell = np.array([struct.lattice.a,struct.lattice.b,struct.lattice.c])
angles = np.array([struct.lattice.alpha,struct.lattice.beta,struct.lattice.gamma])
with open('data.'+outfile,'w') as fid:
      fid.write('LAMMPS input file\n')
      fid.write('\n{} atoms\n'.format(nat))
      fid.write('\n{} atom types\n'.format(len(types)))
      
      xmin = pos[:,1].min()
      xmax = pos[:,1].max()
      bx = (cell[0]-(xmax-xmin))/2
      
      ymin = pos[:,2].min()
      ymax = pos[:,2].max()
      by = (cell[1]-(ymax-ymin))/2
      
      zmin = pos[:,3].min()
      zmax = pos[:,3].max()
      bz = (cell[2]-(zmax-zmin))/2
      
      fid.write('\n{:.6f} {:.6f} xlo xhi\n'.format(xmin-bx,xmax+bx))
      fid.write('{:.6f} {:.6f} ylo yhi\n'.format(ymin-by,ymax+by))
      fid.write('{:.6f} {:.6f} zlo zhi\n'.format(zmin-bz,zmax+bz))

      fid.write('\nMasses\n\n')
      
      for i in range(len(types)):
            fid.write('{} {}\n'.format(i+1,masses[i]))
            
      fid.write('\nAtoms\n\n')
      for i in range(nat-1):
            fid.write('{} {} {:.6f} {:.6f} {:.6f}\n'.format(i+1,int(pos[i,0]),
                      pos[i,1],pos[i,2],pos[i,3]))
      fid.write('{} {} {:.6f} {:.6f} {:.6f}'.format(nat,int(pos[-1,0]),
                      pos[-1,1],pos[-1,2],pos[-1,3]))
      

                  
                  
                  
                  
                  
                  
                  
            
