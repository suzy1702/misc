#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 10:12:53 2019

module containing functions for calculating transmission

@author: Ty Sterling <ty.sterling@colorado.edu>
"""

def readHarmonicKij(forcefile,tol=0.0):
    '''
    This function reads a file, 'forcefile', of lammps force data calculated
    by finite displacements. Some assumptions are made about the structure of 
    the lammps output file. The first line should be 'du ###' where ### is the
    size of the displacement step. The lammps code loops over all atoms, num,
    by displacing them first in the positive x direction then the negative 
    x and then the positive y, negative y, etc. It does this for each atom.
    Each time an atom is displaced, the forces acting on all atoms is printed
    to the file. The difference in forces in each direction, e.g. (F_+x)-(F_-x)
    is calculated and interpreted as the 'finite difference derivative' dF. 
    The force constant d^2V_ij/du_i/du_j is the second order term of the
    taylor expansion of the potential energy calculated as d^2V_ij/du_i/du_j =
    -dF_i/du_j. dF_i is the difference in force on each atom du to displacement
    of atom j in the x,y,z directions. NOTE: its acutally -dF/(2*du) since 
    the derivative is taken about the equilibrium position in +/- minus x,y,z.
    
    The matrix fij this creates is formatted as follows:
        
               i1x      i1y      i1z       i2x ...
        ____|______________________________________
            |
        j1x |  -dFj1x   -dFj1x  -dFj1x   -dFj1x   
            |   ----     ----     ----     ----   ... 
            |  2*dui1x  2*dui1y  2*dui1z  2*dui2x
            |
            |
        j1y |  -dFj1y   -dFj1y   -dFj1y   -dFj1y  
            |   ----     ----     ----     ----   ...
            |  2*dui1x  2*dui1y  2*dui1z  2*dui2x
            |
            |
        j1z |  -dFj1z   -dFj1z   -dFj1z   -dFj1z   
            |   ----     ----     ----     ----   ...
            |  2*dui1x  2*dui1y  2*dui1z  2*dui2x
            |
            |
        j2x |  -dFj2x  -dFj2x   -dFj2x    -dFj2x
          . |   ----    ----     ----      ----   ...
          . |  2*dui1x  2*dui1y  2*dui1z  2*dui2x
          . |     .       .         .        .
                  .       .         .        .
                  .       .         .        .
                  
    The optional argument 'tol' can be set to some small float value.
    if the magnitude of the value in the matrix is lower than tol, it is set
    to 0. This could be used to sparsify the matrix but I dont care...

    this function return: n, nl, and nr = the total number of atoms, the number 
    in the left, and the number in the right groups respectively. idsL, idsR,
    and idsfij correspons to the ids in each group. kii is the force constant
    matrix between atoms in the left region and atoms in the left region - not 
    used in transmission calculation. kij is the force constant matrix between
    atoms in the left region and atoms in the right region
    '''
    import numpy as np
    import sys
    with open(forcefile, 'r') as fid:

        nl = int(fid.readline().strip().split()[1]) #left atoms
        nr = int(fid.readline().strip().split()[1]) #right atoms
        n = nl+nr #total number of atoms
        kii = np.zeros((nl*3,nl*3)) #force constant matrix left w/ left
        kij = np.zeros((nl*3,nr*3)) #force constant matrix left w/ right
        idsL = np.zeros((1,1)) #ids of left atoms
        idsR = np.zeros((1,1)) #ids of right atoms
        idsfij = np.zeros((n,1)) #LAMMPS ids in force file
        
        for i in range(n):
            tmp = fid.readline().strip().split()
            idsfij[i] = int(tmp[0])
            if int(tmp[1]) == 1:
                idsL = np.append(idsL,i)
            else:
                idsR = np.append(idsR,i)
        idsL = np.delete(idsL,0,axis=0).astype(int)
        idsR = np.delete(idsR,0,axis=0).astype(int)
        
        if (nl != nr) or (len(idsL) != len(idsR)): #make sure NL == NR
            sys.exit('Number of atoms on left must equal number on right!')
        
        du = float(fid.readline().strip().split()[1]) #atomic displacement
        #for finite difference derivative
        
        for i in range(nl*3): #move left atoms (i) x, y, and z
            fplus = np.zeros((n,3)) #dFij from positive shift of atom i
            fminus = np.zeros((n,3)) #dFij from negative shift of atom i
            
            for j in range(9): #skip comments
                fid.readline()
            for k in range(n): #dF on all atoms, positive shift
                tmp = fid.readline().strip().split()
                fplus[k,:] = tmp[1:4]
            for j in range(9): #skip comments
                fid.readline()
            for k in range(n): #dF on all atoms, negative shift of atom i
                tmp = fid.readline().strip().split()
                fminus[k,:] = tmp[1:4]
                
            fii = np.subtract(fplus[idsL,:],fminus[idsL,:]) 
            fii = -(np.reshape(fii,nl*3)).transpose()/(2.0*du)
            #dF between atom i and atoms in left side
            fij = np.subtract(fplus[idsR,:],fminus[idsR,:]) 
            fij = -(np.reshape(fij,nl*3)).transpose()/(2.0*du)
            #dF between atoms i and atoms in right ride
            if tol != 0.0:
                fii[np.argwhere(np.abs(fii) < tol)] = 0
                fij[np.argwhere(np.abs(fij) < tol)] = 0
                
            kii[i,:] = fii
            kij[i,:] = fij
        
    #Force constant matricies Kij = [d^2U/duiduj] = [-dFj/dui]
#    kii = -kii/(2*du) #finite difference derivative: dFplus-dFminus/2*du
#    kij = -kij/(2*du) 
    return [n, nl, nr, idsL, idsR, idsfij, kii, kij]

    
