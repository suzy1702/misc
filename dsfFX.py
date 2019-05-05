#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 13:15:13 2019

@author: ty
"""
import numpy as np
import time

def readData(filename):
    """
    Reads a LAMMPS data file and returns the no of atoms, no of types, masses 
    of each types, and an Nx5 array containing ids, type, x, y, z coords
    """
    nlines = sum(1 for line in open(filename,'r'))

    boxvec = np.zeros((3,2))
    with open(filename,'r') as fid:
        for i in range(nlines):
            tmp = fid.readline().strip().split()
            if len(tmp) > 1 and tmp[1] == 'atoms':
                natoms = int(tmp[0])
            if len(tmp) > 2 and tmp[1] == 'atom' and tmp[2] == 'types':
                ntypes = int(tmp[0])
                masses = np.zeros(ntypes)
            if len(tmp) > 3 and tmp[-1] == 'xhi':
                boxvec[0,:] = tmp[0:2]
                boxvec[1,:] = fid.readline().strip().split()[0:2]
                boxvec[2,:] = fid.readline().strip().split()[0:2]
            if len(tmp) > 0 and tmp[0] == 'Masses':
                fid.readline()
                for j in range(ntypes):
                    masses[j] = fid.readline().strip().split()[1]
            if len(tmp) > 0 and tmp[0] == 'Atoms':
                fid.readline()
                tmp = fid.readline().strip().split()
                break
            
        ncol = len(tmp)
        pos = np.zeros((natoms,ncol))
        pos[0,:] = tmp
        for i in range(1,natoms):
            pos[i,:] = fid.readline().strip().split()
      
    return [natoms, ntypes, masses, pos, boxvec]

########################################################################
def tic():
    """
    Same as MATLAB tic and toc functions. Use ty.tic() at the beginning of
    code you want to time and ty.toc() at the end. Once ty.toc() is reached,
    elapsted time will be printed to screen and optionally (by default) written
    to 'log.txt' file.
    """
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    """
    Same as MATLAB tic and toc functions. Use ty.tic() at the beginning of
    code you want to time and ty.toc() at the end. Once ty.toc() is reached,
    elapsted time will be printed to screen and optionally (by default) written
    to 'log.txt' file.
    """
    if 'startTime_for_tictoc' in globals():
        print(("\n\tElapsed time is "+
              str(np.round(time.time()-
                       startTime_for_tictoc,decimals=3))+" seconds."))
    else:
        print("\n\t\tToc: start time not set") 