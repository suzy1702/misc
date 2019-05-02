#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 12:04:35 2019

@author: ty
"""
import numpy as np
#import sys

def getVol(datain):
    """
    Read the log.lammps output file ('login') for a line starting with
    'Simulation Cell', e.g. > Simulation Cell = 10320.5028715619 A^3,
    then copy the number to a variable 'vol' which is the simulation cell
    volume
    """
    fid = open(datain,'r')
    for line in fid:
        if 'xhi' in line:
            xcoord = np.array(line.strip().split()[0:2]).astype(float)
        if 'yhi' in line:
            ycoord = np.array(line.strip().split()[0:2]).astype(float)
        if 'zhi' in line:
            zcoord = np.array(line.strip().split()[0:2]).astype(float)
            break
    vol = ((xcoord[1]-xcoord[0])*(ycoord[1]-ycoord[0])
    *(zcoord[1]-zcoord[0])) #Ang^3
    return vol

def getLines(infile):
    """
    This function returns the number of lines in a file
    """
    nlines = sum(1 for lines in open(infile,'r'))
    return nlines

def trap(inarr,step):
    """
    This function returns the antiderivative of an array using trapezoidal 
    integration along axis 0. The input array is integrated along each column.
    The input array 'inarr' is integrated with step size 'step' and assigned  
    to 'outarr'. outarr is returned.
    """
    outarr = np.zeros((len(inarr[:,0]),len(inarr[0,:])))

    for i in range(len(inarr[0,:])):
        for j in range(len(inarr[:,0])-1):
            outarr[j+1,i] = (outarr[j,i] + (inarr[j,i]+inarr[j+1,i])*step/2.0)
            
    return outarr

                
            
        