#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 09:47:48 2019

Author:---------Ty Sterling
Contact:--------ty.sterling@colorado.edu
Instituion:-----University of Colorado Boulder
--------------------Raman Spectroscopy and Neutron Scattering Laboratory
--------------------Professor Dmitry Reznik

Description: This script reads the output of lm-sensors over a time period
and prints the temperatures of the CPUs to check if they got too hot.
      
"""
import numpy as np
import matplotlib.pyplot as plt

nprocs = 16

class proc_data:
      def __init__(self,temps,cores,maxtemps,critical,dt):
            """
            Writes the real temperatures, cpus,
            """
            self.data = {'Temperature':temps,
                         'Cores':cores,
                         'Max Temps':maxtemps,
                         'Critical Temps':critical,
                         'Time Step':dt}

cores = np.arange(nprocs)
maxtemps = np.zeros(nprocs)
critical = np.zeros(nprocs)
nl = sum(1 for lines in open('results.txt','r')) # number of lines in file
blocks = (nl-2)//31
temps = np.zeros((nprocs,blocks))
with open('results.txt','r') as fid:
      dt = int(fid.readline()) # seconds between data block
      fid.readline()
      
      for i in range(blocks):
            for j in range(9): # skip bullcrap
                  fid.readline()
            for j in range(8): # read first processor
                  tmp = fid.readline().strip().split()
                  temps[j,i] = tmp[2][1:5]
                  if i == 0:
                        maxtemps[j] = tmp[5][1:5]
                        critical[j] = tmp[8][1:5]
            for j in range(4): # skip more bullcrap
                  fid.readline()
            for j in range(8): # read second processor
                  tmp = fid.readline().strip().split()
                  temps[8+j,i] = tmp[2][1:5]
                  if i == 0:
                        maxtemps[8+j] = tmp[5][1:5]
                        critical[8+j] = tmp[8][1:6]
            for j in range(2): # skip more bullcrap
                  fid.readline() 
                        
results = proc_data(temps,cores,maxtemps,critical,dt) # save to class just cause
fig, ax = plt.subplots() # plot the temps
for i in range(blocks):
      d = (i+1)/blocks
      ax.scatter(cores,temps[:,i],c=[[0.0,0.0,0.0+d]],label=None)
ax.plot(cores,np.ones(nprocs)*maxtemps,'r',label='Max design temp')
ax.plot(cores,np.ones(nprocs)*critical,'k',label='Critical design temp')
ax.axis((-0.5,nprocs-0.5,temps.min()-5,critical.max()+5))
fig.legend()
ax.set_ylabel('Temperature (\'C)')
ax.set_xlabel('Cores')

plt.show()

