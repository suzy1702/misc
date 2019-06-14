#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:38:04 2019
@author: ty
"""
import numpy as np
import matplotlib.pyplot as plt
from orbitals import atomic_orbital
from grids import grid

###############################################################################
############################### user input ####################################
###############################################################################

## enter an arbitrary number of states to be superposed here. The required 
## usage is a list of tuples (or lists!) each containing (n, l, m) quantum nums.
#qnums = [(2,1,0),
#         (3,2,1),
#         (4,3,2)]

## or just enter a single state like this (same usage, only one list element)
qnums = [(4,3,2)]

## these parameters should be the same for all states
rmax = 50
dr = 0.0001
nedge = 1000
plane = 'xz'



###############################################################################
########################## build the wavefunctions ############################
###############################################################################
nwfx = len(qnums)
wfxs = ['']*nwfx

for i in range(nwfx):
   # quantum state determined by n, l, and m quantum numbers
   # compute radial and angular part of wavefunctions
   wfxs[i] = atomic_orbital(qnums[i][0],qnums[i][1],qnums[i][2])
   wfxs[i].radial(rmax=rmax,dr=dr,units='Bohr')
   wfxs[i].angular(dtheta=0.001,dphi=0.001)
   
   # map radial and angular parts to a grid for plotting
   wfxs[i].grid = grid(wfxs[i],nedge=nedge)
   wfxs[i].grid.wave_func(plane=plane)

tot = np.zeros((nedge,nedge))
rad = np.zeros((nedge,nedge))
ang = np.zeros((nedge,nedge))
for j in range(nwfx):
   tot = tot+wfxs[j].grid.wavefunc['total wavefunction']
   rad = rad+wfxs[j].grid.wavefunc['radial wavefunction']
   ang = ang+wfxs[j].grid.wavefunc['angular wavefunction']
   
# ignore this. I was using it to set values == 0 to a finite number so the 
# logarithm is defined, i.e. to plot np.log(abs(tot)**2+mask)
mask = ((tot > 1e-6).astype(int)+1e-24)*(tot <= 1e-6).astype(int)
   


###############################################################################
########################### plot the densities ################################
###############################################################################
color = 'jet'
fig, axs = plt.subplots(1,3)

images = ['']*3
images[0] = axs[0].imshow(abs(tot)**2,interpolation='none',
                                         cmap=color,aspect='equal')
images[1] = axs[1].imshow(abs(rad)**2,interpolation='none',
                                         cmap=color,aspect='equal')
images[2] = axs[2].imshow(abs(ang)**2,interpolation='none',
                                         cmap=color,aspect='equal')

fig.set_size_inches(8,4,forward=True)

cbar = fig.colorbar(images[0],ax=axs,orientation='horizontal',fraction=0.1,
                    ticks=[0,(abs(tot)**2).max()])
cbar.ax.set_xticklabels(['-','+'])

axs[0].set_title('Probability Density')
axs[1].set_title('Raidal Contribution')
axs[2].set_title('Angular Contribution')
for ax in axs:
   ax.set_axis_off()

plt.show()