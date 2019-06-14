#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:38:04 2019
@author: ty
"""
import matplotlib.pyplot as plt
from orbitals import atomic_orbital
from grids import grid

# Create wave functions from atomic_orbital class 
wf = atomic_orbital(3,2,0) # n, l, m quantum numbers
wf.radial(rmax=35,dr=0.0001,units='Bohr') # radial part of the wavefunction
wf.angular() # angular part, see docstring for options

wfgrid = grid(wf,nedge=1000)
wfgrid.density(plane='xz')

############################################################
############### plot the densities #########################
############################################################
color = 'hot'
fig, axs = plt.subplots(1,3)

images = ['']*3
images[0] = axs[0].imshow(abs(wfgrid.dens['total dens']**2),interpolation='none',
                                         cmap=color,aspect='equal')
images[1] = axs[1].imshow(abs(wfgrid.dens['radial dens'])**2,interpolation='none',
                                         cmap=color,aspect='equal')
images[2] = axs[2].imshow(abs(wfgrid.dens['angular dens'])**2,interpolation='none',
                                         cmap=color,aspect='equal')

fig.suptitle(r'$|{},{},{}\rangle$'.format(wfgrid.q_num['n'],
                                      wfgrid.q_num['l'],wfgrid.q_num['m']))
fig.set_size_inches(8,4,forward=True)

cbar = fig.colorbar(images[0],ax=axs,orientation='horizontal',fraction=0.1,
                    ticks=[0,(abs(wfgrid.dens['total dens'])**2).max()])
cbar.ax.set_xticklabels(['-','+'])

axs[0].set_title('Probability Density')
axs[1].set_title('Raidal Contribution')
axs[2].set_title('Angular Contribution')
for ax in axs:
   ax.set_axis_off()

plt.show()