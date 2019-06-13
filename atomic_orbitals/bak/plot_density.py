#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:38:04 2019

@author: ty
"""
import matplotlib.pyplot as plt
import numpy as np
from orbital import atomic_orbital


# Create wave functions from atomic_orbital class
# enter different n, l, m numbers. If they aren't physically reasonable,
# it will bitch at you

wf1 = atomic_orbital(8,4,0) # n, l, m quantum numbers

# for l constant and n -> big, rmax -> big 
# for n constant and l -> big, rmax -> big
# m doesn't do anything to the radial part.
wf1.radial(rmax=150,dr=0.1,units='Bohr') # radial part of the Hydrogen wavefunctions
wf1.angular() # angular part, see docstring for options


####################################################################
###### you don't have to change anything below here ################
####################################################################
# Stuff to plot prob. density, still shitty and doesn't include angular part 
# in the density(yet)
ang = wf1.wave_func['Angular']
rad = wf1.wave_func['Radial']
big_r = rad['R(r)']
r = rad['r']
rmax = r.max()
nr = len(r)

ng = 250 # size of grid edge
dens = np.zeros((ng,ng))
edgepos = np.arange(-ng/2,ng/2,1)
edgepos = edgepos/(ng/2.1)*rmax # normalize grid values so prob. fits on grid

for i in range(ng):
   for j in range(ng):
      dist = np.sqrt(edgepos[i]**2+edgepos[j]**2)
      if dist >= rmax:
         dens[i,j] = 0
      else:
         dens[i,j] = big_r[np.argwhere(r <= dist).max()]**2
         # obviously this doens't include the angular part

# Doesn't include angular part yet. The spherical harmonics are calculated by
# the .angular() method, I just haven't included the angular part in the density yet.

# To do (to include angular part):
# calculate the azimuthal angle for each point on the real space grid.
# the basic idea is write a method that takes the azimuthal and polar angles
# as arguments and makes a slice through probability density space along a 
# constant theta-phi plane. The simplest case is theta = 0, i.e the x-y plane

fig, ax = plt.subplots()
ax.imshow(dens,interpolation='none',cmap='inferno',aspect='equal')
ax.set_title('|n,l,m > = |{},{},{} >'.format(wf1.q_num['n'],wf1.q_num['l'],wf1.q_num['m']))
ax.set_axis_off()
fig.tight_layout()
plt.show()
