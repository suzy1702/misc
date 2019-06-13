#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 16:08:15 2019

@author: ty
"""
import numpy as np

class grid:
   def __init__(self,wf,nedge=1000):
      """
      Takes a wavefunction object, wf, as an angurment
      
      nedge is the size of the grids to be plotted, i.e. number of elements 
      on an edge. Default is 1000.
      """
      # build radial and angular grids to map wavefunction components to
      edgepos = np.linspace(-1,1,nedge)
      xgrid = np.array([edgepos]*nedge)
      ygrid = np.copy(np.fliplr(xgrid).T)
      zgrid = np.copy(ygrid)
      xgrid = xgrid+(((xgrid == 0).astype(int))*((xgrid != 0).astype(int)+1e-9))
      # set all the pos. where xgrid == 0 to finite value to avoid divide by 0
      
      self.grids = {'xgrid':xgrid,'ygrid':ygrid,'zgrid':zgrid,'nedge':nedge}
      self.q_num = {'n':wf.q_num['n'],'l':wf.q_num['l'],'m':wf.q_num['m']}
      self.angular = wf.wave_func['Angular']
      self.radial = wf.wave_func['Radial']
      
   def density(self,plane='xz'):
      """
      Enter the plane to compute angular density from
      
      Acceptable usage is 'xz', 'zx', 'xy', or 'yx',
      'XZ', 'ZX', 'XY', or 'YX'
      """
      if plane == 'xz' or plane == 'XZ' or plane == 'zx' or plane == 'ZX':
         plane = 'xz'
      elif plane == 'xy' or plane == 'XY' or plane == 'yx' or plane == 'YX':
         plane = 'xy'
      else:
         print('\tUSAGE ERROR: plane \"{}\" not recognized, using xz plane')
         plane = 'xz'
         
      spharm = self.angular['Spherical Harmonics']
      phi = self.angular['phi']
      nphi = len(phi)
      theta = self.angular['theta']
      ntheta = len(theta)
      
      big_r = self.radial['R(r)']
      r = self.radial['r']
      nr=len(r)
         
      # radial component of wavefunciton
      rgrid = np.sqrt(self.grids['xgrid']**2+self.grids['ygrid']**2)*nr
      rgrid = rgrid.astype(int) # indicies for radial part of wavefunction
      # mask radius beyond rad of wavefunction
      mask = (rgrid <= nr-1).astype(int) # beyond radius of wave function
      rdens = rgrid*mask # set values beyond to 0
      rdens = big_r[rdens]*mask # 0's are now the value of big_r[0], so set to 0 again   
      
      # phi component of spherical harmonic
      phigrid = np.arctan2(self.grids['ygrid'],self.grids['xgrid']) # x-y plane, i.e. theta = 0
      # arctan2 returns angles E [0,+/-pi), I want phi E (0,2*pi]
      phigrid[self.grids['nedge']//2:,:] = phigrid[self.grids['nedge']//2:,:]+2*np.pi 
      phigrid = (phigrid/(2*np.pi)*nphi).astype(int) # noralize to indicies of phi
      
      # theta component of spherical harmonic
      thetagrid = (np.arccos(self.grids['zgrid']/np.sqrt(self.grids['xgrid']**2
                             +self.grids['ygrid']**2))/np.pi*ntheta).astype(int)
      
      self.grids['rgrid'] = rgrid
      self.grids['phigrid'] = phigrid
      self.grids['thetagrid'] = thetagrid
      
      self.dens = {'radial dens':rdens}
      if plane == 'xz':
         self.dens['angular dens'] = spharm[thetagrid,0]
      else:
         self.dens['angular dens'] = spharm[0,phigrid]
      self.dens['total dens'] = self.dens['radial dens']*self.dens['angular dens']
      
         
      