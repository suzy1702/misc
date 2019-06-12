#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 12:22:53 2019

@author: ty

Rewrite of Neogi's NEMD code aimed at learning OOP.

"""
import datetime
import os
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
#plt.switch_backend('TKagg')
import scipy.stats

oops = 0.1*0.1/1e-20

## structure to store all the data
class nemd_data:
   def __init__(self,hflux,hfstep,temp,tstep,tcoord):
      self.hflux = {'steps':hfstep,'flux':hflux}
      self.temp = {'temperature':temp,'steps':tstep,'coords':tcoord}
      
   def get_bounds(self): # get bounds of data to keep
      print('\n\tPick the bounds of the heatflux data to compute transport '
            'from: (i.e. pick the flat region)')
      fig1, ax1 = plt.subplots()
      ax1.plot(self.hflux['steps'],self.hflux['flux'][:,0],'r',
               self.hflux['steps'],abs(self.hflux['flux'][:,1]),'b')
      ax1.set(xlabel='Step',ylabel='Flux, (W/m**2)',title='LAMMPS Heatflux')
      ax1.legend(('flux in','abs(flux out)'))
      ax1.axes.grid(True,linestyle=':',alpha=0.5,c=(0.5,0.5,0.5))
      plt.show()
      
      # get the input
      lowerstep = float(input('\tEnter the begenning of the bounds '
                              '(i.e. the step):\t'))
      upperstep = float(input('\tEnter the end of the bounds:\t'))

      # make sure input is within the bounds of the actual data
      if upperstep <= lowerstep:
         sys.exit('\n\tERROR: upper bound must be larger than lower bound!')
      if lowerstep < self.hflux['steps'].min():
         sys.exit('\n\tERROR: lower bound must be greater than {}!'.
                  format(self.hflux['steps'].min()))
      if upperstep > self.hflux['steps'].max():
         sys.exit('\n\tERROR: upper bound must be less than {}!'.
                  format(self.hflux['steps'].max()))
         
      # bounds for heatflux data determined by input
      hf_lbound = min(np.argwhere(self.hflux['steps'] >= lowerstep))[0]
      hf_ubound = max(np.argwhere(self.hflux['steps'] <= upperstep))[0]
      
      # bounds for temp data determined by input
      t_lbound = min(np.argwhere(self.temp['steps'] >= lowerstep))[0]
      t_ubound = max(np.argwhere(self.temp['steps'] <= upperstep))[0]

      # number of blocks determined by number of times temp is printed
      self.numavg = int(input('\tEnter the number of blocks to average '
                              '(must be between 1 and {}):\t'.
                              format(t_ubound-t_lbound)))
      
#      ### DELETE THIS BLOCK FOR FINAL SCRIPT
#      # bounds for heatflux data determined by input
#      hf_lbound = 1800 
#      hf_ubound = 2500
#      # bounds for temp data determined by input
#      t_lbound = 89
#      t_ubound = 124
#      # number of blocks determined by number of times temp is printed
#      self.numavg = 8
      
      ### NOT THIS ONE :)
      # trim the data and overwrite it
      self.hflux['flux'] = self.hflux['flux'][hf_lbound:hf_ubound,:]
      self.hflux['steps'] = self.hflux['steps'][hf_lbound:hf_ubound]
      self.temp['temperature'] = self.temp['temperature'][t_lbound:t_ubound,:]
      self.temp['steps'] = self.temp['steps'][t_lbound:t_ubound]
      
      hflux = self.hflux['flux']
      flux_stats = {'std':{},'avg':{}}
      avgin = np.mean(hflux[:,0],dtype=np.float64) 
      avgout = np.mean(-hflux[:,1],dtype=np.float64) 
      stdin = np.std(hflux[:,0],dtype=np.float64) 
      stdout = np.std(-hflux[:,1],dtype=np.float64) 
      flux_stats['std']['in'] = stdin
      flux_stats['std']['out'] = stdout
      flux_stats['avg']['in'] = avgin
      flux_stats['avg']['out'] = avgout
         
      self.hflux['flux_stats'] = flux_stats
           
   # pick to good regions from the temp profile for linear fitting
   def get_temp_bins(self):
      fig2, ax2 = plt.subplots()
      ax2.plot(np.arange(1,len(self.temp['temperature'][0,:])+1),
               np.mean(self.temp['temperature'],axis=0),'ko--',markersize=5)
      ax2.set(xlabel='Index (arb.)',ylabel='Temp (K)',
              title='LAMMPS Temperature Profile')
      ax2.axes.grid(True, linestyle=':', alpha=0.5, c=(0.5,0.5,0.5))
      plt.show()         
      
      self.grad_points = np.array(input('\n\tEnter the bounds to fit the temperature '
         'profile as follows:\n\tleftmost point on left, rightmost point on '
         'the left, the center of the interface,\n\tthe leftmost point on the '
         'right, and the rightmost point on the right\n\tUSAGE: enter 5 '
         'integers seperated by spaces (note that numbering starts at 1):\n\t').
         strip().split()).astype(int)-1
      
      if len(self.grad_points) != 5:
         sys.exit('\tERROR: there must be exactly 5 boundary points for '
                  'temperature fitting')
      if (self.grad_points[0] >= self.grad_points[1] or self.grad_points[1] >= 
          self.grad_points[2] or self.grad_points[2] >= self.grad_points[3] or
          self.grad_points[3] >= self.grad_points[4]):
         sys.exit('\tERROR: the 5 boundary points must be ascending and none '
                  'of them can be equal to another')
      
#      ### DELETE THIS LINE FOR FINAL CODE! 
#      self.grad_points = np.array([10, 17, 20, 22, 28])
      
      temp = self.temp['temperature']
      binsize = len(temp[:,0])//self.numavg
      std = np.zeros((self.numavg,len(temp[0,:])))
      avg = np.zeros((self.numavg,len(temp[0,:])))
      for i in range(self.numavg):
         tmptemp = temp[i*binsize:(i+1)*binsize,:]
         std[i,:] = np.std(tmptemp,dtype=np.float64,axis=0)
         avg[i,:] = np.mean(tmptemp,dtype=np.float64,axis=0)
      self.temp['bins'] = {'std':std,'avg':avg}
      
   def get_temp_fits(self):
      bins = self.temp['bins']
      pts = self.grad_points
      coords = self.temp['coords']
      fits = {'lslope':['']*self.numavg,'rslope':['']*self.numavg,
              'deltaT':['']*self.numavg}
      
      if self.numavg%2 != 0:
         nplots = self.numavg+1
      else:
         nplots = self.numavg
      t_fig, t_ax = plt.subplots(nplots//2,2)
      
      for i in range(self.numavg):
         fits['lslope'][i] = np.polyfit(coords[pts[0]:pts[1]],
             bins['avg'][i,pts[0]:pts[1]],1)
         fits['rslope'][i] = np.polyfit(coords[pts[3]:pts[4]],
             bins['avg'][i,pts[3]:pts[4]],1)
         fits['deltaT'][i] = np.polyval(fits['lslope'][i],
             coords[pts[0]:pts[2]])[-1]-np.polyval(fits['rslope'][i],
                   coords[pts[2]:pts[4]])[0]
      
         if i > (nplots//2-1):
            x = i-nplots//2; y = 1
         else: 
            x = i; y = 0
            
         t_ax[x,y].errorbar(coords[:-1],bins['avg'][i,:-1],
             xerr=None,yerr=bins['std'][i,:-1],color='k',marker='o',markersize=2) 
         t_ax[x,y].plot(coords[pts[0]:pts[2]],np.polyval(fits['lslope'][i],
               coords[pts[0]:pts[2]]),'r',
               coords[pts[2]:pts[4]],np.polyval(fits['rslope'][i],
               coords[pts[2]:pts[4]]),'r',zorder=100)
         t_ax[x,y].set_title('fit no. {}'.format(i+1))
         t_ax[x,0].set_ylabel('Temp, (K)')
         t_ax[-1,y].set_xlabel('Coord, (Angstrom)')
         t_fig.tight_layout()
         
      plt.show()
      fit_stats = {'lslope':{},'rslope':{},'deltaT':{}}
      fit_stats['lslope']['std'] = np.std(np.array(fits['lslope'])[:,0],
               dtype=np.float64) 
      fit_stats['lslope']['avg'] = np.mean(np.array(fits['lslope'])[:,0],
               dtype=np.float64) 
      fit_stats['rslope']['std'] = np.std(np.array(fits['rslope'])[:,0],
               dtype=np.float64) 
      fit_stats['rslope']['avg'] = np.mean(np.array(fits['rslope'])[:,0],
               dtype=np.float64) 
      fit_stats['deltaT']['std'] = np.std(np.array(fits['deltaT'])[:],
               dtype=np.float64) 
      fit_stats['deltaT']['avg'] = np.mean(np.array(fits['deltaT'])[:],
               dtype=np.float64) 
      self.temp['fit_stats'] = fit_stats
      self.temp['fits'] = fits

   def compute_kappa(self):
      flux_stats = self.hflux['flux_stats']
      fit_stats = self.temp['fit_stats']
      kappa_l1 = abs(flux_stats['avg']['in']/fit_stats['lslope']['avg']*1e-10)
      k_l1_std = kappa_l1*np.sqrt((flux_stats['std']['in']/flux_stats['avg']['in'])**2+
                                  (fit_stats['lslope']['std']/fit_stats['lslope']['avg'])**2)      
      
      kappa_l2 = abs(flux_stats['avg']['out']/fit_stats['lslope']['avg']*1e-10)
      k_l2_std = kappa_l2*np.sqrt((flux_stats['std']['out']/flux_stats['avg']['out'])**2+
                                  (fit_stats['lslope']['std']/fit_stats['lslope']['avg'])**2)
      
      kappa_r1 = abs(flux_stats['avg']['in']/fit_stats['rslope']['avg']*1e-10)
      k_r1_std = kappa_r1*np.sqrt((flux_stats['std']['in']/flux_stats['avg']['in'])**2+
                                  (fit_stats['rslope']['std']/fit_stats['rslope']['avg'])**2)
      
      kappa_r2 = abs(flux_stats['avg']['out']/fit_stats['rslope']['avg']*1e-10)
      k_r2_std = kappa_r1*np.sqrt((flux_stats['std']['out']/flux_stats['avg']['out'])**2+
                                  (fit_stats['rslope']['std']/fit_stats['rslope']['avg'])**2)
      
      kappa_l = (kappa_l1+kappa_l2)/2
      sumoferror_of_sumofsquares_l = (k_l1_std**2)*(self.numavg-1)+(k_l2_std**2)*(self.numavg-1)
      overall_group_sumofsquares_l = ((kappa_l-kappa_l1)**2)*self.numavg+((kappa_l-kappa_l2)**2)*self.numavg
      kappa_l_std = np.sqrt((sumoferror_of_sumofsquares_l+overall_group_sumofsquares_l)/(2*self.numavg-1))
      
      kappa_r = (kappa_r1+kappa_r2)/2
      sumoferror_of_sumofsquares_r = (k_r1_std**2)*(self.numavg-1)+(k_r2_std**2)*(self.numavg-1)
      overall_group_sumofsquares_r = ((kappa_r-kappa_r1)**2)*self.numavg+((kappa_r-kappa_r2)**2)*self.numavg
      kappa_r_std = np.sqrt((sumoferror_of_sumofsquares_r+overall_group_sumofsquares_r)/(2*self.numavg-1))
      
      g_1 = abs(flux_stats['avg']['in']/fit_stats['deltaT']['avg'])/1e6 
      g_1_std = g_1*np.sqrt((flux_stats['std']['in']/flux_stats['avg']['in'])**2+
                            (fit_stats['deltaT']['std']/fit_stats['deltaT']['avg'])**2)
      
      g_2 = abs(flux_stats['avg']['out']/fit_stats['deltaT']['avg'])/1e6 
      g_2_std = g_2*np.sqrt((flux_stats['std']['out']/flux_stats['avg']['out'])**2+
                            (fit_stats['deltaT']['std']/fit_stats['deltaT']['avg'])**2)
         
      g = (g_1+g_2)/2
      sumoferror_of_sumofsquares_g = (g_1_std**2)*(self.numavg-1)+(g_2_std**2)*(self.numavg-1)
      overall_group_sumofsquares_g = ((g-g_1)**2)*self.numavg+((g-g_2)**2)*self.numavg
      g_std = np.sqrt((sumoferror_of_sumofsquares_g+overall_group_sumofsquares_g)/(2*self.numavg-1))
      
      print('Left Side: {} +/- {} (W/m/K)'.format(kappa_l,kappa_l_std))
      print('Right Side: {} +/- {} (W/m/K)'.format(kappa_r,kappa_r_std))
      print('Interface: {} +/- {} (W/m/m/K)'.format(g,g_std))
               
## get the heatflux data from log.lammps
nlines = sum(1 for lines in open('log.lammps','r'))   
with open('log.lammps','r') as fid:
   for i in range(nlines):
      # find the location and size of the heatflux block, only the first one
      tmp = fid.readline().strip().split()
      if len(tmp) == 13 and tmp[-1] == 'v_fluxin': # column header
         fpos = fid.tell()
         break
         
   j = 0
   while j != -1: # get the size of the block
      if j >= 1e12:
         sys.exit('\tERROR: seems to be stuck in the while loop...')
      j = j+1
      tmp = fid.readline().strip().split()
      try:
         float(tmp[-1])
      except:
         nhflux = j-1
         break
            
   # go back to the start of the block
   hflux = np.zeros((nhflux,13))
   fid.seek(fpos)
   for i in range(nhflux):
      hflux[i,:] = fid.readline().strip().split()
   hfstep = hflux[:,0]
      
## get the temperature data from tmp.profile.0
nlines  = sum(1 for line in open('tmp.profile.0','r'))
with open('tmp.profile.0','r') as fid:
   for i in range(4): # skip comments
      tmp = fid.readline().strip().split()
   fid.seek(0)
   ndat = int(tmp[1]) # number of coordinates with temps in data
   nblocks = (nlines-3)//(ndat+1) # number of blocks (steps) printed
   tstep = np.zeros(nblocks)
   temp = np.zeros((nblocks,ndat))
   tcoord = np.zeros(ndat)
   for i in range(3): 
      fid.readline()
   for i in range(nblocks): # loop over all blocks
      tstep[i] = int(fid.readline().strip().split()[0])
      for j in range(ndat): # loop over block
          tmp = fid.readline().strip().split()
          if i == 0:
             tcoord[j] = float(tmp[1]) # coords same for all blocks
          temp[i,j] = float(tmp[-1]) # temp at each coord
    
      
data = nemd_data(hflux[:,-2:]*oops,hfstep,temp,tstep,tcoord)
data.get_bounds() # get heatflux-converged region from user input
data.get_temp_bins() # get temperature profile fits from user input
data.get_temp_fits() # fit the slope in each region
data.compute_kappa()
