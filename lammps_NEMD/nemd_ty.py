#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 12:22:53 2019

@author: ty

Rewrite of Neogi's NEMD code aimed at learning OOP.
Instead of asking the user to enter the name of the lammps output files,
it is automatically assumed that the lammps thermo data file is named 
'log.lammps' and that the temperature data is named 'tmp.profile.0'. If you 
don't like this functionality, modify the code!

This code reads the FIRST BLOCK of heatflux data from a LAMMPS simulution 
according to my LAMMPS code included in the same directory as this. Note that
it looks for the header of the thermo data and finds the keyword 'v_fluxin' to
find the start of the data. If you change the name of the variable or the number
or order of columns in the thermo data in the LAMMPS script, you'll have to 
modify the 'read lammps data' part of this code. It should be easy.

Once the heatflux data is read in, it is printed to the screen and you have to 
pick the 'stationary state' region. I.e. pick the part that is flat. You enter
the beginning and end of the stationary region into the prompt in units of steps.
Note that is acceptable to type 1e6 instead of 1000000. It then asks you to 
divide up the data into block for averaging. The limits on the number of blocks
is put by how many blocks of temperature profile data there is. It may not work
if you enter less that 3 because of the way I plot it. If it crashes, you can
fix it by just averaging using a few more blocks... or go to the part labeled
'clunky plotting part' and rewrite it to work better.

Next, it reads in the temperature profile and plots it. You then have to pick
the linear regions on either side of the interface and identify the center of 
the interface. Enter the points you choose as 5 integers seperated by spaces.
If you enter anything else, it will complain. Also note that what you enter has
to make sense, i.e. the rightmost point can't be on the left of the leftmost...
You should pick the regions on either side of the interface 

The code then takes over. It computes the average stationary state flux, then 
divides it by the temperature gradient (or discontinuity at the interface) to 
determine thermal transport properties. Uncertainty propogation is tracked so
error bars are included in the results!
"""
import numpy as np
import sys
import matplotlib.pyplot as plt

# the LAMMPS script was supposed to 1/Area in angstrom to meters as 
# 1/1e-10/1e-10. I botched it and accidentally did 1/1e-1/1e-1 this 'oops' 
# factor retroactively corrects it during this script! If it gets fixed in 
# the LAMMPS script later, just set oops = 1!
oops = 0.1*0.1/1e-20

## structure to store all the data
class nemd_data:
   """
   This class defines a structure containing the NEMD data. 
   
   hflux contains all the raw heatflux data, fits, and stats and temp contains
   all the temperature data, fits, and stats. 
   
   The methods are used to get the bounds for fitting.
   """
   def __init__(self,hflux,hfstep,temp,tstep,tcoord):
      self.hflux = {'steps':hfstep,'flux':hflux}
      self.temp = {'temperature':temp,'steps':tstep,'coords':tcoord}
      
   def get_bounds(self): # get bounds of data to keep
      """
      Find the stationary state region of the heatflux. Follow the prompts and
      pick the flat region!
      """
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
      
      # trim the data and overwrite it
      self.hflux['flux'] = self.hflux['flux'][hf_lbound:hf_ubound,:]
      self.hflux['steps'] = self.hflux['steps'][hf_lbound:hf_ubound]
      self.temp['temperature'] = self.temp['temperature'][t_lbound:t_ubound,:]
      self.temp['steps'] = self.temp['steps'][t_lbound:t_ubound]
      
      # average the heatflux across the bounded region
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
      # save the fit stats to the structure for later use
      self.hflux['flux_stats'] = flux_stats
           
   # pick to good regions from the temp profile for linear fitting
   def get_temp_bins(self):
      """
      Again, follow the prompts to pick the good regions of data.
      
      For the temperature profile, you should pick the linear regions on either
      side of the interface and identify the center of the interface.
      
      Acceptable usage is a string of 5 ascending integers seperated by spaces.
      """
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
      
      # average across the blocks so that the fitting is done on better data
      temp = self.temp['temperature']
      binsize = len(temp[:,0])//self.numavg
      std = np.zeros((self.numavg,len(temp[0,:])))
      avg = np.zeros((self.numavg,len(temp[0,:])))
      for i in range(self.numavg):
         tmptemp = temp[i*binsize:(i+1)*binsize,:]
         std[i,:] = np.std(tmptemp,dtype=np.float64,axis=0)
         avg[i,:] = np.mean(tmptemp,dtype=np.float64,axis=0)
      self.temp['bins'] = {'std':std,'avg':avg}
      
      
      # the rest of this function does the linear fitting and averaging.
      # a plot of all the different fits to the different blocks is shown to 
      # verify that the results are as expected.
      bins = self.temp['bins']
      pts = self.grad_points
      coords = self.temp['coords']
      fits = {'lslope':['']*self.numavg,'rslope':['']*self.numavg,
              'deltaT':['']*self.numavg}
      
      ########### CLUNKY PLOTTING PART #######################
      # this logic to divide up the number of plots is a little crappy. If you
      # get crashes, rewrite this part. Its pretty obvious what its plotting.
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
      ########## END OF THE CLUNKY PLOTTING PART ################
      
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



##############################################################################
## MAIN FUNCTION ##
##############################################################################
#
#               
## get the heatflux data from log.lammps
nlines = sum(1 for lines in open('log.lammps','r'))   
with open('log.lammps','r') as fid:
   for i in range(nlines):
      # find the location and size of the heatflux block, only the first one
      # i.e. if there is a second NEMD portion for transmission, just ignore it
      tmp = fid.readline().strip().split()
      if len(tmp) == 13 and tmp[-1] == 'v_fluxin': # column header
         fpos = fid.tell() 
         break
         
   ### FIND HEATFLUX DATA ###
   j = 0
   while j != -1: # get the size of the block
      if j >= 1e12: # loop will run forever if conditions fail so this will 
         # break it eventually....
         sys.exit('\tERROR: seems to be stuck in the while loop\n\tApparently '
                  'the code can\'t find the end of the heatflux data block!\n\t'
                  'Look at the part of the code labeled \'find heatflux data\'!')
      j = j+1
      tmp = fid.readline().strip().split()
      try:
         float(tmp[-1])
      except:
         nhflux = j-1
         break
            
   # go back to the start of the block and read in the data
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
    
   
## create an instance of the class and get all the fit data
## *oops is included because I accidentally convereted angstrom to meters
## as *1e-1 instead of 1e-10 ...
data = nemd_data(hflux[:,-2:]*oops,hfstep,temp,tstep,tcoord)
data.get_bounds() # get heatflux-converged region from user input
data.get_temp_bins() # get temperature profile fits from user input


# compute kappa for each side and G across the interface
# I copied the error propogation part from the code we used in my old group
# (thanks Dr. Sanghamitra Neogi and Paul Salame!) but here is the reference 
# they cave:
#----------Combine two measurements to get composite kappa value---------
# Source: http://www.burtonsys.com/climate/composite_standard_deviations.html

flux_stats = data.hflux['flux_stats']
fit_stats = data.temp['fit_stats']
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
sumoferror_of_sumofsquares_l = (k_l1_std**2)*(data.numavg-1)+(k_l2_std**2)*(data.numavg-1)
overall_group_sumofsquares_l = (((kappa_l-kappa_l1)**2)*data.numavg+
                                ((kappa_l-kappa_l2)**2)*data.numavg)
kappa_l_std = np.sqrt((sumoferror_of_sumofsquares_l+
                       overall_group_sumofsquares_l)/(2*data.numavg-1))

kappa_r = (kappa_r1+kappa_r2)/2
sumoferror_of_sumofsquares_r = (k_r1_std**2)*(data.numavg-1)+(k_r2_std**2)*(data.numavg-1)
overall_group_sumofsquares_r = (((kappa_r-kappa_r1)**2)*data.numavg+
                                ((kappa_r-kappa_r2)**2)*data.numavg)
kappa_r_std = np.sqrt((sumoferror_of_sumofsquares_r+
                       overall_group_sumofsquares_r)/(2*data.numavg-1))

g_1 = abs(flux_stats['avg']['in']/fit_stats['deltaT']['avg'])/1e6 
g_1_std = g_1*np.sqrt((flux_stats['std']['in']/flux_stats['avg']['in'])**2+
                      (fit_stats['deltaT']['std']/fit_stats['deltaT']['avg'])**2)

g_2 = abs(flux_stats['avg']['out']/fit_stats['deltaT']['avg'])/1e6 
g_2_std = g_2*np.sqrt((flux_stats['std']['out']/flux_stats['avg']['out'])**2+
                      (fit_stats['deltaT']['std']/fit_stats['deltaT']['avg'])**2)
   
g = (g_1+g_2)/2
sumoferror_of_sumofsquares_g = (g_1_std**2)*(data.numavg-1)+(g_2_std**2)*(data.numavg-1)
overall_group_sumofsquares_g = ((g-g_1)**2)*data.numavg+((g-g_2)**2)*data.numavg
g_std = np.sqrt((sumoferror_of_sumofsquares_g+
                 overall_group_sumofsquares_g)/(2*data.numavg-1))


# the calculation is over, just write the data!
# print the results to the screen and to a text file
print('Left Side: {} +/- {} (W/m/K)'.format(kappa_l,kappa_l_std))
print('Right Side: {} +/- {} (W/m/K)'.format(kappa_r,kappa_r_std))
print('Interface: {} +/- {} (W/m/m/K)'.format(g,g_std))

with open('thermal_prop.txt','w') as fid:
   fid.write('Thermal conductivity in the left and right sides of a heterostructure '
             'and conductance across the interface!\n\n')
   fid.write('! Left Side: {} +/- {} (W/m/K)\n'.format(kappa_l,kappa_l_std))
   fid.write('! Right Side: {} +/- {} (W/m/K)\n'.format(kappa_r,kappa_r_std))
   fid.write('! Interface: {} +/- {} (W/m/m/K)\n'.format(g,g_std))
   
