#!/usr/bin/env python
#This program creates the average and standard deviation variation of thermal conductivity as a function of correlation time over number of independent runs -- helps to confirm the fit  

import numpy as np
import sys
import math

# Check if all the input is provided 
if len(sys.argv ) !=5 :
    sys.exit ( "Usage: executable column_to_average numruns firstrun diffusionfilesuffix" )
# input file
# input_xvg = open(sys.argv[1])
# column to average
column = int(sys.argv[1])
# number of runs
nruns = int(sys.argv[2])
# first run
firstrun = int(sys.argv[3])
# name of file
diffusionfilesuffix = sys.argv[4]

a = firstrun
b = diffusionfilesuffix

for a in range(nruns):
    tmp = np.loadtxt('heat.' + str(a) + '/diffusion' + str(b) + '.dat')
    if a == firstrun:
        data = tmp[:,(column-1)]
    else:
        data = np.column_stack((data,tmp[:,(column-1)]))

kappamean = np.mean(data, axis=1)
kappastd = np.std(data, axis=1, dtype=np.float64, ddof=1)

#kappaxstd = np.std(kxdata, axis=None, dtype=np.float64,ddof=1) # divided by (N-1)

#kappaerror = kappastd/(math.sqrt(nruns-1))
kappacombi = np.column_stack((tmp[:,0],kappamean,kappastd))

# write kappa matrix
ka = open('run_kappa' +str(b) +'.dat', 'w')

for i in range(kappamean.size):
  kappacombi[i,:].tofile(ka, sep="    ", format="%e")
  ka.write("\n")
ka.close()

