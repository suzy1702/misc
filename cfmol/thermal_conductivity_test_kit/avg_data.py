#!/usr/bin/env python
#This program computes average and standard deviation of a column of data 

import numpy as np
import sys
import math

# Check if all the input is provided 
if len(sys.argv ) !=2 :
    sys.exit ( "Usage: executable inputfile" )

# input file
input_file = open(sys.argv[1])
tmp = np.loadtxt(input_file)

data = tmp
print len(data)

datamean = np.mean(data, axis=None)
#datastd = np.std(data, axis=None, dtype=np.float64)  #divides by N
datastd = np.std(data, axis=None, dtype=np.float64,ddof=1) # divided by (N-1)
dataerror = datastd/np.sqrt(len(data))

print 'data', datamean, '+/-', dataerror, 'std',datastd

