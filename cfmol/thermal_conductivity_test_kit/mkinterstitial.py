#!/usr/bin/env python
#This program creates input configurations for bulk lattice with atoms diffused in the interstitial locations
#
#---Sanghamitra Neogi--------
#---Boulder, January, 2017---
#----------------------------

import sys
import os
import math
import numpy as np
import random

def randpoint_incube(side):
    """ Generate a random point inside a cube.
    """
    radius = side/2
#    x = random.random()*radius*2-radius
    x = random.random()*radius*2
    y = random.random()*radius*2
    z = random.random()*radius*2
    return [x,y,z]

#configFileName = input('Enter configuration filename(in .xyz):')
concR = float(input('Enter interstitials concentration(x/100%):'))
#conc = float(input('Enter interstitials concentration for upper limit(x/100%):'))
atomName2Substitute = input('Enter atom name to be substituted:')
interstitialAtomName = input('Enter interstitial atom name:')
#xLowCutoff = float(input('Enter lower x boundary of alloy region):'))
#xUpCutoff = float(input('Enter upper x boundary of alloy region):'))
#yLowCutoff = float(input('Enter lower y boundary of alloy region):'))
#yUpCutoff = float(input('Enter upper y boundary of alloy region):'))
#zLowCutoff = float(input('Enter lower z boundary of alloy region):'))
#zUpCutoff = float(input('Enter upper z boundary of alloy region):'))

print('Config file name is X.xyz') 

configFileName = 'X.xyz'
xLowCutoff = 5
xUpCutoff = 80
yLowCutoff = 5
yUpCutoff = 80
zLowCutoff = 5
zUpCutoff = 80
cubeBoundary = 86.8960 
interstitialCount = 0

f = open(configFileName,'r')
coordinates = f.readlines()
f.close()

print(coordinates[0])
atomstoInsert = float(coordinates[0])*float(concR)
print(atomstoInsert)

alloyfile = open("Xwithinterstitials.xyz", 'w')


for numlines in range(0,len(coordinates)):
	line = coordinates[numlines].split()
	if (numlines < 2):
		print(coordinates[numlines])
		linestowrite = [coordinates[numlines]]
#		linestowrite = [line[0], "    ", line[1], "    ", line[2], "    ", line[3], "\n"]
		alloyfile.writelines(linestowrite)
	else:
		atomtype = line[0]
		linestowrite = [atomtype, "    ", line[1], "    ", line[2], "    ", line[3], "\n"]
		alloyfile.writelines(linestowrite)

atomCount = 0
while (interstitialCount < atomstoInsert):
#while ((atomCount <= int(coordinates[0])) and (interstitialCount <= atomstoInsert)):
	randPosition = randpoint_incube(cubeBoundary)
	print(randPosition, int(coordinates[0]), interstitialCount)
	if ( (randPosition[0] > 1.0) and (randPosition[1] > 1.0) and (randPosition[2] > 1.0) and (randPosition[0] < xUpCutoff) and (randPosition[1] < yUpCutoff) and (randPosition[2] < zUpCutoff) ):
		for numlines in range(2,len(coordinates)):
			line = coordinates[numlines].split()
			if ( ( abs(float(line[1])-randPosition[0]) < 1.0) and (abs(float(line[2])-randPosition[1]) < 1.0) and (abs(float(line[3])-randPosition[2]) < 1.0)):
				linestowrite = [interstitialAtomName, "    ", str(float(line[1])+1.5), "    ", str(float(line[2])+1.5), "    ", str(float(line[3])), "\n"]
#				linestowrite = [interstitialAtomName, "    ", str(randPosition[0]), "    ", str(randPosition[1]), "    ", str(randPosition[2]), "\n"]
				alloyfile.writelines(linestowrite)
				interstitialCount+= 1
#	atomCount+= 1

#alloyfile.close()

print(interstitialCount) 
print('Change atom numbers in Xwithinterstitials.xyz') 
