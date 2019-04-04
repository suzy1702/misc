#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 10:39:39 2018

@author: ty
"""
import copy as cp
import numpy as np

infile = 'config.xyz'
etc = 'patterns'

dl = 4 #number of monolayers to shift in each direction TRANSVERSE to interface
#dw = 2 #number of monolayers to shift in each direction PARALLEL to interface
##if 0, pattern is at the center of the structure in z, - includes more Si, + includes
##more ge

#run it and see

with open(infile, 'r') as fid:
    
    size = int(fid.readline()) #get system size
    total = cp.deepcopy(size)
    fid.readline() #skip comments
    
    data = np.zeros((size,5))
    for i in range(size):
        tmp = fid.readline().strip().split() #read the coordinates
        
        data[i,0] = (i+1) #id
        if i == 0:
            atom_type = tmp[0]
        if tmp[0] == 'Si':
            data[i,1] = 1 #type
        else:
            data[i,1] = 2
        data[i,2] = float(tmp[1]) #x coord
        data[i,3] = float(tmp[2]) #y coord
        data[i,4] = float(tmp[3]) #z coord
        
###########
size = len(data)

si = data[:size//2,:]
ge = data[size//2:,:]

six = np.sort(np.unique(si[:,2]))[-dl:]
gex = np.sort(np.unique(ge[:,2]))[:dl]

siids = np.argwhere(data[:,2] == six[0])
geids = np.argwhere(data[:,2] == gex[0])
for i in range(1,dl):
    if dl == 1:
        break
    siids = np.append(siids,np.argwhere(data[:,2] == six[i]))
    geids = np.append(geids,np.argwhere(data[:,2] == gex[i]))
    
z = np.sort(np.unique(data[:,4]))
bottom = z[:len(z)//2]
top = z[len(z)//2:]
    
topids = np.argwhere(data[:,4] == top[0])
botids = np.argwhere(data[:,4] == bottom[0])
for i in range(1,len(z)//2):   
    topids = np.append(topids,np.argwhere(data[:,4] == top[i]))
    botids = np.append(botids,np.argwhere(data[:,4] == bottom[i]))
    
siids = np.intersect1d(siids,topids)
geids = np.intersect1d(geids,botids)
    
data[siids,1] = 2
data[geids,1] = 1


###find max and min coords          
maxcoords = data.max(axis=0)
mincoords = data.min(axis=0)
xmax = maxcoords[2]
ymax = maxcoords[3]
zmax = maxcoords[4]
xmin = mincoords[2]
ymin = mincoords[3]
zmin = mincoords[4]
del maxcoords, mincoords


###########  WRITE TO .xyz FILE    ##################
filename = etc + '.xyz'
with open(filename, 'w') as f:
    f.write(str(len(data))+'\n')
    f.write(str(xmax)+'\t0'+'\t0'+'\t0\t'+str(ymax)+'\t0'+'\t0'+'\t0\t'+str(zmax)+'\n')
    for i in range(len(data)-1):
        if data[i][1] == 1:
            f.write('Si\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
        if data[i][1] == 2:
            f.write('Ge\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
    if data[-1][1] == 1:
        f.write('Si\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
    if data[-1][1] == 2:
        f.write('Ge\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
####################################################



###########  WRITE LAMMPS DATA  ###################
xbuff = 5.431/8
yzbuff = 5.431/8

masses = ['28.08556', '72.6400']
 
datafile = 'data.' + etc
with open(datafile, 'w') as f:
    
    f.write(str('LAMMPS DATA FILE\n'))
    
    f.write('\n' + str(len(data)) + ' atoms\n')
    f.write('\n' + str(len(masses)) + ' atom types\n')
    f.write('\n' + str(xmin-xbuff)+' '+str(xmax+xbuff)+' xlo'+' xhi\n')
    f.write(str(ymin-yzbuff)+' '+str(ymax+yzbuff)+' ylo'+' yhi\n')
    f.write(str(zmin-yzbuff)+' '+str(zmax+yzbuff)+' zlo'+' zhi\n')
    f.write('\nMasses\n')
    for i in range(len(masses)):
        f.write('\n' + str(i+1) + ' ' + str(float(masses[i])))
    f.write('\n\nAtoms\n\n')
    for i in range(len(data)-1):
        f.write(str(int(i+1)) + ' ' + str(int(data[i,1])) + ' ' + str(data[i,2]) + ' ' +
                str(data[i,3]) + ' ' + str(data[i,4]) + '\n')
    f.write(str(len(data)) +  ' ' + str(int(data[-1,1])) + ' ' + str(data[-1,2]) + ' ' + 
            str(data[-1,3]) + ' ' + str(data[-1,4]))
##################################################
