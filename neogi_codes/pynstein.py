#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 12:04:50 2019

@author: ty

Calculate thermal conductivity from Einstein relation (Green-Kubo MD)

This code applies to a 3D solid with input file format:
        <time(ps) JX JY JZ(ev*A/ps)>
"""
import numpy as np
import pynsteinMod as pm
import matplotlib.pyplot as plt
import scipy.integrate as sci

##### INPUT PARAMETERS #######
hfluxin = 'HEATFLUX_new' #heatcurrent is in ev*A/ps
datain = 'data.vels' 

dtMD = 0.001 #MD timestep in ps 
dn = 5 #number of MD steps between writing to file
dt = dtMD*dn #effective timestep of heatflux data
temp = 300 #Kelvin

nlines = pm.getLines(hfluxin)-1 #number of data points minus comment line
vol = pm.getVol(datain) #Ang^3

split = 1 #number of times to split a single trajectory for averaging
ncorr = nlines/split #length of the data fit in number of points
tTot = dt*ncorr #total time of one trajectory block in ps
time = np.linspace(0,tTot-dt,ncorr) #time array

kb = 8.6173324e-5 #eV/K
const = 1/(kb*temp**2*vol) #1/(eV*K*Ang^3)
#ev = 1.60217646e-19 #eV to J
##############################

for i in range(1):#split): #debugging
    jxyz = np.zeros((ncorr,3))
    with open(hfluxin,'r') as fid:
        fid.readline() #skip header
        for j in range(ncorr): #read the currents for this block
            jxyz[j,:] = fid.readline().strip().split()[1:5]
            
        avgj = np.mean(jxyz,axis=0)
        stdj = np.std(jxyz,axis=0,ddof=1) #divisor is (N-1)
        jxyz = jxyz-avgj #why do this? 
        
        antij = pm.trap(jxyz,dt) #antiderivative of jxyz



