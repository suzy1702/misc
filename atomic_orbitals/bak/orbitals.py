#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 18:30:23 2019
@author: ty
"""
import numpy as np
from math import factorial as fact
from scipy.special import binom
import sys

class atomic_orbital:
   """
   Given 3 quantum numbers, radial and angular grid size, build Hydrogren
   atom atomic wavefunctions.
   """
   def __init__(self,n,l,m):
      """
      Error check the quantum number input and save the values.
      """
      if (not isinstance(n,int) or not isinstance(l,int) or not isinstance(m,int)):
         sys.exit('\tUSAGE ERROR: all quantum numbers must be integers!')
      if n <= 0:
         sys.exit('USAGE ERROR: principal quantum number n must be n >= 1!')
      if l < 0 or l > n-1:
         sys.exit('\tUSAGE ERROR: angular momentum number l must but 0 <= l <= n-1!')
      if m < -l or m > l:
         sys.exit('USAGE ERROR: magnetic quantum number m must be -l <= m <= l!')
         
      self.q_num = {'n':n,'l':l,'m':m}
      self.wave_func = {'Radial':{},'Angular':{}}
      
   def radial(self,dr=0.001,rmax=10,units='Bohr'):
      """
      Construct the radial part of Hydrogen atom wavefunctions
      
      Default radial grid step size is 0.001, max radius is 10 Bohr
      
      Default radial units are Bohr, set units='Angstrom' for Angstroms.
      
      Spherical polar coordinates:
         - x = r.sin(theta).cos(phi)
         - y = r.sin(theta).sin(phi)
         - z = r.cos(theta)
      """
      n = self.q_num['n']
      l = self.q_num['l']
      
      if (units != 'Bohr' and units != 'Angstrom' and 
          units != 'bohr' and units != 'angstrom'):
         print('\tUSAGE ERROR: radial units \"{}\" not recognized, using Bohr '
               'instead!'.format(units))
         units = 'Bohr'
         
      if units == 'Bohr' or units == 'bohr':
         a_0 = 1 # Bohr
      else:
         a_0 = 0.5291772 # Angstrom
         
      # radial grid
      r = np.arange(0,rmax/a_0,dr)
      nr = len(r)
         
      # David Miller, Quantum Mechanics, eq. (10.65)
      # generalized Laguerre Polynomials
      lpoly = np.zeros(nr)
      p = n-l-1
      j = 2*l+1
      for q in range(p+1):
         lpoly[:] = lpoly[:]+binom(p+j,p-q)*(-1)**q*(2*r/n/a_0)**q/fact(q)
            
      # David Miller, Quantum Mechanics, eq. (10.72)
      big_r = np.zeros(nr)
      norm = np.sqrt(fact(n-l-1)/(2*n)/fact(n+l)*(2/n/a_0)**3) 
      # normalization coefficient
      big_r = norm*(2*r/n/a_0)**l*lpoly*np.exp(-r/n/a_0)
      
      self.wave_func['Radial'] = {'R(r)':big_r,'r':r,'Units':units}
      self.lpoly = lpoly
      
   def angular(self,dtheta=0.001,dphi=0.001):
      """
      Construct the angular part of the Hydrogen atom wavefunctions
      
      Default step size for azimuthal and polar grid is 0.001 
      
      phi E [0,2*pi), theta E [0,pi)
            
      Note that the Spherical harmonic function is given as a complex NxM 
      matrix where the N rows are the theta values and the M columns are the
      phi values.
            
      Spherical polar coordinates:
         - x = r.sin(theta).cos(phi)
         - y = r.sin(theta).sin(phi)
         - z = r.cos(theta)
      """
      l = self.q_num['l']
      m = self.q_num['m']
      
      # azimuthal angle grid (between r and positive x-axis)
      phi = np.arange(0,2*np.pi,dphi)
      nphi = len(phi)
      # polar angle grid (between r and positive z-axis)
      theta = np.arange(0,np.pi,dtheta)
      ntheta = len(theta)
      
      # David Miller, Quantum Mechanics, eq. (9.31)
      # associated Legendre functions
      # according to wikipedia, l and m are required to be non-negative integers
      # <https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Definition ...
      # _for_non-negative_integer_parameters_%E2%84%93_and_m>
      lfunc = np.zeros(ntheta)
      x = np.cos(theta)
      coeff = (-1)**abs(m)*2**l*(1-x**2)**(abs(m)/2)
      for k in range(abs(m),l+1):
         lfunc[:] = lfunc[:]+(coeff * fact(k)/fact(k-abs(m)) * x**(k-abs(m)) * 
              binom(l,k) * binom(0.5*(l+k-1),l))
      
      # David Miller, Quantum Mechanics, eq. (9.34)
      # here in the normalization part, l and m have usual rules.
      norm = (-1)**m*np.sqrt((2*l+1)/(4*np.pi)*fact(l-abs(m))/fact(l+abs(m))) 
      spharm = np.array([lfunc]*nphi).T
      spharm = spharm*np.exp(1j*m*np.array([phi]*ntheta))*norm
      # theta == rows, phi == columns
      self.wave_func['Angular'] = {'Spherical Harmonics':spharm,
                    'theta':theta,'phi':phi,'cos(theta)':x}
      
      