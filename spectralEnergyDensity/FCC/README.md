# Contents
computeSED-FCC.py: python code that calculates phonon spectral energy density from EMD velocity data. This code is applicable to FCC diamond solids only, e.g. C, Si, Ge, and zincblende structures such as GaAs, SiC, GaN, etc. You must include ty.py in your working directory for this code to work. 

in.relax and in.vels: LAMMPS MD codes to relax the structure and collect velocites for FCC diamond sytems with a triclinic simulation box. Use ty.makeTriclinic(n1,n2,n3,lammps='data.relax',element='c') to create a LAMMPS input configuration file for, e.g. carbon, with a triclinic simulation box. 

