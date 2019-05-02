#!/bin/bash

# Calculate VACF from trajectory
awk '{if(NF==4&& $1+0==$1) print $0}' traj.100ps.25fs |vacf.fft.v03 4096 4001 25
mv vacf.dat vacf.100ps.25fs.dat

# Calculate FT[VACF]
width=15
awk '{print $1, $2}' vacf.100ps.25fs.dat |fft.vacf 4001 25 $width
mv vdos.dat vdos.100ps.25fs.dat

awk '{print $1, $3}' vacf.100ps.25fs.dat |fft.vacf 4001 25 $width
mv vdos.dat vdos.x.100ps.25fs.dat

awk '{print $1, $4}' vacf.100ps.25fs.dat |fft.vacf 4001 25 $width
mv vdos.dat vdos.y.100ps.25fs.dat

awk '{print $1, $5}' vacf.100ps.25fs.dat |fft.vacf 4001 25 $width
mv vdos.dat vdos.z.100ps.25fs.dat

xmgrace vdos.*.dat &

