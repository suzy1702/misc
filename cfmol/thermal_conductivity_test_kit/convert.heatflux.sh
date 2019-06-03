#!/bin/bash

source /home/sane3962/.bashrc

j=0
while [ $j -lt 10 ]; do
  echo $j
  subdir="heat.$j"
  echo $subdir
  (cd "$subdir"; cp log.lammps heat1 && vi heat1; heatflux_format heat1 3 0.00025 0; mv HEATFLUX_formatted HEATFLUX && rm heat1;)
  (cd "$subdir"; lammps2dlpCONFIG sisnap.lammpstrj 0.00025 0 3 0 0 1 Si)
  let j=j+1
done
pwd

