#!/bin/bash

source /home/sane3962/.bashrc

for((j=1;j<10;j++))
do
  mkdir heat.${j}
  cp heat.0/in.kappa.gk.bulk.Si heat.${j}/
  cp heat.0/lammps.run.gk.bulk.Si heat.${j}/
done

