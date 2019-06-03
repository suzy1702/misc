#!/bin/sh

source /home/sane3962/.bashrc

#a=0
#for subdir in */; do
#  (cd "$subdir"; ./jobsubmit.sh $a)
#  let a=a+1
#  echo $a
#  cd ..
#done

j=0
while [ $j -lt 10 ]; do
  subdir="heat.$j"
  (cd "$subdir"; sbatch lammps.run.gk.bulk.Si)
  echo $j
  echo $subdir
  let j=j+1
done
pwd
#a=0
#find -name jobsubmit.sh -mindepth 2 -maxdepth 2 -exec sh {} \;
#find -name jobsubmit.sh -mindepth 2 -maxdepth 2 -exec sh {}\;
#echo $a
#a=a+1
#echo "${PWD##*/}"
#exit 0
