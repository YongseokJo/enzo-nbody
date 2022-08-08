#!/bin/sh  

# This makes enzo.exe on Happiness@SNU

echo "Making enzo on"
pwd
make clean
make default
cd ../../
./configure
cd src/enzo/

module add icc
module add mpi

make machine-linux-mpich
make precision-64 integers-32 particle-id-32 max-baryons-30 opt-aggressive lcaperf-no max-tasks-per-node-36 grackle-yes new-problem-types-yes
make show-config
make show-flags
make -j3
echo "Make done!"
pwd
date
