#!/bin/sh  

# This makes enzo.exe on Happiness@SNU

echo "Making enzo on"
pwd
make clean
make default
cd ../../
./configure
cd src/enzo/

module add gcc/7.5.0
module add openmpi4
module add ucx
#module add openmpi/1.10.7
#module add intel-parallel-studio
module add hdf5/1.8.22

make machine-rusty
make precision-64 integers-32 particle-id-32 max-baryons-30 opt-aggressive lcaperf-no max-tasks-per-node-36 grackle-yes new-problem-types-yes photon-no nbody-yes
make show-config
make show-flags
make -j3
cp enzo.exe enzo_spare.exe
echo "Make done!"
pwd
date
