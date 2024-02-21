#!/bin/sh  

# This makes enzo.exe on Rusty@Flatiron

echo "Making enzo on"
pwd
make clean
make default
cd ../../
./configure
cd src/enzo/

#module add openmpi4
#module add modules/2.1.1-20230405
#module add gcc/10.4.0
module add intel-oneapi-compilers
module add intel-oneapi-mpi
#module add ucx
module add cuda/12.1.1
#module add openmpi/1.10.7
#module add intel-parallel-studio
module add hdf5/1.8.22

make machine-rusty
make precision-64 integers-32 particle-id-32 max-baryons-30 lcaperf-no max-tasks-per-node-36 grackle-yes new-problem-types-yes photon-no nbody-yes opt-aggressive cuda-no
make show-config
make show-flags
make -j3
#cp enzo.exe enzo_spare.exe
#cp enzo.exe enzo_spare_v100.exe
cp enzo.exe enzo_test.exe
#cp enzo.exe enzo_const.exe
#cp enzo.exe enzo_1e6.exe
echo "Make done!"
pwd
date
