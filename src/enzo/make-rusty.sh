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
#module add intel-oneapi-compilers
#module add intel-oneapi-mpi
#module add modules/2.1-20220630
#module add modules/2.2-20230808 /2.1.1-20230405
#2.1.1-20230405
#module add gcc/10.3.0
#module add openmpi

module add modules/2.0-20220630 
module add gcc/7.5.0
#module add openmpi/1.10.7
module add openmpi/cuda-4.0.7
module add cuda/11.4.4
module add hdf5/1.8.22
#module add ucx
#module add cuda/12.1.1
#module add openmpi/1.10.7
#module add intel-parallel-studio

make machine-rusty
make precision-64 integers-32 particle-id-32 max-baryons-30 lcaperf-no max-tasks-per-node-36 grackle-yes new-problem-types-yes photon-no nbody-yes opt-aggressive cuda-no
make show-config
make show-flags
##make -j3
make -j16
#cp enzo.exe enzo_spare.exe
#cp enzo.exe enzo_spare_v100.exe
#cp enzo.exe enzo_test_openmpi_escape.exe
#cp enzo.exe enzo_test_openmpi.exe
#cp enzo.exe enzo_sf_test.exe
cp enzo.exe enzo_debug.exe
#cp enzo.exe enzo_cosmo.exe
#cp enzo.exe enzo_orbit.exe
#cp enzo.exe enzo_nbn.exe
#cp enzo.exe enzo_test_irr.exe
#cp enzo.exe enzo_const.exe
#cp enzo.exe enzo_1e6.exe
echo "Make done!"
pwd
date
