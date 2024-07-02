#!/bin/sh  

# This makes enzo.exe on Happiness@SNU

echo "Making enzo on"
pwd
make clean
make default
cd ../../
./configure
cd src/enzo/

module add icc/latest
module add mpi/latest

make machine-linux-mpich
make precision-64 integers-32 particle-id-32 max-baryons-30 opt-aggressive lcaperf-no max-tasks-per-node-36 grackle-yes new-problem-types-yes nbody-yes photon-no cuda-no
make show-config
make show-flags
make -j8
#cp enzo.exe enzo_spare.exe
#cp enzo.exe /data1/wispedia/enzo_nfw/lessSFB
echo "Make done!"
pwd
date
