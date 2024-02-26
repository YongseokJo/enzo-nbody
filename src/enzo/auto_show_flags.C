#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"CPP = /mnt/sw/nix/store/6560mkp838syd8jpp6gdyyisczwcvs67-gcc-11.4.0/bin/cpp\n");
   fprintf (fp,"CC  = /mnt/sw/nix/store/5gzzqpkmbz1rjmlwxa3f9mvmzalkzqlk-intel-oneapi-mpi-2021.9.0/mpi/2021.9.0/bin/mpicc\n");
   fprintf (fp,"CXX = /mnt/sw/nix/store/5gzzqpkmbz1rjmlwxa3f9mvmzalkzqlk-intel-oneapi-mpi-2021.9.0/mpi/2021.9.0/bin/mpicxx\n");
   fprintf (fp,"FC  = /mnt/sw/nix/store/5gzzqpkmbz1rjmlwxa3f9mvmzalkzqlk-intel-oneapi-mpi-2021.9.0/mpi/2021.9.0/bin/mpif77\n");
   fprintf (fp,"F90 = /mnt/sw/nix/store/5gzzqpkmbz1rjmlwxa3f9mvmzalkzqlk-intel-oneapi-mpi-2021.9.0/mpi/2021.9.0/bin/mpif90\n");
   fprintf (fp,"LD  = /mnt/sw/nix/store/5gzzqpkmbz1rjmlwxa3f9mvmzalkzqlk-intel-oneapi-mpi-2021.9.0/mpi/2021.9.0/bin/mpicxx\n");
   fprintf (fp,"\n");
   fprintf (fp,"DEFINES = -DLINUX -DH5_USE_16_API   -D__max_subgrids=100000 -D__max_baryons=30 -D__max_cpu_per_node=36 -D__memory_pool_size=100000 -DINITS64 -DSMALL_INTS -DCONFIG_PINT_4 -DIO_32   -DNEW_PROBLEM_TYPES -DUSE_MPI   -DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8  -DUSE_HDF5_GROUPS      -DNEW_GRID_IO -DFAST_SIB   -DNBODY -D GPU     -DENZO_PERFORMANCE  -DUSE_GRACKLE  -DUSE_UUID -DSAB\n");
   fprintf (fp,"\n");
   fprintf (fp,"INCLUDES = -I/mnt/sw/nix/store/bjxnws41q81ccp0lybgzacy5wdv2n82g-hdf5-1.8.22/include  -I/mnt/sw/nix/store/5gzzqpkmbz1rjmlwxa3f9mvmzalkzqlk-intel-oneapi-mpi-2021.9.0/include       -I/mnt/home/yjo10/packages/grackle/build/include    -I.\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPPFLAGS = -P -traditional\n");
   fprintf (fp,"CFLAGS   =  -O3 -g\n");
   fprintf (fp,"CXXFLAGS =  -O3 -g\n");
   fprintf (fp,"FFLAGS   = -lmpi -lrt -ldl -W -O3 -g\n");
   fprintf (fp,"F90FLAGS =  -O3 -g\n");
   fprintf (fp,"LDFLAGS  =  -O3 -g\n");
   fprintf (fp,"\n");
   fprintf (fp,"LIBS     = -L/mnt/sw/nix/store/bjxnws41q81ccp0lybgzacy5wdv2n82g-hdf5-1.8.22/lib -lhdf5 -lz  -lgfortran   -L/mnt/sw/nix/store/5gzzqpkmbz1rjmlwxa3f9mvmzalkzqlk-intel-oneapi-mpi-2021.9.0/lib -lmpi        -L/mnt/sw/nix/store/zi2wc26znf75csf5hhz77p0d2bbz53ih-cuda-11.8.0/lib64 -lcudart -L/mnt/home/yjo10/packages/grackle/build/lib -lgrackle\n");
   fprintf (fp,"\n");
}
