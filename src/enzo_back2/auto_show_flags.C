#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"CPP = /usr/bin/cpp\n");
   fprintf (fp,"CC  = /appl/intel/oneapi/mpi/2021.4.0/bin/mpicc\n");
   fprintf (fp,"CXX = /appl/intel/oneapi/mpi/2021.4.0/bin/mpicxx\n");
   fprintf (fp,"FC  = /appl/intel/oneapi/mpi/2021.4.0/bin/mpifc\n");
   fprintf (fp,"F90 = /appl/intel/oneapi/mpi/2021.4.0/bin/mpif90\n");
   fprintf (fp,"LD  = /appl/intel/oneapi/mpi/2021.4.0/bin/mpicxx\n");
   fprintf (fp,"\n");
   fprintf (fp,"DEFINES = -DLINUX -DH5_USE_16_API   -D__max_subgrids=100000 -D__max_baryons=30 -D__max_cpu_per_node=36 -D__memory_pool_size=100000 -DINITS64 -DSMALL_INTS -DCONFIG_PINT_4 -DIO_32   -DNEW_PROBLEM_TYPES -DUSE_MPI   -DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8  -DUSE_HDF5_GROUPS   -DTRANSFER   -DNEW_GRID_IO -DFAST_SIB      -DENZO_PERFORMANCE  -DUSE_GRACKLE  -DUSE_UUID -DSAB\n");
   fprintf (fp,"\n");
   fprintf (fp,"INCLUDES = -I/home/kerex/packages/hdf5-1.8.20/build/include  -I/appl/intel/oneapi/mpi/2021.4.0/include       -I/home/kerex/packages/grackle/build/include    -I.\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPPFLAGS = -P -traditional \n");
   fprintf (fp,"CFLAGS   =  -O3 -g\n");
   fprintf (fp,"CXXFLAGS =  -O3 -g\n");
   fprintf (fp,"FFLAGS   =  -O3 -g\n");
   fprintf (fp,"F90FLAGS =  -O3 -g\n");
   fprintf (fp,"LDFLAGS  =  -O3 -g\n");
   fprintf (fp,"\n");
   fprintf (fp,"LIBS     = -L/home/kerex/packages/hdf5-1.8.20/build/lib -lhdf5 -lz  -lgfortran   -L/appl/intel/oneapi/mpi/2021.4.0/lib -lmpi        -L/home/kerex/packages/grackle/build/lib -lgrackle\n");
   fprintf (fp,"\n");
}
