#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,"make[1]: Warning: File `DEPEND' has modification time 123 s in the future\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPP = /bin/cpp\n");
   fprintf (fp,"CC  = /appl/intel/oneapi/mpi/2021.4.0/bin/mpicc\n");
   fprintf (fp,"CXX = /appl/intel/oneapi/mpi/2021.4.0/bin/mpicxx\n");
   fprintf (fp,"FC  = /appl/intel/oneapi/mpi/2021.4.0/bin/mpifc\n");
   fprintf (fp,"F90 = /appl/intel/oneapi/mpi/2021.4.0/bin/mpif90\n");
   fprintf (fp,"LD  = /appl/intel/oneapi/mpi/2021.4.0/bin/mpicxx\n");
   fprintf (fp,"\n");
   fprintf (fp,"DEFINES = -DLINUX -DH5_USE_16_API   -D__max_subgrids=100000 -D__max_baryons=30 -D__max_cpu_per_node=36 -D__memory_pool_size=100000 -DINITS64 -DSMALL_INTS -DCONFIG_PINT_4 -DIO_32   -DNEW_PROBLEM_TYPES -DUSE_MPI   -DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8  -DUSE_HDF5_GROUPS      -DNEW_GRID_IO -DFAST_SIB   -DNBODY -D GPU     -DENZO_PERFORMANCE  -DUSE_GRACKLE  -DUSE_UUID -DSAB\n");
   fprintf (fp,"\n");
   fprintf (fp,"INCLUDES = -I/home/sykim/local/hdf5-1.12.1/include  -I/appl/intel/oneapi/mpi/latest/include       -I/home/sykim/local/grackle/include    -I.\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPPFLAGS = -P -traditional \n");
   fprintf (fp,"CFLAGS   =  -O3 -g\n");
   fprintf (fp,"CXXFLAGS =  -O3 -g\n");
   fprintf (fp,"FFLAGS   = -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -O3 -g\n");
   fprintf (fp,"F90FLAGS =  -O3 -g\n");
   fprintf (fp,"LDFLAGS  =  -O3 -g\n");
   fprintf (fp,"\n");
   fprintf (fp,"LIBS     = -L/home/sykim/local/hdf5-1.12.1/lib -lhdf5 -lz  -lgfortran   -L/appl/intel/oneapi/mpi/latest/lib -lmpi        -L/usr/local/cuda-11.8/lib64 -lcudart -L/home/sykim/local/grackle/lib -lgrackle\n");
   fprintf (fp,"\n");
   fprintf (fp,"make[1]: warning:  Clock skew detected.  Your build may be incomplete.\n");
}
