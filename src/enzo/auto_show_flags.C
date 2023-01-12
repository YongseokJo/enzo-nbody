#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"CPP = /bin/cpp\n");
   fprintf (fp,"CC  = /home/sykim/local/openmpi-4.1.3/bin/mpicc\n");
   fprintf (fp,"CXX = /home/sykim/local/openmpi-4.1.3/bin/mpic++\n");
   fprintf (fp,"FC  = /bin/gfortran\n");
   fprintf (fp,"F90 = /bin/gfortran\n");
   fprintf (fp,"LD  = /home/sykim/local/openmpi-4.1.3/bin/mpic++\n");
   fprintf (fp,"\n");
   fprintf (fp,"DEFINES = -DLINUX -DH5_USE_16_API   -D__max_subgrids=100000 -D__max_baryons=30 -D__max_cpu_per_node=36 -D__memory_pool_size=100000 -DINITS64 -DSMALL_INTS -DCONFIG_PINT_4 -DIO_32   -DNEW_PROBLEM_TYPES -DUSE_MPI   -DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8  -DUSE_HDF5_GROUPS      -DNEW_GRID_IO -DFAST_SIB   -DNBODY    -DENZO_PERFORMANCE  -DUSE_GRACKLE  -DUSE_UUID -DSAB\n");
   fprintf (fp,"\n");
   fprintf (fp,"INCLUDES = -I/home/sykim/local/hdf5-1.12.1/include  -I/home/sykim/local/openmpi-4.1.3/include       -I/home/sykim/local/grackle/include    -I.\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPPFLAGS = \n");
   fprintf (fp,"CFLAGS   =  -O3 -g\n");
   fprintf (fp,"CXXFLAGS =  -O3 -g\n");
   fprintf (fp,"FFLAGS   = -fno-second-underscore -ffixed-line-length-132 -O3 -g\n");
   fprintf (fp,"F90FLAGS = -fno-second-underscore -O3 -g\n");
   fprintf (fp,"LDFLAGS  =  -O3 -g\n");
   fprintf (fp,"\n");
   fprintf (fp,"LIBS     = -L/home/sykim/local/hdf5-1.12.1//lib -lhdf5  -lgfortran  -L/home/sykim/local/openmpi-4.1.3/lib -lmpi         -L/home/sykim/local/grackle/lib -lgrackle\n");
   fprintf (fp,"\n");
}
