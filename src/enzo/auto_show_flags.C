#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"CPP = /mnt/sw/nix/store/494gvfa9gj2ibqg8210v89h2iljfgqj8-gcc-7.5.0/bin/cpp\n");
   fprintf (fp,"CC  = /mnt/sw/nix/store/gxy45l70ghs65xpwq73ayg7wb0a6sp2s-openmpi-4.0.7/bin/mpicc\n");
   fprintf (fp,"CXX = /mnt/sw/nix/store/gxy45l70ghs65xpwq73ayg7wb0a6sp2s-openmpi-4.0.7/bin/mpicxx\n");
   fprintf (fp,"FC  = /mnt/sw/nix/store/gxy45l70ghs65xpwq73ayg7wb0a6sp2s-openmpi-4.0.7/bin/mpif77\n");
   fprintf (fp,"F90 = /mnt/sw/nix/store/gxy45l70ghs65xpwq73ayg7wb0a6sp2s-openmpi-4.0.7/bin/mpif90\n");
   fprintf (fp,"LD  = /mnt/sw/nix/store/gxy45l70ghs65xpwq73ayg7wb0a6sp2s-openmpi-4.0.7/bin/mpicxx\n");
   fprintf (fp,"\n");
   fprintf (fp,"DEFINES = -DLINUX -DH5_USE_16_API   -D__max_subgrids=100000 -D__max_baryons=30 -D__max_cpu_per_node=36 -D__memory_pool_size=100000 -DINITS64 -DSMALL_INTS -DCONFIG_PINT_4 -DIO_32   -DNEW_PROBLEM_TYPES -DUSE_MPI   -DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8  -DUSE_HDF5_GROUPS      -DNEW_GRID_IO -DFAST_SIB   -DNBODY -D GPU     -DENZO_PERFORMANCE  -DUSE_GRACKLE  -DUSE_UUID -DSAB\n");
   fprintf (fp,"\n");
   fprintf (fp,"INCLUDES = -I/mnt/sw/nix/store/5lm2q7hm14wv3dnc09144zi5bh7v7p42-hdf5-1.8.22/include  -I/include       -I/mnt/home/yjo10/packages/grackle/build/include    -I.\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPPFLAGS = -P -traditional\n");
   fprintf (fp,"CFLAGS   =  -O3 -g\n");
   fprintf (fp,"CXXFLAGS = -lcudart -O3 -g\n");
   fprintf (fp,"FFLAGS   = -lmpi -lrt -ldl -Wall  -O3 -g\n");
   fprintf (fp,"F90FLAGS =  -O3 -g\n");
   fprintf (fp,"LDFLAGS  =  -O3 -g\n");
   fprintf (fp,"\n");
   fprintf (fp,"LIBS     = -L/mnt/sw/nix/store/5lm2q7hm14wv3dnc09144zi5bh7v7p42-hdf5-1.8.22/lib -lhdf5 -lz  -lgfortran   -L/lib -lmpi        -L/mnt/sw/nix/store/bdhdh478f6slibd9zpgmgw8grnqq78im-cuda-11.4.4/lib64 -lcudart -L/mnt/home/yjo10/packages/grackle/build/lib -lgrackle\n");
   fprintf (fp,"\n");
}
