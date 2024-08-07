#include <stdio.h>
void auto_show_config(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"   MACHINE: Generic Ubuntu 8.04\n");
   fprintf (fp,"   MACHINE-NAME: rusty\n");
   fprintf (fp,"\n");
   fprintf (fp,"   PARAMETER_MAX_SUBGRIDS  [max-subgrids-###]                : 100000\n");
   fprintf (fp,"   PARAMETER_MAX_BARYONS  [max-baryons-###]                  : 30\n");
   fprintf (fp,"   PARAMETER_MAX_TASKS_PER_NODE  [max-tasks-per-node-###]    : 36\n");
   fprintf (fp,"   PARAMETER_MEMORY_POOL_SIZE  [memory-pool-###]             : 100000\n");
   fprintf (fp,"\n");
   fprintf (fp,"   CONFIG_PRECISION  [precision-{32,64}]                     : 64\n");
   fprintf (fp,"   CONFIG_PARTICLES  [particles-{32,64,128}]                 : 64\n");
   fprintf (fp,"   CONFIG_INTEGERS  [integers-{32,64}]                       : 32\n");
   fprintf (fp,"   CONFIG_PARTICLE_IDS  [particle-id-{32,64}]                : 32\n");
   fprintf (fp,"   CONFIG_INITS  [inits-{32,64}]                             : 64\n");
   fprintf (fp,"   CONFIG_IO  [io-{32,64}]                                   : 32\n");
   fprintf (fp,"   CONFIG_USE_MPI  [use-mpi-{yes,no}]                        : yes\n");
   fprintf (fp,"   CONFIG_TASKMAP  [taskmap-{yes,no}]                        : no\n");
   fprintf (fp,"   CONFIG_PACKED_AMR  [packed-amr-{yes,no}]                  : yes\n");
   fprintf (fp,"   CONFIG_PACKED_MEM  [packed-mem-{yes,no}]                  : no\n");
   fprintf (fp,"   CONFIG_LCAPERF  [lcaperf-{yes,no}]                        : no\n");
   fprintf (fp,"   CONFIG_PAPI  [papi-{yes,no}]                              : no\n");
   fprintf (fp,"   CONFIG_PYTHON  [python-{yes,no}]                          : no\n");
   fprintf (fp,"   CONFIG_NEW_PROBLEM_TYPES  [new-problem-types-{yes,no}]    : yes\n");
   fprintf (fp,"   CONFIG_ECUDA  [cuda-{yes,no}]                             : no\n");
   fprintf (fp,"   CONFIG_NBODY  [nbody-{yes,no}]                            : yes\n");
   fprintf (fp,"   CONFIG_OOC_BOUNDARY  [ooc-boundary-{yes,no}]              : no\n");
   fprintf (fp,"   CONFIG_ACCELERATION_BOUNDARY  [acceleration-boundary-{yes,no}] : yes\n");
   fprintf (fp,"   CONFIG_OPT  [opt-{warn,debug,cudadebug,high,aggressive}]  : aggressive\n");
   fprintf (fp,"   CONFIG_TESTING  [testing-{yes,no}]                        : no\n");
   fprintf (fp,"   CONFIG_PHOTON  [photon-{yes,no}]                          : no\n");
   fprintf (fp,"   CONFIG_HYPRE  [hypre-{yes,no}]                            : no\n");
   fprintf (fp,"   CONFIG_EMISSIVITY  [emissivity-{yes,no}]                  : no\n");
   fprintf (fp,"   CONFIG_USE_HDF4  [use-hdf4-{yes,no}]                      : no\n");
   fprintf (fp,"   CONFIG_NEW_GRID_IO  [newgridio-{yes,no}]                  : yes\n");
   fprintf (fp,"   CONFIG_BITWISE_IDENTICALITY  [bitwise-{yes,no}]           : no\n");
   fprintf (fp,"   CONFIG_FAST_SIB  [fastsib-{yes,no}]                       : yes\n");
   fprintf (fp,"   CONFIG_GRAVITY_4S  [gravity-4s-{yes,no}]                  : no\n");
   fprintf (fp,"   CONFIG_ENZO_PERFORMANCE  [enzo-performance-{yes,no}]      : yes\n");
   fprintf (fp,"   CONFIG_GRACKLE  [grackle-{yes,no}]                        : yes\n");
   fprintf (fp,"   CONFIG_LOG2ALLOC  [log2alloc-{yes,no}]                    : no\n");
   fprintf (fp,"   CONFIG_UUID  [uuid-{yes,no}]                              : yes\n");
   fprintf (fp,"\n");
}
