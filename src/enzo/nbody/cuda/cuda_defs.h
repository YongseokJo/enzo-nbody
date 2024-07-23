#define _PROFILE
#define _six 6
#define _two 2
#define _seven 7
#define new_size(A) ((A > 1024) ? int(pow(2,ceil(log(A)/log(2.0)))) : 1024)


typedef double CUDA_REAL;
