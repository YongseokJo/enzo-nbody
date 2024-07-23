#ifndef KERNEL_H
#define KERNEL_H
#include "../defs.h"


__global__	void initialize(CUDA_REAL* result, int* neighbor, int* num_neighbor, CUDA_REAL* diff, CUDA_REAL *magnitudes, int n, int m, int* subset);
__global__ void compute_pairwise_diff_subset(const CUDA_REAL* ptcl, CUDA_REAL* diff, int n, int m, const int* subset);
__global__ void compute_magnitudes_subset(const CUDA_REAL *r2, const CUDA_REAL* diff, CUDA_REAL* magnitudes, int n, int m, int* subset);
__global__ void compute_forces_subset(const CUDA_REAL* ptcl, CUDA_REAL *diff, const CUDA_REAL* magnitudes, int n, int m, const int* subset);
__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const CUDA_REAL* r2, const CUDA_REAL* magnitudes, int n, int m, const int *subset);
__global__ void reduce_forces(const CUDA_REAL *diff, CUDA_REAL *result, int n, int m);

/*
__device__ void _addition(Result &result, const Result res);
__device__ void _copy(Result &result, const Result res);
*/

#endif
