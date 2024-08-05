#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cassert>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "cuda_defs.h"
#include "../defs.h"
#include "cuda_kernels.h"



// CUDA kernel to compute the forces for a subset of particles
__global__ void print_forces_subset(CUDA_REAL* result, int m) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m) {
		printf("acc: (%d) = %e\n", idx, result[_six*idx]);
	}
}



__global__	void initialize(CUDA_REAL* result, CUDA_REAL* diff, CUDA_REAL *magnitudes, int n, int m, int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m * n) {
		diff[_six*idx    ] = 0.;
		diff[_six*idx + 1] = 0.;
		diff[_six*idx + 2] = 0.;
		diff[_six*idx + 3] = 0.;
		diff[_six*idx + 4] = 0.;
		diff[_six*idx + 5] = 0.;

		magnitudes[_two*idx    ] = 0.;
		magnitudes[_two*idx + 1] = 0.;

		int i = idx / n;
		int j = idx % n;

		if (j == 0) {
			result[_six*i] = 0.;
			result[_six*i + 1] = 0.;
			result[_six*i + 2] = 0.;
			result[_six*i + 3] = 0.;
			result[_six*i + 4] = 0.;
			result[_six*i + 5] = 0.;
		}
	}
}

// CUDA kernel to compute pairwise differences for a subset of particles
__global__ void compute_pairwise_diff_subset(const CUDA_REAL* ptcl, CUDA_REAL* diff, int n, int m, const int* subset, int start) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m * n) {
		int i = subset[idx / n + start];
		int j = idx % n;
		idx *= _six;
		i *= _seven;
		j *= _seven;

		diff[idx]   = ptcl[j]   - ptcl[i];
		diff[idx+1] = ptcl[j+1] - ptcl[i+1];
		diff[idx+2] = ptcl[j+2] - ptcl[i+2];
		diff[idx+3] = ptcl[j+3] - ptcl[i+3];
		diff[idx+4] = ptcl[j+4] - ptcl[i+4];
		diff[idx+5] = ptcl[j+5] - ptcl[i+5];

		//printf("(%d,%d) = %e, %e, %e\n", i/_seven, j/_seven,  ptcl[i], ptcl[j], diff[idx]);
	}
}


__global__ void compute_magnitudes_subset(const CUDA_REAL *r2, const CUDA_REAL* diff, CUDA_REAL* magnitudes, int n, int m, int* subset, bool* neighbor2, int start) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < n * m) {
		int i = subset[idx / n + start];
		int j = idx % n;
		int six_idx = _six*idx;
		int two_idx = _two*idx;


		magnitudes[two_idx]   += diff[(six_idx)]    *diff[(six_idx)];
		magnitudes[two_idx]   += diff[(six_idx) + 1]*diff[(six_idx) + 1];
		magnitudes[two_idx]   += diff[(six_idx) + 2]*diff[(six_idx) + 2];
		magnitudes[two_idx+1] += diff[(six_idx)]    *diff[(six_idx) + 3];
		magnitudes[two_idx+1] += diff[(six_idx) + 1]*diff[(six_idx) + 4];
		magnitudes[two_idx+1] += diff[(six_idx) + 2]*diff[(six_idx) + 5];

		//printf("(%d,%d) = %e, %e\n", i, j,  magnitudes[two_idx], r2[i]);

		if (magnitudes[two_idx] <= r2[i]) {
			//printf("(%d, %d): %e, %e\n",subset[i], j, magnitudes[two_idx], r2[i]);
			//magnitudes[two_idx]   = -magnitudes[two_idx];
			magnitudes[two_idx]   = 0;
			neighbor2[idx] = true;
		}
		else {
			neighbor2[idx] = false;
		}
	}
}





// CUDA kernel to compute the forces for a subset of particles
__global__ void compute_forces_subset(const CUDA_REAL* ptcl, CUDA_REAL *diff, const CUDA_REAL* magnitudes, int n, int m, const int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m * n) {
		//int i = subset[idx / n];
		int i = idx / n;
		int j = idx % n;
		int six_idx = idx*_six;
		idx *= _two;
		CUDA_REAL acc[Dim], adot[Dim];

		if (magnitudes[idx] <= 0.) {
			acc[0]  = 0.;
			acc[1]  = 0.;
			acc[2]  = 0.;
			adot[0] = 0.;
			adot[1] = 0.;
			adot[2] = 0.;
		}
		else {
			CUDA_REAL scale = ptcl[_seven*j+6] / (magnitudes[idx] * sqrtf(magnitudes[idx]));
			acc[0]  = scale * diff[six_idx];
			acc[1]  = scale * diff[six_idx + 1];
			acc[2]  = scale * diff[six_idx + 2];

			adot[0] = scale * (diff[six_idx + 3] - 3*magnitudes[idx+1]*diff[six_idx    ]/magnitudes[idx]);
			adot[1] = scale * (diff[six_idx + 4] - 3*magnitudes[idx+1]*diff[six_idx + 1]/magnitudes[idx]);
			adot[2] = scale * (diff[six_idx + 5] - 3*magnitudes[idx+1]*diff[six_idx + 2]/magnitudes[idx]);
		}

		diff[six_idx]   = acc[0];
		diff[six_idx+1] = acc[1];
		diff[six_idx+2] = acc[2];
		diff[six_idx+3] = adot[0];
		diff[six_idx+4] = adot[1];
		diff[six_idx+5] = adot[2];

		//printf("compute_forces: (%d, %d) = %e\n", i, j,  diff[six_idx]);
	}
}


void reduce_forces_cublas(cublasHandle_t handle, const CUDA_REAL *diff, CUDA_REAL *result, int n, int m) {

	CUDA_REAL *d_matrix;
    cudaMalloc(&d_matrix, m * n * sizeof(CUDA_REAL));

    // Create a vector of ones for the summation
    double *ones;
    cudaMalloc(&ones, n * sizeof(double));
    double *h_ones = new double[n];
    for (int i = 0; i < n; ++i) {
        h_ones[i] = 1.0;
    }
    cudaMemcpy(ones, h_ones, n * sizeof(double), cudaMemcpyHostToDevice);
    // Initialize result array to zero
    cudaMemset(result, 0, m * 6 * sizeof(double));

    const double alpha = 1.0;
    const double beta = 0.0;

    // Sum over the second axis (n) for each of the 6 elements
    for (int i = 0; i < _six; ++i) {

		cublasDcopy(handle, m * n, diff + i, _six, d_matrix, 1);
        cublasDgemv(
            handle,
            CUBLAS_OP_T,  // Transpose
            n,            // Number of rows of the matrix A
            m,            // Number of columns of the matrix A
            &alpha,       // Scalar alpha
            d_matrix, // Pointer to the first element of the i-th sub-matrix
            n,     // Leading dimension of the sub-matrix
            ones,         // Pointer to the vector x
            1,            // Increment between elements of x
            &beta,        // Scalar beta
            result + i, // Pointer to the first element of the result vector
            _six             // Increment between elements of the result vector
        );
    }
    // Cleanup
    delete[] h_ones;
    cudaFree(ones);
	cudaFree(d_matrix);
}


/*
__device__ CUDA_REAL warpReduce(CUDA_REAL val) {
	val += __shfl_down_sync(0xffffffff, val, 16);
	val += __shfl_down_sync(0xffffffff, val, 8);
	val += __shfl_down_sync(0xffffffff, val, 4);
	val += __shfl_down_sync(0xffffffff, val, 2);
	val += __shfl_down_sync(0xffffffff, val, 1);
	return val;
}
*/

__inline__ __device__ CUDA_REAL warpReduce(CUDA_REAL val)
{
	for (int offset = warpSize/2; offset > 0; offset /= 2) 
		val += __shfl_down_sync(0xffffffff, val, offset);
	return val;
}


#define NEW_FORCE
#ifdef NEW_FORCE
__global__ void reduce_forces(const CUDA_REAL *diff, CUDA_REAL *result, int n, int m) {

	__shared__ CUDA_REAL warpSum[64]; // Assumes max 32 warps per block
	__shared__ CUDA_REAL res[_six]; //  this is for storing the results
	int lane = threadIdx.x % warpSize;
	int wid = threadIdx.x / warpSize;
	int bdim = blockDim.x;
	CUDA_REAL sum;
	int six_idx;
	int k,l, i = blockIdx.x, j;
	int a = (n+bdim-1)/(bdim);

	if (threadIdx.x < _six) 
		res[threadIdx.x] = 0.;
	__syncthreads();


	for (l=0;l<a*bdim;l+=bdim) {
		j = threadIdx.x + l;
		//printf("(%d,%d,%d)\n", i, j, l);
		six_idx = _six*(i*n+j);
		warpSum[wid] = 0.;
		__syncthreads();
#pragma unroll 
		for (k=0;k<_six;k++) { // ax ay az adotx adoty adotz

			sum = (i < m && j < n) ? diff[six_idx+k] : 0;

			/*
			if (k == 0)
				if (i < m && j < n)
					printf("(%d,%d) = %e, %e\n", i, j, sum, diff[six_idx+k]);
					*/

			// Warp reduce
			sum = warpReduce(sum);
			/*
			if (k == 0)
				if (i < m && j < n)
					printf("first reduction (%d,%d) = %e\n", i, j, sum);
					*/

			// Block reduce
			if (lane == 0) warpSum[wid] = sum;
			__syncthreads();

			if (wid == 0)
			{
				sum = (threadIdx.x < blockDim.x / warpSize) ? warpSum[lane] : 0;
				/*
				if (k == 0)
					if (i < m && j < n)
						printf("before second reduction (%d,%d) = %e\n", i, j, sum);
						*/

				sum = warpReduce(sum);

				/*
				if (k == 0)
					if (i < m && j < n)
						printf("second reduction (%d,%d) = %e\n", i, j, sum);
						*/

				if (lane == 0 && i < m) {
					res[k] += sum;
					//printf("%d = (%e,%e)\n", i, sum, res[k]);
				}
			}
			__syncthreads();
		} // reduce across threads
	}
	if (wid == 0 && lane == 0 && i < m) {
#pragma unroll
		for (k=0; k<_six;k++) {
			//printf("%d = (%e)\n", threadIdx.x, res[k]);
			result[_six*i+k] = res[k];
		}
	}
	__syncthreads();
}

#else
__global__ void reduce_forces(const CUDA_REAL *diff, CUDA_REAL *result, int n, int m) {
	int idx = blockIdx.x * n + threadIdx.x;
	__shared__ CUDA_REAL warpSum[64]; // Assumes max 32 warps per block
	int lane = threadIdx.x % warpSize;
	int wid = threadIdx.x / warpSize;
	CUDA_REAL sum;
	int i = blockIdx.x;
	int j = threadIdx.x;
	int six_idx = _six*(i*n+j);
	int k;

	//printf("old version\n");
#pragma unroll 
	for (k=0;k<_six;k++) {
		sum = (i < m && j < n) ? diff[six_idx+k] : 0;
		/*
		if (k == 0)
			if (i < m && j < n)
				printf("(%d,%d) = %e\n", blockIdx.x,threadIdx.x, diff[six_idx+k]);
				*/

		// Warp reduce
		sum = warpReduce(sum);

		// Block reduce
		if (lane == 0) warpSum[wid] = sum;
		__syncthreads();

		if (wid == 0)
		{
			sum = (threadIdx.x < blockDim.x / warpSize) ? warpSum[lane] : 0;
			sum = warpReduce(sum);
			if (lane == 0) result[_six*i+k] = sum;
		}
	}
}
#endif




#define MAX_SIZE 9 // maximum size of int array  blockDim.x*MaxSize = total size of int array 
#define NEW_V2 
#ifdef NEW_V2 // this works fine as :)
__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const CUDA_REAL* r2, const CUDA_REAL* magnitudes, int n, int m, const int *subset) {

	//int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	int bid = blockIdx.x;

		//printf("%d's bdim=%d, n=%d, blockdim=%d\n", tid,bid,n,blockDim.x);
	if (tid < n && bid < m) {
		extern __shared__ int sdata[];
		__shared__ int offset;

		sdata[tid+1]=0;
		if (tid==0)  {
			sdata[tid]=0;
			offset = 0;
		}
		__syncthreads();

		int gdim = min(m,gridDim.x);
		int bdim = min(n,blockDim.x);
		int list[MAX_SIZE];
		int n_num;
		int a = (n+bdim*MAX_SIZE-1)/(bdim*MAX_SIZE);
		int start, end;
		int i = subset[bid]; //target particle id
		int j, k;
		int idx = 0;


		//printf("assign_neighbor: %d\n",l);

		// background particles I
		for (k=0; k<a; k++) {

			start = tid+(k*MAX_SIZE*bdim);
			end   = min((k+1)*MAX_SIZE*bdim, n);
			sdata[tid+1] = 0;
			n_num = 0;

			//printf("tid=%d: (start, end) =(%d, %d)\n", tid, start, end);

			// background particles II
			for (j=start; j<end; j+=bdim) {
				if (i != j) {
					//printf("(l,j)=(%d,%d)\n",l,j);
					idx = _two*(n*bid+j);
					//printf("(%d, %d,%d) = %d, %e, %e\n", l, i, j, num_neighbor[l], magnitudes[idx], r2[l]);
					if (magnitudes[idx] < 0) {
						list[n_num] = j;
						n_num++;
						//printf("(%d,%d,%d) = %d, %e, %e\n", l, i, j, n_num, magnitudes[idx], r2[l]);
					}
				}
			} // endfor bk ptcl II
			sdata[tid+1] = n_num;

			__syncthreads();

			if (tid == 0) {

				for (j=2; j<=bdim; j++)
					sdata[j] += sdata[j-1];

				if ((offset+sdata[bdim]) > NumNeighborMax) {
					printf("blockid=%d, Too many neighbors (%d, %d)\n", bid, offset, sdata[bdim]);
					assert(offset+sdata[bdim] < NumNeighborMax);
				}

				/*
					 if (l == 11) {
					 printf("\n(%d, %d) = sdata[bdim]=%d\n", l, i, sdata[bdim]);
				//for (j=0; j<=bdim; j++) 
				//printf("%d, ",sdata[j]);
				//printf("\n");
				}
				 */
			}
			__syncthreads();

			/*
				 if (l==0)
				 printf("(%d,%d), (num, sdata) =%d, %d\n", l, tid, n_num, sdata[tid]);
			 */

			for (j=0;j<n_num;j++) {
				neighbor[NumNeighborMax*bid+offset+sdata[tid]+j] = list[j];
				//printf("(%d,%d), j=%d\n", l, tid, list[j]);
			}
			__syncthreads();

			if (tid == 0) {
				//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
				offset += sdata[bdim];
				//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
			}
			__syncthreads();
		} //endfor bk ptcl I

		if (tid == 0) {
			num_neighbor[bid] = offset; // bid shoud be modified
			offset = 0;
			sdata[0] = 0;
		}
	} // m*n stuff
}

#elif defined(NEW_V1) // works well

__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const CUDA_REAL* r2, const CUDA_REAL* magnitudes, int n, int m, const int *subset) {

	//int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	int bid = blockIdx.x;

		//printf("%d's bdim=%d, n=%d, blockdim=%d\n", tid,bid,n,blockDim.x);
	if (tid < n && bid < m) {
		extern __shared__ int sdata[];
		__shared__ int offset;

		sdata[tid+1]=0;
		if (tid==0)  {
			sdata[tid]=0;
			offset = 0;
		}
		__syncthreads();

		int gdim = min(m,gridDim.x);
		int bdim = min(n,blockDim.x);
		int list[MAX_SIZE];
		int n_num;
		int a = (n+bdim*MAX_SIZE-1)/(bdim*MAX_SIZE);
		int start, end;
		int i; //target particle id
		int j, k, l;
		int idx = 0;


		//printf("%d's bdim=%d, n=%d, blockdim=%d\n", tid,bdim,n,blockDim.x);
		// target particles
		for (l=bid; l<m; l+=gdim) {
			i = subset[l];
			//printf("assign_neighbor: %d\n",l);

			// background particles I
			for (k=0; k<a; k++) {

				start = tid+(k*MAX_SIZE*bdim);
				end   = min((k+1)*MAX_SIZE*bdim, n);
				sdata[tid+1] = 0;
				n_num = 0;

				//printf("tid=%d: (start, end) =(%d, %d)\n", tid, start, end);

				// background particles II
				for (j=start; j<end; j+=bdim) {
					if (i != j) {
						//printf("(l,j)=(%d,%d)\n",l,j);
						idx = _two*(n*l+j);
						//printf("(%d, %d,%d) = %d, %e, %e\n", l, i, j, num_neighbor[l], magnitudes[idx], r2[l]);
						if (magnitudes[idx] < 0) {
							list[n_num] = j;
							n_num++;
							//printf("(%d,%d,%d) = %d, %e, %e\n", l, i, j, n_num, magnitudes[idx], r2[l]);
						}
					}
				} // endfor bk ptcl II
				sdata[tid+1] = n_num;

				__syncthreads();

				if (tid == 0) {

					for (j=2; j<=bdim; j++)
						sdata[j] += sdata[j-1];

					if ((offset+sdata[bdim]) > NumNeighborMax) {
						printf("blockid=%d, Too many neighbors (%d, %d)\n", bid, offset, sdata[bdim]);
						assert(offset+sdata[bdim] < NumNeighborMax);
					}

					/*
					if (l == 11) {
						printf("\n(%d, %d) = sdata[bdim]=%d\n", l, i, sdata[bdim]);
						//for (j=0; j<=bdim; j++) 
							//printf("%d, ",sdata[j]);
						//printf("\n");
					}
					*/
				}
				__syncthreads();

				/*
					 if (l==0)
					 printf("(%d,%d), (num, sdata) =%d, %d\n", l, tid, n_num, sdata[tid]);
				 */

				for (j=0;j<n_num;j++) {
					neighbor[NumNeighborMax*l+offset+sdata[tid]+j] = list[j];
					//printf("(%d,%d), j=%d\n", l, tid, list[j]);
				}
				__syncthreads();

				if (tid == 0) {
					//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
					offset += sdata[bdim];
					//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
				}
				__syncthreads();
			} //endfor bk ptcl I
			if (tid == 0) {
				num_neighbor[l] = offset; // bid shoud be modified
				offset = 0;
				sdata[0] = 0;
			}
			__syncthreads();
		} // endfor target paticles
	} // m*n stuff
}





#else

__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const REAL* r2, const REAL* magnitudes, int n, int m, const int *subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m) {
		int i = subset[idx];
		int k = 0;

		for (int j = 0; j < n; j++) {
			if (i != j) {
				k = _two*(n*idx+j);
				if (magnitudes[k] < 0) {
					//printf("(%d, %d,%d) = %d, %e, %e\n", idx, i, j, num_neighbor[idx], magnitudes[k], r2[i]);
					neighbor[NumNeighborMax*idx+num_neighbor[idx]] = j;
					num_neighbor[idx]++;
					if (num_neighbor[idx] > 100)  {
						//printf("Error: (%d, %d,%d) = %d, %e, %e\n", idx, i, j, num_neighbor[idx], magnitudes[k], r2[i]);
						assert(num_neighbor[idx] < 100);
						return;
					}
				}
			}
		}
	}
}

#endif
