#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <cassert>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "../defs.h"
#include "cuda_global.h"
#include "cuda_defs.h"
#include "cuda_kernels.h"
#include "cuda_routines.h"


extern int _NNB;
int _NNB = 0;
static int nbodymax;
static int devid, numGPU;
static bool is_open = false;
static bool devinit = false;
static bool first   = true;
static int variable_size;


extern CUDA_REAL *h_ptcl, *d_ptcl; //, *background;
extern CUDA_REAL *h_result, *d_result;
extern CUDA_REAL *d_diff, *d_magnitudes, *d_r2;
extern int *h_neighbor, *d_neighbor, *h_num_neighbor, *d_num_neighbor;
extern int *d_target;

CUDA_REAL *h_ptcl=nullptr, *d_ptcl=nullptr;; //, *background;
CUDA_REAL *h_result=nullptr, *d_result=nullptr;
CUDA_REAL *d_diff=nullptr,*d_magnitudes=nullptr, *d_r2=nullptr;
int *h_neighbor=nullptr, *d_neighbor=nullptr, *d_num_neighbor=nullptr, *h_num_neighbor=nullptr;
int *d_target=nullptr;


extern cudaStream_t stream;
cudaStream_t stream;

extern CUDA_REAL *h_diff, *h_magnitudes;
CUDA_REAL *h_diff, *h_magnitudes;

/*************************************************************************
 *	 Computing Acceleration
 *************************************************************************/

void GetAcceleration(
		int NumTarget,
		int h_target_list[],
		CUDA_REAL acc[][3],
		CUDA_REAL adot[][3],
		int NumNeighbor[],
		int **NeighborList
		) {

	assert(is_open);
	assert((NumTarget > 0) && (NumTarget <= _NNB));

	int minGridSize, blockSize, gridSize;
	int sharedMemSize;

	//cudaStreamCreate(&stream);

	cublasHandle_t handle;
	initializeCudaAndCublas(&handle);

	/*
	for(int i=0; i<NumTarget; i++) {
		d_result[i].clear();
		d_neighbor[i].clear();
		d_dist = 0.;
	}
	*/
	/*
	fprintf(stderr,"\ntargets=");
	for(int i=0; i<NumTarget; i++) {
		fprintf(stderr,"%d, ", h_target_list[i]);
	}
	fprintf(stderr,"\n");
	*/


	//toDevice(h_target_list, d_target, NumTarget, stream);
	toDevice(h_target_list, d_target, NumTarget, stream);

	// Kernel launch parameters
	//dim3 blockSize(variable_size);
	//dim3 gridSize(NumTarget);
	//dim3 gridSize((NumTarget * NNB + blockSize.x - 1) / blockSize.x);

	// Compute pairwise differences for the subset

	//blockSize = variable_size;
	//gridSize = NumTarget;
	int total_data_num = new_size(_NNB*NumTarget);
	/******* Initialize *********/
	checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
			 	initialize, 0, 0));	
	gridSize = (total_data_num + blockSize - 1) / blockSize;

	initialize<<<gridSize, blockSize, 0, stream>>>\
		(d_result, d_neighbor, d_num_neighbor, d_diff, d_magnitudes, _NNB, NumTarget, d_target);
	cudaDeviceSynchronize();


	/******* Differencese *********/
	checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
			 	compute_pairwise_diff_subset, 0, 0));	
	gridSize = (total_data_num + blockSize - 1) / blockSize;

	compute_pairwise_diff_subset<<<gridSize, blockSize, 0, stream>>>\
		(d_ptcl, d_diff, _NNB, NumTarget, d_target);
	cudaDeviceSynchronize();

	/******* Magnitudes *********/
	checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
			 	compute_magnitudes_subset, 0, 0));	
	gridSize = (total_data_num + blockSize - 1) / blockSize;

	compute_magnitudes_subset<<<gridSize, blockSize, 0, stream>>>\
		(d_r2, d_diff, d_magnitudes, _NNB, NumTarget, d_target);
	cudaDeviceSynchronize();

	/******* Force *********/
	checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
			 	compute_forces_subset, 0, 0));
	gridSize = (total_data_num + blockSize - 1) / blockSize;

	compute_forces_subset<<<gridSize, blockSize, 0, stream>>>\
		(d_ptcl, d_diff, d_magnitudes, _NNB, NumTarget, d_target);

	/******* Neighborhood *********/
#define MAX_SIZE 9
	sharedMemSize = ((MAX_SIZE+1)*512) * sizeof(int);
	checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
				assign_neighbor, sharedMemSize, 0));
	gridSize = (total_data_num + blockSize - 1) / blockSize;

	
	//blockSize = std::min(blockSize, 512);
	//gridSize = (_NNB * NumTarget + blockSize - 1) / blockSize;

	//blockSize = variable_size;
	//gridSize = NumTarget;

	assign_neighbor<<<gridSize, blockSize, sharedMemSize, stream>>>\
		(d_neighbor, d_num_neighbor, d_r2, d_magnitudes, _NNB, NumTarget, d_target);
	cudaDeviceSynchronize();

	/******* Reduction *********/
	checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
			 	reduce_forces, 0, 0));
	gridSize = (total_data_num + blockSize - 1) / blockSize;
	//blockSize = _NNB;
	//gridSize  = NumTarget;
	//blockSize = 128;
	//blockSize = variable_size;
	//gridSize = NumTarget;


	//	sharedMemSize = 256 * sizeof(double);
	reduce_forces<<<gridSize, blockSize, 0, stream>>>\
		(d_diff, d_result, _NNB, NumTarget);
	cudaDeviceSynchronize();

	/*
	print_forces_subset<<<gridSize, blockSize>>>\
		(d_result, NumTarget);
		*/




	/*
	toHost(h_diff, d_diff, _six*NumTarget*NNB);
	for (int i = 0; i < NumTarget; ++i) {
		//std::cerr << "PID=" << h_target_list[i] << std::endl;
		for (int j = 0; j < NNB; ++j) {
			std::cerr << h_diff[_six*(i * NNB + j)] << " ";
		}
		std::cerr << std::endl;
	}
	*/

	//toHost(h_result  , d_result  , variable_size, stream);
	//toHost(h_neighbor, d_neighbor, variable_size, stream);


	cudaStreamSynchronize(stream); // Wait for all operations to finish

	toHost(h_result      , d_result      ,           _six*NumTarget);
	toHost(h_neighbor    , d_neighbor    , NumNeighborMax*NumTarget);
	toHost(h_num_neighbor, d_num_neighbor,                NumTarget);
	//printf("CUDA: transfer to host done\n");


	//cudaStreamSynchronize(stream); // Wait for all operations to finish


	for (int i=0;i<NumTarget;i++) {
		for (int j=0;j<h_num_neighbor[i];j++) {
			NeighborList[i][j] = h_neighbor[NumNeighborMax*i+j];
		}
		NumNeighbor[i] = h_num_neighbor[i];

		/*
		fprintf(stderr, "%d (%d) neighbors of %d = ", i, h_target_list[i], h_num_neighbor[i]);
		for (int j=0;j<h_num_neighbor[i];j++) {
			fprintf(stderr, "%d, ", NeighborList[i][j]);
		}
		fprintf(stderr, "\n");
		*/

		/*
		fprintf(stderr, "PID=%d: a=(%.4e,%.4e,%.4e), adot=(%.4e,%.4e,%.4e)\n",
				h_target_list[i],
				h_result[_six*i],
				h_result[_six*i+1],
				h_result[_six*i+2],
				h_result[_six*i+3],
				h_result[_six*i+4],
				h_result[_six*i+5]
				);
				*/
	}
	//fprintf(stderr, "\n");

	// out data
	for (int i=0; i<NumTarget; i++) {
		acc[i][0]  = h_result[_six*i];
		acc[i][1]  = h_result[_six*i+1];
		acc[i][2]  = h_result[_six*i+2];
		adot[i][0] = h_result[_six*i+3];
		adot[i][1] = h_result[_six*i+4];
		adot[i][2] = h_result[_six*i+5];
	}

	cublasDestroy(handle);
	/*
	my_free(h_background , d_background);
	my_free(h_result     , d_result);
	my_free(h_target     , d_target);
	my_free(h_neighbor   , d_neighbor);
	*/
	//cudaStreamDestroy(stream);
	//my_free_d(do_neighbor);
	//printf("CUDA: done?\n");
}










/*************************************************************************
 *	 Communication with HOST
 *************************************************************************/
void _ReceiveFromHost(
		int __NNB,
		CUDA_REAL m[],
		CUDA_REAL x[][3],
		CUDA_REAL v[][3],
		CUDA_REAL r2[],
		CUDA_REAL mdot[]
		){

	nbodymax = 100000000;
	_NNB     = __NNB;  // this values can be different from NNB 
	assert(_NNB <= nbodymax);
	cudaError_t cudaStatus;

	//printf("CUDA: receive starts\n");

	if ((first) || (new_size(_NNB) > variable_size )) {
		variable_size = new_size(_NNB);
		if (!first) {
			my_free(h_ptcl				 , d_ptcl);
			my_free(h_result       , d_result);
			my_free(h_neighbor     , d_neighbor);
			my_free(h_num_neighbor , d_num_neighbor);
			cudaFree(d_target);
			cudaFree(d_r2);
			cudaFree(d_diff);
			cudaFree(d_magnitudes);
		}
		else {
			first = false;
		}
		my_allocate(&h_ptcl         , &d_ptcl        ,         _seven*variable_size); // x,v,m
		my_allocate(&h_result       , &d_result      ,           _six*variable_size);
		my_allocate(&h_num_neighbor , &d_num_neighbor,                variable_size);
		my_allocate(&h_neighbor     , &d_neighbor    , NumNeighborMax*variable_size);
		cudaMalloc((void**)&d_r2        ,        variable_size * sizeof(CUDA_REAL));
		cudaMalloc((void**)&d_target    ,        variable_size * sizeof(int));
		cudaMalloc((void**)&d_diff      , _six * variable_size * variable_size * sizeof(CUDA_REAL));
		cudaMalloc((void**)&d_magnitudes, _two * variable_size * variable_size * sizeof(CUDA_REAL));
		//cudaMallocHost((void**)&h_diff          , _six * variable_size * variable_size * sizeof(CUDA_REAL));
		//cudaMallocHost((void**)&h_magnitudes    , _two * variable_size * variable_size * sizeof(CUDA_REAL));
	}


	for (int j=0; j<_NNB; j++) {
		for (int dim=0; dim<Dim; dim++) {
			h_ptcl[_seven*j+dim]   = x[j][dim];
			h_ptcl[_seven*j+dim+3] = v[j][dim];
		}
		h_ptcl[_seven*j+6] = m[j];
		//h_particle[j].setParticle(m[j], x[j], v[j], r2[j], mdot[j]);
	}

	//toDevice(h_background,d_background,variable_size);
	toDevice(h_ptcl,d_ptcl, _seven*_NNB, stream);
	toDevice(r2    ,d_r2  ,        _NNB, stream);
	//fprintf(stdout, "CUDA: receive done\n");
}





void _InitializeDevice() {

	if(is_open){
		fprintf(stderr, "it is already open\n");
		return;
	}
	is_open = true;

	std::cout << "Initializing CUDA ..." << std::endl;
	int device = 0; 
	int deviceCount;

	cudaGetDeviceCount(&deviceCount);
	cudaStreamCreate(&stream);

	std::cout << "There are " << deviceCount << " GPUs." << std::endl;


	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, devid);
	//  char *hostname = getenv("HOSTNAME");


	fprintf(stderr, "-------- GPU initialization : device: %s (%d)----------\n", prop.name, devid);


	cudaSetDevice(device);


	is_open = true;

	std::cout << "CUDA initialized successfully!" << std::endl;
}



void _CloseDevice() {
	if(!is_open) {
		fprintf(stderr, "gpunb: it is already close\n");
		return;
	}
	is_open = false;


	cudaError_t error;


	printf("CUDA: ?!! ...\n");
	//my_free(&h_result    , &d_result);
	fprintf(stderr, "result ...\n");
	//my_free(&h_target    , &d_target);
	fprintf(stderr, "target ...\n");
	//my_free(&h_neighbor  , &d_neighbor);
	fprintf(stderr, "neighbor ...\n");
	//my_free(&h_background, &d_background);

}





extern "C" {
	void InitializeDevice() {
		_InitializeDevice();
	}
	void CloseDevice(){
		_CloseDevice();
	}
	void SendToDevice(int *__NNB, CUDA_REAL m[], CUDA_REAL x[][3], CUDA_REAL v[][3], CUDA_REAL r2[], CUDA_REAL mdot[]) {
		_ReceiveFromHost(*__NNB, m, x, v, r2, mdot);
	}
	void CalculateAccelerationOnDevice(int *NumTarget, int *h_target_list, CUDA_REAL acc[][3], CUDA_REAL adot[][3], int NumNeighbor[], int **NeighborList) {
		GetAcceleration(*NumTarget, h_target_list, acc, adot, NumNeighbor, NeighborList);
	}
}

