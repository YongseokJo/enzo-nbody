#include "../defs.h"
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
//#include "cuda_routines.h"





#ifdef unuse
void initializeCudaAndCublas(cublasHandle_t* handle) {
	cudaError_t cudaStat = cudaSetDevice(0);
	if (cudaStat != cudaSuccess) {
		std::cerr << "cudaSetDevice failed!" << std::endl;
		exit(1);
	}

	cublasStatus_t stat = cublasCreate(handle);
	if (stat != CUBLAS_STATUS_SUCCESS) {
		std::cerr << "CUBLAS initialization failed!" << std::endl;
		exit(1);
	}
}

void checkCudaError(cudaError_t result) {
	if (result != cudaSuccess) {
		std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << std::endl;
		exit(EXIT_FAILURE);
	}
}

template <typename T>
void my_allocate(T **host, T **device, const int size) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(device, size*sizeof(T));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	//*host = (T*)calloc(size, sizeof(T));
	//*host = (T*)malloc(size*sizeof(T));
	/*
		 if (host == NULL) {
		 fprintf(stderr, "Memory allocation failed\n");
		 }
	 */
	//host = new T[size]();
	cudaStatus = cudaMallocHost(host, size*sizeof(T));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMallocHost failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}

template <typename T>
void my_free(T &host, T &device) {
	cudaError_t cudaStatus;
	cudaStatus = cudaFree(device);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaFreeHost(host);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	//free(host);
}



template <typename T>
void my_allocate_d(T **device, const int size) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(device, size*sizeof(T));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}

template <typename T>
void my_free_d(T &device) {
	cudaError_t cudaStatus;
	cudaStatus = cudaFree(device);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}




template <typename T>
void toDevice(T *host, T *device, const int size, cudaStream_t &stream) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpyAsync(device, host, size * sizeof(T), cudaMemcpyHostToDevice, stream);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyHostToDevice failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}

template <typename T>
void toHost(T *host, T *device, const int size, cudaStream_t &stream) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpyAsync(host, device, size * sizeof(T), cudaMemcpyDeviceToHost, stream);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyDeviceToHost failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}


template <typename T>
void toDevice(T *host, T *device, const int size) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(device, host, size * sizeof(T), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyHostToDevice failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}

template <typename T>
void toHost(T *host, T *device, const int size) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(host, device, size * sizeof(T), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyDeviceToHost failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}
#endif


