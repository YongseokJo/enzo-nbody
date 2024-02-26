#pragma once
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <cutil.h>
#include <cuda_runtime.h>

//#define MaxNeighbor 100
#define NAN_CHECK(val) assert((val) == (val));
//#define CUDA_SAFE_CALL checkCudaErrors


template <typename T>
void my_allocate(T **host, T **device, const int size) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(device, size*sizeof(T));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	*host = (T*)calloc(size, sizeof(T));
	//*host = (T*)malloc(size*sizeof(T));
	if (host == NULL) {
		fprintf(stderr, "Memory allocation failed\n");
	}
	//host = new T[size]();
	//cudaStatus = cudaMallocHost(&host, size*sizeof(T));
	//if (cudaStatus != cudaSuccess) {
		//fprintf(stderr, "cudaMallocHost failed: %s\n", cudaGetErrorString(cudaStatus));
	//}
}

template <typename T>
void my_free(T *host, T *device) {
	cudaError_t cudaStatus;
	cudaStatus = cudaFree(device);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	free(host);
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






struct TargetParticle{
	float3 pos;
	float  mdot;
	float3 vel;
	float  r2; // for AC neighbor

	void setParticle(double _mdot, double x[3], double v[3], double _r2){
		mdot   = _mdot;
		pos.x  = x[0];
		pos.y  = x[1];
		pos.z  = x[2];
		vel.x  = v[0];
		vel.y  = v[1];
		vel.z  = v[2];
		r2 = _r2;

		NAN_CHECK(x[0]);
		NAN_CHECK(x[1]);
		NAN_CHECK(x[2]);
		NAN_CHECK(_mdot);
		NAN_CHECK(v[0]);
		NAN_CHECK(v[1]);
		NAN_CHECK(v[2]);
		NAN_CHECK(_r2);
  }
};



struct BackgroundParticle{
  float3 pos;
  float  mass;
  float3 vel;
	float  mdot;

	//BackgroundParticle(int) {}
	BackgroundParticle(double m, double x[3], double v[3], double _mdot){
		mass  = m;
		pos.x = x[0];
    pos.y = x[1];
    pos.z = x[2];
    vel.x = v[0];
    vel.y = v[1];
    vel.z = v[2];
		mdot  = _mdot;

    NAN_CHECK(x[0]);
    NAN_CHECK(x[1]);
    NAN_CHECK(x[2]);
    NAN_CHECK(m);
    NAN_CHECK(v[0]);
    NAN_CHECK(v[1]);
    NAN_CHECK(v[2]);
    NAN_CHECK(_mdot);
  }

  void setParticle(double m, double x[3], double v[3], double _mdot){
		mass  = m;
		pos.x = x[0];
		pos.y = x[1];
		pos.z = x[2];
		vel.x = v[0];
		vel.y = v[1];
		vel.z = v[2];
		mdot  = _mdot;

		NAN_CHECK(x[0]);
		NAN_CHECK(x[1]);
		NAN_CHECK(x[2]);
		NAN_CHECK(m);
		NAN_CHECK(v[0]);
		NAN_CHECK(v[1]);
		NAN_CHECK(v[2]);
		NAN_CHECK(_mdot);
	}
  //__device__ BackgroundParticle() {}
};



struct Result{
	float3 acc;
	float3 adot;
	//unsigned short num_ac;          //  8 words
	//unsigned short ac_list[MaxNeighbor];

	void clear_h(void) {
		acc.x  = acc.y  = acc.z  = 0.f;
		adot.x = adot.y = adot.z = 0.f;
		//nnb = 0;
	}

	__device__  void clear() {
		acc.x  = acc.y  = acc.z  = 0.f;
		adot.x = adot.y = adot.z = 0.f;
	}

	/*
	__device__ void operator+=(const Result &rhs){
		acc.x  += rhs.acc.x;
		acc.y  += rhs.acc.y;
		acc.z  += rhs.acc.z;
		adot.x += rhs.adot.x;
		adot.y += rhs.adot.y;
		adot.z += rhs.adot.z;
	}
	*/
};

struct  Neighbor{
	int NumNeighbor;
	int NeighborList[100]; // this needs to be modified.


	__device__ void clear() {
		NumNeighbor = 0;
#pragma unroll
		for (int i=0; i<100; i++) {
			NeighborList[i] = 0; 
		}
	}
	/*
  Neighbor(){
		NumNeighbor = 0;
		for (int i=0; i<100; i++) {
			NeighborList[i] = -1;
		}
  }
	*/
};


/*
struct  Neighbor{
	int width = 2;
	int height = 100;

	int NumNeighbor[2];
	int NeighborList[2][100]; // this needs to be modified.

	int* NeighborList_d;
	int* NumNeighbor_d;
	size_t pitch;

  Neighbor(){
		cudaError_t cudaStatus;
		cudaStatus = cudaMallocPitch(&NeighborList_d, &pitch, width * sizeof(int), height);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaMallocPitch failed: %s\n", cudaGetErrorString(cudaStatus));
		}
		cudaStatus = cudaMalloc(&NumNeighbor_d, width * sizeof(int));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
		}
  }
	void toHost() {
		cudaError_t cudaStatus;
		cudaStatus = cudaMemcpy(NumNeighbor, NumNeighbor_d, width * sizeof(int), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaMemcpy failed: %s\n", cudaGetErrorString(cudaStatus));
		}
		cudaStatus = cudaMemcpy2D(NeighborList, width * sizeof(int), NeighborList_d, pitch, width * sizeof(int), height, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaMemcpy2d failed: %s\n", cudaGetErrorString(cudaStatus));
		}
	}

	void toDevice() {
		cudaError_t cudaStatus;
		cudaStatus = cudaMemcpy(NumNeighbor_d, NumNeighbor, width * sizeof(int), cudaMemcpyHostToDevice);
		cudaStatus = cudaMemcpy2D(NeighborList_d, pitch, NeighborList, width * sizeof(int), width * sizeof(int), height, cudaMemcpyHostToDevice);
	}

	void free() {
		cudaError_t cudaStatus;
		cudaStatus = cudaFree(NeighborList_d);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaFree list failed: %s\n", cudaGetErrorString(cudaStatus));
		}
		cudaStatus = cudaFree(NumNeighbor_d);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaFree num failed: %s\n", cudaGetErrorString(cudaStatus));
		}
	}
};
*/
