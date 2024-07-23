#pragma once
#include "../defs.h"
#include "cuda_defs.h"
extern "C" {
	void InitializeDevice();
	void CloseDevice();
	void ProfileDevice(int *irank);
	void SendToDevice(int *_NNB, CUDA_REAL m[], CUDA_REAL x[][3], CUDA_REAL v[][3], CUDA_REAL r[], CUDA_REAL mdot[]);
	void CalculateAccelerationOnDevice(int *NumTarget, int *h_target_list, CUDA_REAL acc[][3], CUDA_REAL adot[][3], int NumNeighbor[], int **NeighborList);
}
