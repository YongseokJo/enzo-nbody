#pragma once
extern "C" {
	void InitializeDevice(int *irank);
	void OpenDevice(const int *irank);
	void CloseDevice();
	void ProfileDevice(int *irank);
	void SendToDevice(int *_NNB, double m[], double x[][3], double v[][3], double mdot[], int *_NumNeighborMax);
	void CalculateAccelerationOnDevice(int *NumTarget, int *offset, double x[][3], double v[][3], double acc[][3], double adot[][3], double mdot[], double radius[], int NumNeighbor[], int **NeighborList);
}
