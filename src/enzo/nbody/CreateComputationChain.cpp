#include "global.h"
#include <iostream>
#include <vector>
#define no_debug

void merge(int* index, double *timesteps, int left, int mid, int right);
void mergeSort(int *index, double *timesteps, int left, int right);


int createComputationChain(std::vector<Particle*> &particle) {
	// This stores the first particle of each level in the particle chain

	int *index;
	double *timesteps;
	timesteps = new double[NNB];
	index = new int[NNB];



	for (int i=0; i<NNB; i++) {
		index[i] = i;
		timesteps[i] = particle[i]->TimeStepIrr;
	}



#ifdef debug
	std::cout << "Index" << '\n';
	for (int i=0; i<NNB; i++) {
		std::cout << index[i] << ' ';
	}
	std::cout << '\n';

	std::cout << "Timesteps" << '\n';
	for (int i=0; i<NNB; i++) {
		std::cout << timesteps[i] << ' ';
	}
	std::cout << '\n';
#endif


	mergeSort(index, timesteps, 0, NNB-1);

	LevelList.push_back(index[0]);
	for (int i=0; i<NNB-1; i++) {
		if (timesteps[i] != timesteps[i+1])
			LevelList.push_back(index[i+1]);
	}



#ifdef debug
	std::cout << "Index" << '\n';
	for (int i=0; i<NNB; i++) {
		std::cout << index[i] << ' ';
	}
	std::cout << '\n';

	std::cout << "LevelList" << '\n';
	for (int element : LevelList) {
		    std::cout << element << " ";
	}
	std::cout << '\n';

	std::cout << "Timesteps" << '\n';
	for (int i=0; i<NNB; i++) {
		std::cout << timesteps[i] << ' ';
	}
	std::cout << '\n';

#endif


	/*
	for (Particle* elem: particle) {
		for (int i=0; i<NNB-1; i++) {
			if (elem->getPID() == index[i])
				elem->NextParticle = particle[index[i+1]];
		}
	}
	*/

	delete [] index;
	//delete [] timesteps;
	return DONE;
}

void merge(int* index, double *timesteps, int left, int mid, int right) {
	int n1 = mid - left + 1;
	int n2 = right - mid;

	//Create temporary arrays
	std::vector<int> L1(n1), R1(n2);
	std::vector<double> L2(n1), R2(n2);

	// Copy data to temp arrays L[] and R[]
	for (int i = 0; i < n1; ++i) {
		L1[i] = index[left + i];
		L2[i] = timesteps[left + i];
	}
	for (int j = 0; j < n2; ++j) {
		R1[j] = index[mid + 1 + j];
		R2[j] = timesteps[mid + 1 + j];
	}

	// Merge the temp arrays back into arr[left..right]
	int i = 0, j = 0, k = left;
	while (i < n1 && j < n2) {
		if (L2[i] <= R2[j]) {
			index[k]     = L1[i];
			timesteps[k] = L2[i];
			++i;
		} else {
			index[k]     = R1[j];
			timesteps[k] = R2[j];
			++j;
		}
		++k;
	}

	// Copy the remaining elements of L[], if any
	while (i < n1) {
		index[k]     = L1[i];
		timesteps[k] = L2[i];
		++i;
		++k;
	}

	// Copy the remaining elements of R[], if any
	while (j < n2) {
		index[k]     = R1[j];
		timesteps[k] = R2[j];
		++j;
		++k;
	}
}

void mergeSort(int *index, double *timesteps, int left, int right) {
	if (left < right) {
		int mid = left + (right - left) / 2;

		//Sort first and second halves
		mergeSort(index, timesteps, left, mid);
		mergeSort(index, timesteps, mid + 1, right);

		// Merge the sorted halves
		merge(index, timesteps, left, mid, right);
	}
}
