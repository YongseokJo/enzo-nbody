#include <algorithm>
#include <iostream>
#include <vector>
#include "global.h"
#define debug

void Merge(std::vector<int> index, std::vector<double> timesteps, int left, int mid, int right);
void MergeSort(std::vector<int> index, std::vector<double> timesteps, int left, int right);


// Function to perform argsort on a vector
bool SortComputationChain(std::vector<Particle*> ComputationChainTmp) {

	std::vector<int> index{};
	std::vector<double> time{};
	double NextIrrTime = 0.0;

	int i=0;
	for (Particle *ptcl : ComputationChainTmp)
	{
		NextIrrTime = ptcl->CurrentTimeIrr + ptcl->TimeStepIrr;
		//std::cout << NextIrrTime << " ";
		if ((ptcl->NumberOfAC != 0) && (NextIrrTime <= NextRegTime)) {
			index.push_back(i);
			time.push_back(NextIrrTime);
		}
		i++;
	}
	//std::cout << std::endl;

	if (index.size() == 0) {
		return false;
	}

	// Sort the index array based on the values in the original vector
	std::sort(index.begin(), index.end(), [&time](int i1, int i2) {return time[i1] < time[i2];});

	/*
	ComputationChain.clear();
	for (int ind: index) {
		ComputationChain.push_back(ComputationChainTmp[ind]);
	}
	*/

	return true;
}


// Function to perform argsort on a vector
Particle *SortComputationChain(Particle* ptcl) {

	Particle *NextParticle, *PreviousParticle, *NextComputation;
	double NextIrrTime = 0.0, NextParticleNextIrrTime=0.0;
	
	NextComputation = ptcl->NextParticleForComputation;

	NextIrrTime = ptcl->CurrentTimeIrr + ptcl->TimeStepIrr;

	// this paticle's reached NextRegTime.
	if ((ptcl->NumberOfAC == 0) || (NextIrrTime > NextRegTime)) {
		//return false;
		return NextComputation;
	}

	// if there's only one particle left
	if (NextComputation == nullptr && NextIrrTime <= NextRegTime) {
		return ptcl;
	}

	PreviousParticle = ptcl;
	NextParticle     = NextComputation;
	while (NextParticle != nullptr) {
		NextParticleNextIrrTime = NextParticle->CurrentTimeIrr + NextParticle->TimeStepIrr;
		if (NextIrrTime <= NextParticleNextIrrTime) {
			// This part should be improved.
			if (PreviousParticle == ptcl) {
				ptcl->NextParticleForComputation = NextParticle;
				NextComputation = ptcl;
			} 
			else 
			{
				PreviousParticle->NextParticleForComputation = ptcl;
				ptcl->NextParticleForComputation = NextParticle;
			}
			break;
		}
		PreviousParticle = NextParticle;
		NextParticle     = NextParticle->NextParticleForComputation;
	}

	if ((NextParticle == nullptr) && NextIrrTime <= NextRegTime) {
		PreviousParticle->NextParticleForComputation = ptcl; 
		ptcl->NextParticleForComputation = nullptr;
	}
		



	NextParticle = NextComputation; 
	//std::cout << "Time:";
	while (NextParticle != nullptr) {
		NextParticleNextIrrTime = NextParticle->CurrentTimeIrr + NextParticle->TimeStepIrr;
		//std::cout << NextParticleNextIrrTime << '(' << NextParticle->PID << ')' << ' ';
		NextParticle = NextParticle->NextParticleForComputation;
	}
	//std::cout << std::endl;


	return NextComputation;
}


bool CreateComputationChain(std::vector<Particle*> &particle) {

	std::vector<int> index{};
	std::vector<int> sorted_index{};
	std::vector<double> time{};
	double NextIrrTime = 0.0;

	int i=0;
	//std::cout << "NextIrrTime:\n" << std::endl;
	for (Particle *ptcl : particle)
	{
		// advance irregular time without irregular routine
		NextIrrTime = ptcl->CurrentTimeIrr + ptcl->TimeStepIrr;
		//std::cout << NextIrrTime << " ";
		if ((ptcl->NumberOfAC != 0) && (NextIrrTime <= NextRegTime)) {
			index.push_back(i);
			time.push_back(NextIrrTime);
		}
		i++;
	}
	//std::cout << std::endl;

	for (i=0; i<index.size(); i++) {
		sorted_index.push_back(i);
	}

	if (index.size() == 0) {
		return false;
	}

	/*
	std::cout << "PID" << '\n';
	for (int ind: sorted_index) {
		//std::cout << ind << ' ';
		std::cout << particle[index[ind]]->PID << ' ';
	}
	std::cout << '\n';
	std::cout << "Index" << '\n';
	for (int ind: sorted_index) {
		std::cout << ind << ' ';
	}
	std::cout << '\n';
	std::cout << "Timesteps" << '\n';
	for (int ind: sorted_index) {
		NextIrrTime = particle[index[ind]]->CurrentTimeIrr + particle[index[ind]]->TimeStepIrr;
		std::cout << NextIrrTime << ' ';
	}
	std::cout << '\n';
	*/

	// Sort the index array based on the values in the original vector
	std::sort(sorted_index.begin(), sorted_index.end(), [&time](int i1, int i2) {return time[i1] < time[i2];});

	/*
	std::cout << "ordered PID" << '\n';
	for (int ind: sorted_index) {
		//std::cout << ind << ' ';
		std::cout << particle[index[ind]]->PID << ' ';
	}
	std::cout << '\n';
	std::cout << "ordered Index" << '\n';
	for (int ind: sorted_index) {
		std::cout << ind << ' ';
	}
	std::cout << '\n';
	std::cout << "ordered Timesteps" << '\n';
	*/
	/*
	for (int ind: sorted_index) {
		NextIrrTime = particle[index[ind]]->CurrentTimeIrr + particle[index[ind]]->TimeStepIrr;
		std::cout << NextIrrTime << ' ';
	}
	*/
	//std::cout << '\n';


	//std::cout << "Ordered NextIrrTime: " << std::flush;
	Particle* NextParticle = nullptr;
	int ind;
	for (int i=index.size()-1; i>=0; i--) {
		ind = index[sorted_index[i]];
		NextIrrTime = particle[ind]->CurrentTimeIrr + particle[ind]->TimeStepIrr;
		//std::cout << NextIrrTime << "(" << particle[ind]->PID << ")" <<  " ";
		particle[ind]->NextParticleForComputation = NextParticle;
		NextParticle = particle[ind];
	}
	FirstComputation = NextParticle;
	//std::cout << std::endl;

	/*
		 ComputationChain.clear();
		 for (int ind: index) {
		 ComputationChain.push_back(particle[ind]);
		 }
		 */

	index.clear();
	time.clear();

	return true;
}

/*
bool CreateComputationChain(std::vector<Particle*> &particle, std::vector<Particle*> &ComputationChainTmp) {
	// This stores the first particle of each level in the particle chain

	double NextIrrTime = 0.0;
	std::vector<double >timesteps{};
	std::vector<int> index{};

	for (int i=0; i<NNB; i++) {
		NextIrrTime = particle[i]->CurrentTimeIrr+particle[i]->TimeStepIrr;
		if (particle[i]->NumberOfAC == 0|| NextIrrTime > NextRegTime)
			continue;
		timesteps.push_back(NextIrrTime);
		index.push_back(i);
	}

	if (index.size() == 0) {
		return false;
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


  MergeSort(index, timesteps, 0, index.size()-1);

	ComputationChainTmp.push_back(particle[index[0]]);
	for (int i=0; i<index.size()-1; i++) {
		if (timesteps[i] != timesteps[i+1])
			ComputationChainTmp.push_back(particle[index[i+1]]);
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


	for (Particle* elem: particle) {
		for (int i=0; i<NNB-1; i++) {
			if (elem->getPID() == index[i])
				elem->NextParticle = particle[index[i+1]];
		}
	}

	return true;
}
*/


void Merge(std::vector<int> index, std::vector<double> timesteps, int left, int mid, int right) {
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

void MergeSort(std::vector<int> index, std::vector<double> timesteps, int left, int right) {
	if (left < right) {
		int mid = left + (right - left) / 2;

		//Sort first and second halves
		MergeSort(index, timesteps, left, mid);
		MergeSort(index, timesteps, mid + 1, right);

		// Merge the sorted halves
		Merge(index, timesteps, left, mid, right);
	}
}
