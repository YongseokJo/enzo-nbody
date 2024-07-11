#include <algorithm>
#include <iostream>
#include <vector>
#include "global.h"
#define debug

void Merge(std::vector<int> index, std::vector<double> timesteps, int left, int mid, int right);
void MergeSort(std::vector<int> index, std::vector<double> timesteps, int left, int right);


// Function to perform argsort on a vector
bool CreateComputationList(Particle* ptcl) {

	if (ptcl == nullptr) {
		return false;
	}

	Particle *NextParticle=ptcl->NextParticleForComputation;
	ULL ThisParticleNextIrrBlock = 0;
	ULL NextParticleNextIrrBlock = 0;

	ComputationList.clear();
	ComputationList.push_back(ptcl);
	ThisParticleNextIrrBlock = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr;

	if (ThisParticleNextIrrBlock <= global_time_irr) {
		fprintf(stderr, "Particle %d's next time (%.3e) is smaller than global time of %.3e.\n",
				ptcl->PID, ThisParticleNextIrrBlock*time_step, global_time_irr*time_step);
		throw std::runtime_error("Fatal Error in CreateComputationList\n");	
	}

	while (NextParticle != nullptr) {

		NextParticleNextIrrBlock = NextParticle->CurrentBlockIrr + NextParticle->TimeBlockIrr;

		if (NextParticleNextIrrBlock <= global_time_irr) {
			fprintf(stderr, "Particle %d's next time (%.3e) is smaller than global time of %.3e.\n",
				 	NextParticle->PID, NextParticle->CurrentTimeIrr+NextParticle->TimeStepIrr, global_time_irr*time_step);
			throw std::runtime_error("Fatal Error in CreateComputationList\n");	
		}

		if (NextParticleNextIrrBlock == ThisParticleNextIrrBlock) {
			ComputationList.push_back(NextParticle);
		}
		else if (NextParticleNextIrrBlock > ThisParticleNextIrrBlock){
			FirstComputation = NextParticle;
			break;
		}
		else if (NextParticleNextIrrBlock < ThisParticleNextIrrBlock) {
			throw std::runtime_error("NextParticleNextIrrBlock is smaller then ThisParticleNextIrrBlock.\nThis might indicate that the CompuationChain must be ill-ordered.\n");	
		}
		else {
			throw std::runtime_error("NextIrrBlock is contaminated somehow.\n");	
		}

		NextParticle = NextParticle->NextParticleForComputation;
	}

	if (NextParticle == nullptr) {
		//std::cout << "NextParticle is null." << std::endl;
		FirstComputation = nullptr;
	}

	return true;
}

// Function to perform argsort on a vector
bool UpdateComputationChain(Particle* ptcl) {

	Particle *ThisParticle, *NextParticle, *PreviousParticle;
	ULL ThisParticleNextIrrBlock= 0, NextParticleNextIrrBlock =0;

	ThisParticle = ptcl;
	NextParticle = FirstComputation;
	ThisParticleNextIrrBlock = ThisParticle->CurrentBlockIrr + ThisParticle->TimeBlockIrr; // of this particle

	//std::cout << "in sort, ("<< ThisParticle->PID<<") NextIrrBlock=" << NextIrrBlock << std::endl;


	// this paticle's reached NextRegTime and this particle will be removed from computationlist
	if ((ThisParticle->NumberOfAC == 0) || (ThisParticleNextIrrBlock > NextRegTimeBlock)) { 
		ThisParticle->NextParticleForComputation = nullptr;
		return true;
	}

	// if there's only one particle left and if that should go on
	if (FirstComputation == nullptr && ThisParticleNextIrrBlock <= NextRegTimeBlock) {
		//return nullptr;
		FirstComputation = ThisParticle;
		ThisParticle->NextParticleForComputation = nullptr;
		return true;
	}


	//PreviousParticle = ThisParticle;
	PreviousParticle = nullptr;
	while (NextParticle != nullptr) {
		NextParticleNextIrrBlock = NextParticle->CurrentBlockIrr + NextParticle->TimeBlockIrr; // of Next particle
		// this is where this particle should be located
		if (ThisParticleNextIrrBlock <= NextParticleNextIrrBlock) {
			// where this is first loop
			if (NextParticle == FirstComputation) {
				ThisParticle->NextParticleForComputation = FirstComputation;
				FirstComputation = ThisParticle;
			}
			else
			{
				PreviousParticle->NextParticleForComputation = ThisParticle;
				ThisParticle->NextParticleForComputation = NextParticle;
			}
			break;
		}
		PreviousParticle = NextParticle;
		NextParticle     = NextParticle->NextParticleForComputation;
	}

	if ((NextParticle == nullptr) && ThisParticleNextIrrBlock <= NextRegTimeBlock) {
		PreviousParticle->NextParticleForComputation = ThisParticle; 
		ThisParticle->NextParticleForComputation = nullptr;
	}

	/*
	NextParticle = FirstComputation; 
	//std::cout << "Time:";

	while (NextParticle != nullptr) {
		NextParticleNextIrrBlock = NextParticle->CurrentBlockIrr + NextParticle->TimeBlockIrr;
		//std::cout << NextParticleNextIrrBlock << '(' << NextParticle->PID << ')' << ' ';
		NextParticle = NextParticle->NextParticleForComputation;
	}
	//std::cout << std::endl;
	*/

	return true;
}




bool CreateComputationChain(std::vector<Particle*> &particle) {

	fprintf(nbpout, "CreateComputationList ...\n");
	fflush(nbpout);
	std::vector<int> index{};
	std::vector<int> sorted_index{};
	std::vector<double> time{};
	ULL NextIrrBlock = 0;

	int i=0;
	//std::cout << "NextIrrTime:\n" << std::endl;
	for (Particle *ptcl : particle)
	{
		// advance irregular time without irregular routine
		NextIrrBlock = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr;
		//std::cout << NextIrrBlock << "(" << ptcl->PID << ")" <<  " ";
		if ((ptcl->NumberOfAC != 0) && (NextIrrBlock <= NextRegTimeBlock)) {
			index.push_back(i);
			time.push_back(NextIrrBlock);
		}
		i++;
	}
	//std::cout << std::endl;
	
	if (index.size() == 0) {
		return false;
	}

	for (i=0; i<index.size(); i++) {
		sorted_index.push_back(i);
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


	//std::cerr << "Ordered NextIrrTime: " << std::flush;
	Particle* NextParticle = nullptr;
	int ind;
	for (int i=index.size()-1; i>=0; i--) {
		ind = index[sorted_index[i]];
		NextIrrBlock = particle[ind]->CurrentBlockIrr + particle[ind]->TimeBlockIrr;
		//std::cout << NextIrrBlock << "(" << particle[ind]->PID << ")" <<  " ";
		particle[ind]->NextParticleForComputation = NextParticle;
		NextParticle = particle[ind];
	}
	//std::cout<<std::endl;

	FirstComputation = NextParticle;
	while (NextParticle != nullptr) {
		NextIrrBlock = NextParticle->CurrentBlockIrr + NextParticle->TimeBlockIrr;
		//std::cerr << NextIrrBlock*time_step << "(" << NextParticle->PID << ")" <<  " ";
		NextParticle = NextParticle->NextParticleForComputation;
	}
	//std::cerr << std::endl;

	//fprintf(stderr, "Number of Irregular Particles = %lu\n", index.size());

	/*
		 ComputationChain.clear();
		 for (int ind: index) {
		 ComputationChain.push_back(particle[ind]);
		 }
		 */

	index.clear();
	time.clear();


	fprintf(nbpout, "CreateComputationList done (size=%lu).\n",ComputationList.size());
	fflush(nbpout);

	return true;
}



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
