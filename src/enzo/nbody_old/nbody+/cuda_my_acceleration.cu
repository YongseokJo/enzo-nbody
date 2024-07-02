#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <cassert>
#include <cuda_runtime.h>
#include "cuda_types.h"
//#include "cuda_global.h"
//#include "cuda_functions.h"
//#include "cuda_routines.h"

#define _PROFILE

//#define THREAD 1024 // 2048 for A100
//#define BLOCK 32    // 32 for A100 

#define THREAD 128 // 2048 for A100
#define BLOCK 128    // 32 for A100 

#define ESP2 1e-6
#define new_size(A) (A > 1024) ? int(pow(2,ceil(log(A)/log(2.0)))) : 1024


static int NNB;
static int NumNeighborMax;
static double time_send, time_grav, time_out, time_nb;
static long long numInter;
static int icall,ini,isend;
static int nbodymax;
static int devid, numGPU;
static bool is_open = false;
static bool devinit = false;
static bool first   = true;
static int variable_size;
const int memory_size = 512;
//BackgroundParticle *h_background, *d_background;
BackgroundParticle *h_background; //, *background;
BackgroundParticle *d_background;
Result *h_result, *d_result;
TargetParticle *h_target, *d_target;
Neighbor *h_neighbor, *d_neighbor;
Neighbor *do_neighbor;


/*************************************************************************
 *	 Computing Acceleration
 *************************************************************************/
__global__ void CalculateAcceleration(
		const int NNB,
		const int NumTarget,
		const int tg_offset,
		const TargetParticle target[],
		const BackgroundParticle background[],
		Result result[],
		Neighbor neighbor[]
		);
__device__ void kernel(
		const TargetParticle &i,
		const BackgroundParticle &j,
		Result &res,
		Neighbor &neighbor,
		int &bg_index,
		const int &tg_index
		);

__global__ void OrganizeNeighbor(const Neighbor do_neighbor[], Neighbor d_neighbor[],
		const int offset, const int NNB);

__device__ void initializeResult(Result &res);
__device__ void _addition(Result &result, const Result res);
__device__ void _copy(Result &result, const Result res);

void GetAcceleration(
		int NumTarget,
		double x[][3],
		double v[][3],
		double acc[][3],
		double adot[][3],
		double mdot[],
		double r2[],
		int NumNeighbor[],
		int **NeighborList
		) {
	icall++;
	assert(is_open);
	assert((NumTarget > 0) && (NumTarget <= NNB));

	cudaError_t cudaStatus;
	cudaError_t error;
	int NumTarget_local = 0;



	for(int i=0; i<NumTarget; i++) {
		//printf("CUDA: x=%.3e, y=%.3e, r=%.3e\n", x[offset+i][0], x[offset+i][1], r2[offset+i]);
		h_result[i].clear_h();
		h_neighbor[i].clear_h();
		h_target[i].setParticle(mdot[i], x[i], v[i], r2[i]);
		//printf("CUDA: res acc x=%.3e, y=%.3e\n", h_result[i].acc.x, h_result[i].acc.y);
	}

	// actually all the data is on the memory already we can just pass the indices
	toDevice(h_target  , d_target  , variable_size); 
	////printf("CUDA: transfer done\n");


	dim3 thread_size(THREAD, 1, 1);
	dim3 block_size(BLOCK,1,1);

	//for (int offset=0; offset<NumTarget; offset+=memory_size) {
	//	NumTarget_local = std::min(memory_size,NumTarget-offset);

	for (int tg_offset = 0; tg_offset < NumTarget; tg_offset += BLOCK) {
		CalculateAcceleration <<< block_size, thread_size >>>
			(NNB, NumTarget, tg_offset, d_target, d_background, d_result, do_neighbor);
	} // endfor i, target particles

	for (int nn_offset=0; (nn_offset)*BLOCK*THREAD<NumTarget; nn_offset++) {
		OrganizeNeighbor <<< block_size, thread_size >>>
			(do_neighbor, d_neighbor, nn_offset, NumTarget);
	}
		//printf("CUDA: calculation done\n");
	//} // endfor offset


	//printf("CUDA: neighbor post processing done\n");

	toHost(h_result  , d_result  , variable_size);
	toHost(h_target  , d_target  , variable_size);
	toHost(h_neighbor, d_neighbor, variable_size);
	cudaDeviceSynchronize();
	//printf("CUDA: transfer to host done\n");

	for (int i=0;i<NumTarget;i++) {
		for (int j=0;j<h_neighbor[i].NumNeighbor;j++) {
			NeighborList[i][j] = h_neighbor[i].NeighborList[j];
		}
		NumNeighbor[i] = h_neighbor[i].NumNeighbor;
	}

	// out data
	for (int i=0; i<NumTarget; i++) {
		acc[i][0]  = h_result[i].acc.x;
		acc[i][1]  = h_result[i].acc.y;
		acc[i][2]  = h_result[i].acc.z;
		adot[i][0] = h_result[i].adot.x;
		adot[i][1] = h_result[i].adot.y;
		adot[i][2] = h_result[i].adot.z;
	}

	//printf("CUDA: done?\n");
}


// each block is assigned with one target particles and N_thread background particles
__global__ void CalculateAcceleration(
		const int NNB,
		const int NumTarget,
		const int tg_offset,
		const TargetParticle target[],
		const BackgroundParticle background[],
		Result result[],
		Neighbor neighbor[]
		) {

	int tg_index = blockIdx.x + tg_offset;
	if (tg_index >= NumTarget)
		return;

	int bg_index;
	__shared__ Result res[THREAD];
	res[threadIdx.x].clear();
	Neighbor nb;
	nb.clear();
	neighbor[tg_index*blockDim.x+threadIdx.x].clear();
	result[tg_index].clear();


	// looping over N*BlockDim.x+threadId.x;
	for (int j = 0; j < NNB; j += blockDim.x) {
		bg_index = threadIdx.x + j;
		//printf("CUDA: 1. (%d,%d), res=%e\n", tg_index, bg_index, res[threadIdx.x].acc.x);

		if (bg_index < NNB) {
			//printf("CUDA: 2. (%d,%d), res=%e\n", tg_index, bg_index, res[threadIdx.x].acc.x);
			kernel(target[tg_index], background[bg_index], res[threadIdx.x], nb, bg_index, tg_index);

			neighbor[tg_index*blockDim.x+threadIdx.x].NumNeighbor = nb.NumNeighbor;
			for (int i=0; i<nb.NumNeighbor; i++) 
				neighbor[tg_index*blockDim.x+threadIdx.x].NeighborList[i] = nb.NeighborList[i];

			//neighbor_num[tg_index*blockDim.x+threadIdx.x].NumNeighbor++;
			//printf("CUDA: 3. (%d,%d), res=%e\n", threadIdx.x, blockIdx.x, res[threadIdx.x].acc.x);
			//tg_index, bg_index, res[threadIdx.x].acc.x);
		}
		else
		{
			break;
		}
	} //endfor j, backgroun particles

	//printf("CUDA: 4. (%d,%d), res=%e\n", tg_index, bg_index, res[threadIdx.x].acc.x);
	// Reduction in shared memory
	for (int stride = 1; stride < blockDim.x; stride *= 2) {
		int index = 2 * stride * threadIdx.x;
		if (index < blockDim.x && bg_index < NNB) {
			_addition(res[index], res[index+stride]);
		}
		__syncthreads();
	}

	//printf("CUDA: 5. (%d,%d), res=%e\n", tg_index, bg_index, res[threadIdx.x].acc.x);
	if (threadIdx.x == 0 && tg_index < NumTarget) _addition(result[tg_index], res[0]);

	//printf("CUDA: 6. (%d,%d), res=%e\n", tg_index, bg_index, res[threadIdx.x].acc.x);
	// res=(%.3e,%.3e,%.3e)\n",
	/*
		 if (threadIdx.x == 0 && tg_index < NumTarget) {
		printf("CUDA: (%d,%d), result=(%.3e,%.3e,%.3e)\n",  
			 	threadIdx.x, tg_index, result[tg_index].acc.x, result[tg_index].acc.y, result[tg_index].acc.z);
			 	//res[0].acc.x, res[0].acc.y, res[0].acc.z);
	}
	*/

}


// _ means inverse
__device__ void kernel(
		const TargetParticle &i,
		const BackgroundParticle &j,
		Result &res,
		Neighbor &neighbor,
		int &bg_index,
		const int &tg_index) {

	float dx  = j.pos.x - i.pos.x;
	float dy  = j.pos.y - i.pos.y;
	float dz  = j.pos.z - i.pos.z;
	float dvx = j.vel.x - i.vel.x;
	float dvy = j.vel.y - i.vel.y;
	float dvz = j.vel.z - i.vel.z;

	float dr2 = dx*dx + dy*dy + dz*dz;

	if (dr2 == 0) return;

	// neighbor
	if(dr2 < i.r2) {
		if (neighbor.NumNeighbor < 100)
			neighbor.NeighborList[neighbor.NumNeighbor++] = bg_index;
		//neighbor.NeighborList[0] = bg_index;
		//neighbor.NumNeighbor += 1;
		//[tg_index*blockDim.x+threadIdx.x]
		//*(neighbor_num+tg_index*blockDim.x+threadIdx.x) += 1;
		//*(neighbor_list+tg_index*blockDim.x*100+threadIdx.x*100+neighbor_num) = bg_index;
	}

	/*
	if (dr2 < ESP2) {
		dr2 +=  ESP2;
	}
	*/

	float drdv      = dx*dvx + dy*dvy + dz*dvz;
	float drdv3_dr2 = 3*drdv/dr2;
	float _dr3      = rsqrtf(dr2)/dr2;
	float m_dr3     = j.mass*_dr3;
	//float mdot_dr3  = j.mdot*_dr3;

	res.acc.x  += m_dr3 * dx;
	res.acc.y  += m_dr3 * dy;
	res.acc.z  += m_dr3 * dz;

	res.adot.x += m_dr3 * (dvx - drdv3_dr2 * dx); //+ mdot_dr3 * dx;
	res.adot.y += m_dr3 * (dvy - drdv3_dr2 * dy); //+ mdot_dr3 * dy;
	res.adot.z += m_dr3 * (dvz - drdv3_dr2 * dz); //+ mdot_dr3 * dz;

}

__device__ void initializeResult(Result &res) {
	res.acc.x  = 0;
	res.acc.y  = 0;
	res.acc.z  = 0;
	res.adot.x = 0;
	res.adot.y = 0;
	res.adot.z = 0;
}


__device__ void _addition(Result &result, const Result res) {
	result.acc.x += res.acc.x;
	result.acc.y += res.acc.y;
	result.acc.z += res.acc.z;

	result.adot.x += res.adot.x;
	result.adot.y += res.adot.y;
	result.adot.z += res.adot.z;
}

__device__ void _copy(Result &result, const Result res) {
	result.acc.x = res.acc.x;
	result.acc.y = res.acc.y;
	result.acc.z = res.acc.z;

	result.adot.x = res.adot.x;
	result.adot.y = res.adot.y;
	result.adot.z = res.adot.z;
}


__global__ void OrganizeNeighbor(const Neighbor do_neighbor_t[], Neighbor d_neighbor_t[],
		const int nn_offset, const int NumTarget) {

	const int index = blockIdx.x*blockDim.x+threadIdx.x + nn_offset*THREAD*BLOCK;
	if (index >= NumTarget)
		return;

	d_neighbor_t[index].clear();
	for (int j=0;j<THREAD;j++) { // for threads
		for (int k=0; k<do_neighbor_t[index*THREAD+j].NumNeighbor; k++) { // for neighbors
			if (k+d_neighbor_t[index].NumNeighbor >= 100) {
				d_neighbor_t[index].NumNeighbor = 100;
				return;
			}
			d_neighbor_t[index].NeighborList[d_neighbor_t[index].NumNeighbor+k]
				= do_neighbor_t[index*THREAD+j].NeighborList[k];
		}
		d_neighbor_t[index].NumNeighbor += do_neighbor_t[index*THREAD+j].NumNeighbor;
	}
	//printf("CUDA: index:%d, do N = %d\n", index, do_neighbor_t[index*THREAD].NumNeighbor);
	//printf("CUDA: index:%d, d  N = %d\n", index, d_neighbor_t[index].NumNeighbor);
	// I can add sorting process in a case where the number of neighbors exceeds MAX_NUM_NEIGHBOR
	// so that closest neighbors can be included consistently.

}





/*************************************************************************
 *	 Communication with HOST
 *************************************************************************/
void _ReceiveFromHost(
		int _NNB,
		double m[],
		double x[][3],
		double v[][3],
		double mdot[],
		int _NumNeighborMax){
	//time_send -= get_wtime();
	nbodymax       = 100000000;
	NNB            = _NNB;
	NumNeighborMax = _NumNeighborMax;
	isend++;
	assert(NNB <= nbodymax);
	cudaError_t cudaStatus;

	//my_allocate(&h_background, &d_background_tmp, new_size(NNB));
	//cudaMemcpyToSymbol(d_background, &d_background_tmp, new_size(NNB)*sizeof(BackgroundParticle));
	if ((first) || (new_size(NNB) > variable_size )) {
		variable_size = new_size(NNB);
		if (!first) {
			my_free(&h_background , &d_background);
			my_free(&h_result     , &d_result);
			my_free(&h_target     , &d_target);
			my_free(&h_neighbor   , &d_neighbor);
			my_free_d(&do_neighbor);
		}
		else {
			first = false;
		}
		my_allocate(&h_background , &d_background, variable_size);
		my_allocate(&h_result     , &d_result    , variable_size);
		my_allocate(&h_target     , &d_target    , variable_size);
		my_allocate(&h_neighbor   , &d_neighbor  , variable_size);
		my_allocate_d(&do_neighbor, variable_size*THREAD);
	}

	//fprintf(stdout, "CUDA: receive starts\n");
	//printf("CUDA: new size of NNB=%d\n", variable_size);
	/*
	cudaStatus = cudaMallocManaged(&background, new_size(NNB)*sizeof(BackgroundParticle)); //NNB*sizeof(BackgroundParticle));

	*/

	for (int j=0; j<NNB; j++) {
		h_background[j].setParticle(m[j], x[j], v[j], mdot[j]);
	}
	toDevice(h_background,d_background,variable_size);

	//fprintf(stdout, "CUDA: receive done\n");
}



void _InitializeDevice(int irank){

	std::cout << "Initializing CUDA ..." << std::endl;
	// Select CUDA device (optional)
	int device = 0; // Choose GPU device 0
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	std::cout << "There are " << deviceCount << " GPUs." << std::endl;
	if (device < 0 || device >= deviceCount) {
		    // Handle invalid device index
	}

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, devid);
	//  char *hostname = getenv("HOSTNAME");
	char hostname[150];
	memset(hostname,0,150);
	gethostname(hostname,150);
	fprintf(stderr, "# GPU initialization - rank: %d; HOST %s; NGPU %d; device: %d %s\n", irank, hostname,numGPU, devid, prop.name);

	cudaSetDevice(device);

	// Initialize CUDA context
	/*
	cudaError_t cudaStatus = cudaFree(0);
	if (cudaStatus != cudaSuccess) {
		std::cerr << "CUDA initialization failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		return;
	}
	*/

	is_open = true;
	// CUDA is now initialized and ready to be used
	std::cout << "CUDA initialized successfully!" << std::endl;

	/*
	if(devinit) return;

	cudaGetDeviceCount(&numGPU);
	assert(numGPU > 0);
	char *gpu_list = getenv("GPU_LIST");
	if(gpu_list)
	{
		numGPU = 0;
		char *p = strtok(gpu_list, " ");
		if (p) {
			devid = atoi(p);
			numGPU++;
		}
		assert(numGPU > 0);
	}else{
		devid=irank%numGPU;
	}
	cudaSetDevice(devid);

#ifdef PROFILE
	//  if(!irank)fprintf(stderr, "***********************\n");
	//  if(!irank)fprintf(stderr, "Initializing NBODY6/GPU library\n");
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, devid);
	//  char *hostname = getenv("HOSTNAME");
	char hostname[150];
	memset(hostname,0,150);
	gethostname(hostname,150);
	fprintf(stderr, "# GPU initialization - rank: %d; HOST %s; NGPU %d; device: %d %s\n", irank, hostname,numGPU, devid, prop.name);
	//  if(!irank)fprintf(stderr, "***********************\n");
#endif
	devinit = true;
	*/
}



void _OpenDevice(const int irank){
	time_send = time_grav = time_nb = time_out = 0.0;
	numInter = 0;
	icall = ini = isend = 0;

	//select GPU========================================//
	_InitializeDevice(irank);

	if(is_open){
		fprintf(stderr, "gpunb: it is already open\n");
		return;
	}
	is_open = true;


#ifdef PROFILE
	//	fprintf(stderr, "RANK: %d ******************\n",irank);
	//	fprintf(stderr, "Opened NBODY6/GPU library\n");
	fprintf(stderr, "# Open GPU regular force - rank: %d\n", irank);
	//fprintf(stderr, "***********************\n");
#endif
}



void _CloseDevice() {
	if(!is_open) {
		fprintf(stderr, "gpunb: it is already close\n");
		return;
	}
	is_open = false;


	cudaError_t error;

	printf("CUDA: ?!! ...\n");
	my_free(h_result    , d_result);
	fprintf(stderr, "result ...\n");
	my_free(h_target    , d_target);
	fprintf(stderr, "target ...\n");
	my_free(h_neighbor  , d_neighbor);
	fprintf(stderr, "neighbor ...\n");
	my_free(h_background, d_background);

	error = cudaGetLastError();
	if (error != cudaSuccess) {
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		// Handle error
	}

#ifdef PROFILE
	fprintf(stderr, "Closed NBODY6/GPU library\n");
	fprintf(stderr, "rank: %d***************\n",devid);
	fprintf(stderr, "time send : %f sec\n", time_send);
	fprintf(stderr, "time grav : %f sec\n", time_grav);
	fprintf(stderr, "time nb   : %f sec\n", time_nb);
	fprintf(stderr, "time out  : %f sec\n", time_out);
	fprintf(stderr, "%f Gflops (gravity part only)\n", 60.e-9 * numInter / time_grav);
	fprintf(stderr, "***********************\n");
#endif
}



void _ProfileDevice(int irank) {
#ifdef PROFILE
	if(icall) {
		fprintf(stderr,"[R.%d-D.%d GPU Reg.F ] Nsend %d  Ngrav %d  <Ni> %d   send(s) %f grav(s) %f  nb(s) %f  out(s) %f  Perf.(Gflops) %f\n",irank,devid,isend,icall,ini/isend,time_send,time_grav,time_nb,time_out,60.e-9*numInter/time_grav);
	}
	time_send = time_grav = time_nb = time_out = 0.0;
	numInter = 0;
	icall = ini = isend= 0;
#else
	return;
#endif
}


#define mexPrintf printf

inline void gpuMemReport(size_t * avail, size_t * total, 
		        const char * title = 0, const size_t * free = 0, const bool sense = true) 
{
	char tstring[32] = { '\0' };
	cudaMemGetInfo(avail, total);  

	if (free) {
		if (title) {
			strncpy(tstring, title, 31);
		}
		mexPrintf("%s Memory avaliable: Free: %zu, Total: %zu, %s: %zu\n",
				tstring, *avail, *total, (sense) ? "Allocated\0" : "Freed\0", 
				(sense) ? (*free - *avail) : (*avail - *free));
	} else {
		mexPrintf("Memory avaliable: Free: %zu, Total: %zu\n", *avail, *total);  
	}
}



extern "C" {
	void InitializeDevice(int *irank){
		_InitializeDevice(*irank);
	}
	void OpenDevice(const int *irank){
		_OpenDevice(*irank);
	}
	void CloseDevice(){
		_CloseDevice();
	}
	void SendToDevice(int *_NNB, double m[], double x[][3], double v[][3], double mdot[], int *_NumNeighborMax) {
		_ReceiveFromHost(*_NNB, m, x, v, mdot, *_NumNeighborMax);
	}
	void ProfileDevice(int *irank){
		_ProfileDevice(*irank);
	}
	void CalculateAccelerationOnDevice(int *NumTarget, double x[][3], double v[][3], double acc[][3], double adot[][3], double mdot[], double radius[], int NumNeighbor[], int **NeighborList) {
		GetAcceleration(*NumTarget, x, v, acc, adot, mdot, radius, NumNeighbor, NeighborList);
	}
}
