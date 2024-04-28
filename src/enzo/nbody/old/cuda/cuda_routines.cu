#include <stdio.h>
#include "cuda_global.h"

#define NTHREAD 64 // 64, 96, 128 or 192
#define NJBLOCK 28 // 8800GTS/512 has 16
#define NIBLOCK 16 // 16 or 32 
#define NIMAX (NTHREAD * NIBLOCK) // 1024

double time_send, time_grav, time_out, time_nb;
long long numInter;
int icall,ini,isend;
int nbodymax;
int devid, numGPU;
bool is_open;
bool devinit;
cudaPointer <BackgroundParticle> background;


void _ReceiveFromHost(
		int _NNB,
		double m[],
		double x[][3],
		double v[][3],
		double mdot[]){
	fprintf(stdout, "CUDA: Send Starting ...\n");
	//time_send -= get_wtime();
	nbodymax = 100000000;
	isend++;
	NNB = _NNB;
	assert(NNB <= nbodymax);
	background.allocate(NNB);

	for(int j=0; j<NNB; j++){
		background[j] = BackgroundParticle(m[j], x[j], v[j],mdot[j]);
	}

	// size_t jpsize = nj * sizeof(Jparticle);
	// cudaMemcpy(jp_dev, jp_host, jpsize, cudaMemcpyHostToDevice);
	fprintf(stdout, "CUDA: send1 ...\n");
	background.toDevice(NNB);
	fprintf(stdout, "CUDA: send2 ...\n");
	//time_send += get_wtime();
}



void _InitializeDevice(int irank){
	//if(devinit) return;

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
	//devinit = true;
}



void _OpenDevice(int irank){
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

	printf("cuda: OpenDevice: %s\n", is_open ? "true" : "false");
	/*
	Background.allocate(nbmax);
	Target.allocate(NIMAX);
	Acc.allocate(NIMAX);
	fobuf.allocate(NIMAX);
	nbpart.allocate(NIMAX);
	nblist.allocate(NB_BUF_SIZE);
	nboff.allocate(NIMAX+1);
	nbodymax = nbmax;
	*/

	//    nblist.reserve(nbmax);
#ifdef PROFILE
	//	fprintf(stderr, "RANK: %d ******************\n",irank);
	//	fprintf(stderr, "Opened NBODY6/GPU library\n");
	fprintf(stderr, "# Open GPU regular force - rank: %d; nbmax: %d\n", irank, nbmax);
	//fprintf(stderr, "***********************\n");
#endif
}



void _CloseDevice(){
	if(!is_open){
		fprintf(stderr, "gpunb: it is already close\n");
		return;
	}
	is_open = false;

	background.free();
	/*
	jpbuf.free();
	ipbuf.free();
	fopart.free();
	fobuf.free();
	nbpart.free();
	nblist.free();
	nboff.free();
	*/

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
		// R: rank; D: GPU device number;
		// Nsend: number of call gpunb_send between two profile check points (adjust time interval);
		// Ngrav: number of call gpunb_ref between two profile check points;
		// <Ni>:  averaged i particle number per regular block;
		// send: j particle sending time;
		// grav: force calculation time;
		// nb:   gather neighbor and force time;
		// out:  output data time;
		// Perf: performance for gpu regular calculation
		fprintf(stderr,"[R.%d-D.%d GPU Reg.F ] Nsend %d  Ngrav %d  <Ni> %d   send(s) %f grav(s) %f  nb(s) %f  out(s) %f  Perf.(Gflops) %f\n",irank,devid,isend,icall,ini/isend,time_send,time_grav,time_nb,time_out,60.e-9*numInter/time_grav);
	}
	time_send = time_grav = time_nb = time_out = 0.0;
	numInter = 0;
	icall = ini = isend= 0;
#else
	return;
#endif
}


extern "C" {
	void InitializeDevice(int *irank){
		_InitializeDevice(*irank);
	}
	void OpenDevice(int *irank){
		_OpenDevice(*irank);
	}
	void CloseDevice(){
		_CloseDevice();
	}
	void SendToDevice(int *_NNB, double m[], double x[][3], double v[][3], double mdot[]) {
		_ReceiveFromHost(*_NNB, m, x, v, mdot);
	}
	void ProfileDevice(int *irank){
		_ProfileDevice(*irank);
	}
}
