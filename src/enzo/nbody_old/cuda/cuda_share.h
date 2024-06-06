#pragma once

#define NAN_CHECK(val) assert((val) == (val));

#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#else
static double get_wtime(){
	return 0.0;
}
#endif

struct Jparticle{
  float3 pos;
  float  mass;
  float3 vel;
  float  currentt;
  float3 acc;
  float3 accdot;
  
  Jparticle(int) {}
  Jparticle(double mj, double xj[3], double vj[3]){
    mass  = mj;
    pos.x = xj[0];
    pos.y = xj[1];
    pos.z = xj[2];
    vel.x = vj[0];
    vel.y = vj[1];
    vel.z = vj[2];


    NAN_CHECK(xj[0]);
    NAN_CHECK(xj[1]);
    NAN_CHECK(xj[2]);
    NAN_CHECK(mj);
    NAN_CHECK(vj[0]);
    NAN_CHECK(vj[1]);
    NAN_CHECK(vj[2]);
  }
  Jparticle(double mj, double ctj, double xj[3], double vj[3], double aj[3], double adotj[3]){
    pos.x = xj[0];
    pos.y = xj[1];
    pos.z = xj[2];
    mass  = mj;
    vel.x = vj[0];
    vel.y = vj[1];
    vel.z = vj[2];
    currentt = ctj;
    acc.x = aj[0];
    acc.y = aj[1];
    acc.z = aj[2];
    accdot.x = adotj[0];
    accdot.y = adotj[1];
    accdot.z = adotj[2];

    NAN_CHECK(xj[0]);
    NAN_CHECK(xj[1]);
    NAN_CHECK(xj[2]);
    NAN_CHECK(mj);
    NAN_CHECK(vj[0]);
    NAN_CHECK(vj[1]);
    NAN_CHECK(vj[2]);
    NAN_CHECK(ctj);
    NAN_CHECK(aj[0]);
    NAN_CHECK(aj[1]);
    NAN_CHECK(aj[2]);
    NAN_CHECK(adotj[0]);
    NAN_CHECK(adotj[1]);
    NAN_CHECK(adotj[2]);
  }
  __device__ Jparticle() {}
};

static float2 float2_split(double x){
	const int shift = 20;
	float2 ret;
	x *= (1<<shift);
	double xi = (int)x;
	double xf = x - xi;
	ret.x = xi * (1./(1<<shift));
	ret.y = xf * (1./(1<<shift));
	return ret;
}

static __device__ float2 float2_accum(float2 acc, float x){
  float tmp = acc.x + x;
  acc.y -= (tmp - acc.x) - x;
  acc.x = tmp;
  return acc;
}

static  __device__ float2 float2_regularize(float2 acc){
  float tmp = acc.x + acc.y;
  acc.y = acc.y -(tmp - acc.x);
  acc.x = tmp;
  return acc;
}

static __device__ float2 float2_add(float2 a, float2 b){
  float tmp = a.x + b.x;
  a.y -= (tmp - a.x) - b.x - b.y;
  a.x = tmp;
  // a.x = a.x + b.x;
  // a.y = a.y + b.y;
  return a;
}

extern cudaPointer <Jparticle> jpbuf;
