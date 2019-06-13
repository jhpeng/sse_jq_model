#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <curand_kernel.h>

#include <gsl/gsl_rng.h>

typedef struct placeholder{
    int nstream;
    int length;
    int noo;
    int* sequence;
    int* linkv;
    int nsite;
    int **sigma0;
    int **sigmap;
    int **vfirst;
    int **vlast;
    int nj;
    int nq;
    int* bond2index;
    double* bondst;
    double beta;
    curandState* cu_rng;
} placeholder;

placeholder* device_create_placeholder(int nstream, int length, int nsite, int nj, int nq, double beta){
    
}

/* ---------------------------------------------------- **
** ------------------- set curand --------------------- **
** ---------------------------------------------------- */

__global__ void set_curand(curandState* s, unsigned long long* seed, int size){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id<size){
        curand_init(seed[id],0,0,&s[id]);
    }
}

curandState* get_device_rand(int size, gsl_rng* rng){
    curandState* cu_rng;
    cudaMalloc(&cu_rng,size*sizeof(curandState));
    unsigned long long* h_seed;
    unsigned long long* d_seed;
    h_seed = (unsigned long long*)malloc(size*sizeof(unsigned long long));
    cudaMalloc(&d_seed,size*sizeof(unsigned long long));

    for(int i=0;i<size;++i){
        h_seed[i] = (unsigned long long)(gsl_rng_uniform_pos(rng)*ULLONG_MAX);
    }
    cudaMemcpy(d_seed,h_seed,size*sizeof(unsigned long long),cudaMemcpyHostToDevice);
    int nblock = 20;
    int nthread = size/20+1;
    set_curand<<<nblock,nthread>>>(cu_rng,d_seed,size);

    cudaDeviceSynchronize();
    free(h_seed);
    cudaFree(d_seed);

    return cu_rng;
}

int main()
{
    int size=1024;
    int seed=219871;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    curandState* cu_rng = get_device_rand(size,rng);

    gsl_rng_free(rng);
    cudaFree(cu_rng);

    printf("%s\n",cudaGetErrorName(cudaGetLastError()));
}
