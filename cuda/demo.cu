#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <cuda_runtime_api.h>

__global__ void test_random(double* rand_d, int length, curandState* s){
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if(id<length){
        curand_init(0,id,0,&s[id]);
        rand_d[id] = curand_uniform_double(&s[id]);
    }
}

int main(void)
{
    cudaSetDeviceFlags(cudaDeviceMapHost);
    int length=100000;

    double* data_d;
    double* data_h;
    cudaMalloc(&data_d, length*sizeof(double));
    data_h =  (double*)malloc(length*sizeof(double));

	curandState *s;
    cudaMalloc(&s,length*sizeof(curandState));

    int const n_thread = length/200;
    int const n_block = 200;
    test_random<<<n_block,n_thread>>>(data_d,length,s);
    cudaError_t err;
    err = cudaThreadSynchronize();

    cudaMemcpy(data_h,data_d,length*sizeof(double),cudaMemcpyDeviceToHost);

    for(int i=0;i<length;++i) printf("%.5f \n",data_h[i]);

    cudaFree(data_d);
    free(data_h);
    cudaFree(s);
}
