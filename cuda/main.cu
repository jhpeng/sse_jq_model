#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <curand_kernel.h>

#include <gsl/gsl_rng.h>

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

/* ---------------------------------------------------- **
** ----------------- allocate memory ------------------ **
** ---------------------------------------------------- */

static int Nstream;
static int Nsite;
static int Nj,Nq;
static double Beta;
static int* d_Length;
static int* d_Noo;
static int* d_Sequence;
static int* d_Linkv;
static int* d_Sigma0;
static int* d_Sigmap;
static int* d_Vfirst;
static int* d_Vlast;
static int* d_Bond2index;
static double* d_Bondst;
static curandState* d_Curng;

void device_create_placeholder(int nstream, int length, int nsite, int nj, int nq, double beta, gsl_rng* rng){
    Nstream = nstream;
    Nsite = nsite;
    Nj = nj;
    Nq = Nq;
    Beta = beta;
    int noo=0;

    cudaMalloc(&d_Length,sizeof(int));
    cudaMalloc(&d_Noo,sizeof(int));
    cudaMalloc(&d_Sequence,length*sizeof(int));
    cudaMalloc(&d_Linkv,8*length*sizeof(int));
    cudaMalloc(&d_Sigma0,nstream*nsite*sizeof(int));
    cudaMalloc(&d_Sigmap,nstream*nsite*sizeof(int));
    cudaMalloc(&d_Vfirst,nstream*nsite*sizeof(int));
    cudaMalloc(&d_Vlast,nstream*nsite*sizeof(int));
    cudaMalloc(&d_Bond2index,4*(nj+nq)*sizeof(int));
    cudaMalloc(&d_Bondst,(nj+nq)*sizeof(double));

    int* sequence = (int*)malloc(length*sizeof(int));
    int* linkv    = (int*)malloc(8*length*sizeof(int));
    int* sigma0   = (int*)malloc(nstream*nsite*sizeof(int));
    int* sigmap   = (int*)malloc(nstream*nsite*sizeof(int));
    int* vfirst   = (int*)malloc(nstream*nsite*sizeof(int));
    int* vlast    = (int*)malloc(nstream*nsite*sizeof(int));
    for(int i=0;i<length;++i) sequence[i]=-1;
    for(int i=0;i<8*length;++i) linkv[i]=-1;
    for(int i=0;i<nsite;++i){
        int spin;
        if(gsl_rng_uniform_pos(rng)<0.5) spin=1;
        else spin=-1;
        for(int j=0;j<nstream;++j){
            sigma0[j*nsite+i] = spin;
            sigmap[j*nsite+i] = spin;
            vfirst[j*nsite+i] = -1;
            vlast[j*nsite+i] = -1;
        }
    }
    

    cudaMemcpy(d_Length,&length,sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_Noo,&noo,sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_Sequence,sequence,length*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_Linkv,linkv,8*length*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_Sigma0,sigma0,nstream*nsite*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_Sigmap,sigmap,nstream*nsite*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_Vfirst,vfirst,nstream*nsite*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_Vlast,vlast,nstream*nsite*sizeof(int),cudaMemcpyHostToDevice);

    d_Curng = get_device_rand(nstream,rng);

    free(sequence);
    free(linkv);
    free(sigma0);
    free(sigmap);
    free(vfirst);
    free(vlast);
}

void device_destroy_placeholder(){
    if(d_Length!=NULL) cudaFree(d_Length);
    if(d_Noo!=NULL) cudaFree(d_Noo);
    if(d_Sequence!=NULL) cudaFree(d_Sequence);
    if(d_Linkv!=NULL) cudaFree(d_Linkv);
    if(d_Sigma0!=NULL) cudaFree(d_Sigma0);
    if(d_Sigmap!=NULL) cudaFree(d_Sigmap);
    if(d_Vfirst!=NULL) cudaFree(d_Vfirst);
    if(d_Vlast!=NULL) cudaFree(d_Vlast);
    if(d_Bond2index!=NULL) cudaFree(d_Bond2index);
    if(d_Bondst!=NULL) cudaFree(d_Bondst);
}

/* ---------------------------------------------------- **
** ------------------ SSE algorithm ------------------- **
** ---------------------------------------------------- */

__device__ void device_propagate_state(int* d_sigmap, int* d_bond2index, int sp, int nsite, int id){
    int i_bond = sp/6;
    int type   = sp%6;

    if(sp==-1 || type==0 || type==2) return;
    else if(type==1 || type==4){
        d_sigmap[id*nsite+d_bond2index[i_bond*4+0]] *= -1;
        d_sigmap[id*nsite+d_bond2index[i_bond*4+1]] *= -1;
    }
    else if(type==3){
        d_sigmap[id*nsite+d_bond2index[i_bond*4+2]] *= -1;
        d_sigmap[id*nsite+d_bond2index[i_bond*4+3]] *= -1;
    }
    else if(type==5){
        d_sigmap[id*nsite+d_bond2index[i_bond*4+0]] *= -1;
        d_sigmap[id*nsite+d_bond2index[i_bond*4+1]] *= -1;
        d_sigmap[id*nsite+d_bond2index[i_bond*4+2]] *= -1;
        d_sigmap[id*nsite+d_bond2index[i_bond*4+3]] *= -1;
    }
    return;
}

__global__ void device_diagonal_update(int nstream, int nsite, int nj, int nq, double beta, int is, int* d_length, int* d_noo, int* d_sequence, int* d_sigma0, int* d_sigmap, int* d_bond2index, double* d_bondst, curandState* d_curng){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id<nstream){
        int length = *d_length;
        int noo = *d_noo;
        int p_max = (int)(length*(id+1)/nstream);
        int p_min = (int)(length*id/nstream);
        int i,p,i_bond,s1,s2,s3,s4;
        double dis;

        for(i=0;i<nsite;++i) d_sigmap[id*nsite+i]=d_sigma0[id*nsite+i];

        if(id==is){
            for(p=p_min;p<p_max;++p){
                if(d_sequence[p]==-1){
                    i_bond = (int)((nj+nq)*curand_uniform_double(&d_curng[id]));
                    s1 = d_sigmap[id*nsite+d_bond2index[i_bond*4+0]];
                    s2 = d_sigmap[id*nsite+d_bond2index[i_bond*4+1]];
                    s3 = d_sigmap[id*nsite+d_bond2index[i_bond*4+2]];
                    s4 = d_sigmap[id*nsite+d_bond2index[i_bond*4+3]];
                    if(i_bond<nj && s1!=s2){
                        dis = curand_uniform_double(&d_curng[id]);
                        if(dis*2*(length-noo)<beta*d_bondst[i_bond]*(nj+nq)){
                            d_sequence[p] = i_bond*6;
                            noo++;
                        }
                    }
                    else if(s1!=s2 && s3!=s4){
                        dis = curand_uniform_double(&d_curng[id]);
                        if(dis*4*(length-noo)<beta*d_bondst[i_bond]*(nj+nq)){
                            d_sequence[p] = i_bond*6+2;
                            noo++;
                        }
                    }
                }
                else if(d_sequence[p]%6==0){
                    i_bond = d_sequence[p]/6;
                    dis = curand_uniform_double(&d_curng[id]);
                    if(beta*(nj+nq)*d_bondst[i_bond]*dis<2*(length-noo+1)){
                        d_sequence[p] = -1;
                        noo--;
                    }
                }
                else if(d_sequence[p]%6==2){
                    i_bond = d_sequence[p]/6;
                    dis = curand_uniform_double(&d_curng[id]);
                    if(beta*(nj+nq)*d_bondst[i_bond]*dis<4*(length-noo+1)){
                        d_sequence[p] = -1;
                        noo--;
                    }
                }
                else device_propagate_state(d_sigmap,d_bond2index,d_sequence[p],nsite,id);
            }
            *d_noo = noo;
            printf("%d\n",noo);
        }
        else{
            for(p=p_min;p<p_max;++p){
                if(d_sequence[p]%6==0){
                    i_bond = (int)((nj+nq)*curand_uniform_double(&d_curng[id]));
                    s1 = d_sigmap[id*nsite+d_bond2index[i_bond*4+0]];
                    s2 = d_sigmap[id*nsite+d_bond2index[i_bond*4+1]];
                    s3 = d_sigmap[id*nsite+d_bond2index[i_bond*4+2]];
                    s4 = d_sigmap[id*nsite+d_bond2index[i_bond*4+3]];
                    if(i_bond<nj && s1!=s2){
                        dis = curand_uniform_double(&d_curng[id]);
                        if(dis*d_bondst[d_sequence[p]/6]<d_bondst[i_bond]){
                            d_sequence[p] = i_bond*6;
                        }
                    }
                    else if(s1!=s2 && s3!=s4){
                        dis = curand_uniform_double(&d_curng[id]);
                        if(dis*d_bondst[d_sequence[p]/6]*2<d_bondst[i_bond]){
                            d_sequence[p] = i_bond*6+2;
                        }
                    }
                }
                else if(d_sequence[p]%6==2){
                    i_bond = (int)((nj+nq)*curand_uniform_double(&d_curng[id]));
                    s1 = d_sigmap[id*nsite+d_bond2index[i_bond*4+0]];
                    s2 = d_sigmap[id*nsite+d_bond2index[i_bond*4+1]];
                    s3 = d_sigmap[id*nsite+d_bond2index[i_bond*4+2]];
                    s4 = d_sigmap[id*nsite+d_bond2index[i_bond*4+3]];
                    if(i_bond<nj && s1!=s2){
                        dis = curand_uniform_double(&d_curng[id]);
                        if(dis*d_bondst[d_sequence[p]/6]<d_bondst[i_bond]*2){
                            d_sequence[p] = i_bond*6;
                        }
                    }
                    else if(s1!=s2 && s3!=s4){
                        dis = curand_uniform_double(&d_curng[id]);
                        if(dis*d_bondst[d_sequence[p]/6]<d_bondst[i_bond]){
                            d_sequence[p] = i_bond*6+2;
                        }
                    }
                }
                else device_propagate_state(d_sigmap,d_bond2index,d_sequence[p],nsite,id);
            }
        }
    }
}

/* ---------------------------------------------------- **
** ------------------ set the model ------------------- **
** ---------------------------------------------------- */


void device_set_lattice_jq_isotropy_2d(int nstream, int length, int nx, int ny, double jbond, double beta, gsl_rng* rng){
    int i,j,t,q;
    Beta = beta;
    Nsite = nx*ny;
    Nj = 2*Nsite;
    Nq = 2*Nsite;

    int* bond2index = (int*)malloc((Nj+Nq)*4*sizeof(int));
    double* bondst = (double*)malloc((Nj+Nq)*sizeof(double));

    for(int i_bond=0;i_bond<(Nj+Nq);++i_bond){
        t = i_bond%Nsite;
        q = i_bond/Nsite;
        i = t%nx;
        j = t/nx;

        if(q==0){
            bond2index[i_bond*4+0] = i+nx*j;
            bond2index[i_bond*4+1] = ((i+1)%nx)+nx*j;
            bond2index[i_bond*4+2] = -1;
            bond2index[i_bond*4+3] = -1;
            bondst[i_bond] = jbond;
        }
        else if(q==1){
            bond2index[i_bond*4+0] = i+nx*j;
            bond2index[i_bond*4+1] = i+nx*((j+1)%ny);
            bond2index[i_bond*4+2] = -1;
            bond2index[i_bond*4+3] = -1;
            bondst[i_bond] = jbond;
        }
        else if(q==2){
            bond2index[i_bond*4+0] = i+nx*j;
            bond2index[i_bond*4+1] = ((i+1)%nx)+nx*j;
            bond2index[i_bond*4+2] = i+nx*((j+1)%ny);
            bond2index[i_bond*4+3] = ((i+1)%nx)+nx*((j+1)%ny);
            bondst[i_bond] = 1.0;
        }
        else if(q==3){
            bond2index[i_bond*4+0] = i+nx*j;
            bond2index[i_bond*4+1] = i+nx*((j+1)%ny);
            bond2index[i_bond*4+2] = ((i+1)%nx)+nx*j;
            bond2index[i_bond*4+3] = ((i+1)%nx)+nx*((j+1)%ny);
            bondst[i_bond] = 1.0;
        }
    }

    device_create_placeholder(nstream,length,Nsite,Nj,Nq,Beta,rng);
    cudaMemcpy(d_Bond2index,bond2index,4*(Nj+Nq)*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_Bondst,bondst,(Nj+Nq)*sizeof(double),cudaMemcpyHostToDevice);
}

int main()
{
    int nstream=1024;
    int nblock =128;
    int nthread=(int)(nstream/nblock);
    int length=1000000;
    int nx=256;
    int ny=256;
    double beta=20.0;
    double jbond=1.0;
    int seed=219871;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    device_set_lattice_jq_isotropy_2d(nstream,length,nx,ny,jbond,beta,rng);
    for(int i=0;i<1000;++i){
        int is = (int)(gsl_rng_uniform_pos(rng)*Nstream);
        device_diagonal_update<<<nblock,nthread>>>(nstream,Nsite,Nj,Nq,Beta,is,d_Length,d_Noo,d_Sequence,d_Sigma0,d_Sigmap,d_Bond2index,d_Bondst,d_Curng);
        cudaDeviceSynchronize();
    }

    gsl_rng_free(rng);
    device_destroy_placeholder();

    printf("%s\n",cudaGetErrorName(cudaGetLastError()));
}
