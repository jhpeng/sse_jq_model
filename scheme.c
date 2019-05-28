/* application/scheme.c
*/
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <gsl/gsl_rng.h>

#include "mid_level.h"
#include "lattice.h"
#include "monte_carlo.h"
#include "measurement.h"

void normal_scheme_isotropy_square_lattice(int Nx, int Ny, double J, double beta, int Nsample, int Nblock, int Nther, int seed, const char* filename){
    int Nobs=3,L=1000;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    lattice_struct* las = lattice_struct_create_model_isotropy_2d(Nx,Ny,J,rng);
    operator_sequence* ops = operator_sequence_alloc(L);
    operator_sequence_init(ops);
    int Nsite = lattice_struct_get_Nsite(las);
    link_vertex* lv = link_vertex_alloc(L,Nsite);
    observable* obs = observable_alloc(Nobs,Nsample);
    observable_init(obs);

    monte_carlo_thermalization(las,&ops,&lv,beta,rng,Nther);

    for(int j=0;j<Nblock;++j){
        for(int i=0;i<Nsample;++i){
            monte_carlo_single_sweep(las,ops,lv,beta,rng);
            /*observable_measure_ms_2d(obs,las,ops,i);*/
            observable_measure_mz_2d(obs,las,i);
            /*observable_measure_uniform_sus_2d(obs,las,beta,i);*/
        }
        observable_result_fileout(obs,beta,filename);
        observable_init(obs);
    }

    operator_sequence_free(ops);
    lattice_struct_free(las);
    link_vertex_free(lv);
    observable_free(obs);
}

static int help=0;
static int mode=0;
static int Nx=8;
static int Ny=8;
static double J=1.0;
static double beta=4.0;
static int Nblock=50;
static int Nsample=2000;
static int Nther=2000;
static int seed=1;
static char filename[128]="test.txt";

int main(int argc, char **argv)
{
    int c;
    while((c=getopt(argc,argv,"hm:x:y:j:b:n:k:t:s:f:"))!=-1){
        switch(c){
            case 'h':
                help=1;
                printf("usage:\n");
                printf("\t-h help\n");
                printf("\t-x <length of x> default 8\n");
                printf("\t-y <length of y> default 8\n");
                printf("\t-m <mode> default 0\n");
                printf("\t-j <J/Q ratio> default 1.0\n");
                printf("\t-b <beta> default 4.0\n");
                printf("\t-n <Nsample> default 2000\n");
                printf("\t-k <Nblock> default 50\n");
                printf("\t-t <Nther> default 2000\n");
                printf("\t-s <seed of random number generator> default 1\n");
                printf("\t-f <the file name of output data> default \"test.txt\"\n");
                break;
            case 'x':
                Nx=atoi(optarg);
                break;
            case 'y':
                Ny=atoi(optarg);
                break;
            case 'm':
                mode=atoi(optarg);
                break;
            case 'j':
                J=atof(optarg);
                break;
            case 'b':
                beta=atof(optarg);
                break;
            case 'n':
                Nsample=atoi(optarg);
                break;
            case 'k':
                Nblock=atoi(optarg);
                break;
            case 't':
                Nther=atoi(optarg);
                break;
            case 's':
                seed=atoi(optarg);
                break;
            case 'f':
                strcpy(filename,optarg);
                break;

        }
    }

    if(help) return 0;

    if(mode==0) normal_scheme_isotropy_square_lattice(Nx,Ny,J,beta,Nsample,Nblock,Nther,seed,filename);
}
