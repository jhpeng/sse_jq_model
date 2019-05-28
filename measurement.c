/* algorithm/measurement.c
*/

#include <stdio.h>
#include <math.h>

#include "low_level.h"
#include "mid_level.h"
#include "lattice.h"
#include "monte_carlo.h"

static void observable_calc_mean(double* mean, const observable* obs){
    int i_obs,i_sample;
    int Nobs = observable_get_Nobs(obs);
    int Nsample = observable_get_Nsample(obs);
    
    for(i_obs=0;i_obs<Nobs;++i_obs){
        mean[i_obs]=0;
        for(i_sample=0;i_sample<Nsample;++i_sample){
            if(observable_get_init(obs,i_obs,i_sample)==0){
                printf("observable_calc_mean : data missing!\n");
                exit(-1);
            }
            mean[i_obs]+=observable_get_data(obs,i_obs,i_sample);
        }
        mean[i_obs]=mean[i_obs]/Nsample;
    }
}

void observable_result_stdout(const observable* obs, double beta){
    int Nobs = observable_get_Nobs(obs);
    int Nsample = observable_get_Nsample(obs);
    double mean[Nobs];

    observable_calc_mean(mean,obs);
    printf("Nsample : %d\n",Nsample);
    printf("beta    : %.5e\n",beta);
    for(int i=0;i<Nobs;++i){
        printf("%.5e ",mean[i]);
    }
    printf("\n---------------------------------------\n");
}

void observable_result_fileout(const observable* obs, double beta, const char* filename){
    int Nobs = observable_get_Nobs(obs);
    double mean[Nobs];

    observable_calc_mean(mean,obs);

    FILE* f = fopen(filename,"a");
    fprintf(f,"%e ",beta);
    for(int i=0;i<Nobs;++i){
        fprintf(f,"%e ",mean[i]);
    }
    fprintf(f,"\n");
    fclose(f);
}

static int_array* antiferro_ms_2d=NULL;

void observable_measure_ms_2d(observable* obs, lattice_struct* las, const operator_sequence* ops, int i_sample){
    int Nobs = observable_get_Nobs(obs);
    int Nsample = observable_get_Nsample(obs);
    if(Nobs!=3){
        printf("observable_measure_ms_2d : Nobs must be 3\n");
        exit(-1);
    }
    else if(i_sample<0){
        printf("observable_measure_ms_2d : i_sample must be larger or equal to zero!\n");
        exit(-1);
    }
    int Nsite = lattice_struct_get_Nsite(las);
    int i,p,sp,L = operator_sequence_get_length(ops);
    if(antiferro_ms_2d==NULL){
        antiferro_ms_2d = int_array_alloc(Nsite);
        int_array_set_all(antiferro_ms_2d,0);
        int_array_set(antiferro_ms_2d,0,1);
        int s1,s2,index[4];
        for(int i_bond=0;i_bond<(2*Nsite);++i_bond){
            lattice_struct_get_bond2index(index,las,i_bond);
            s1 = int_array_get(antiferro_ms_2d,index[0]);
            s2 = int_array_get(antiferro_ms_2d,index[1]);
            if(s1!=0) s2 = -1*s1;
            else s1 = -1*s2;
            int_array_set(antiferro_ms_2d,index[0],s1);
            int_array_set(antiferro_ms_2d,index[1],s2);
        }
    }

    lattice_struct_sync_sigmap(las);
    double m,m1=0,m2=0,m4=0;
    for(p=0;p<L;++p){
        m=0;
        for(i=0;i<Nsite;++i)m+=lattice_struct_get_sigmap(las,i)*int_array_get(antiferro_ms_2d,i);
        sp = operator_sequence_get_sequence(ops,p);
        lattice_struct_propagate_state(las,sp);

        m1 += fabs(m);
        m2 += m*m;
        m4 += m*m*m*m;
    }


    m1 = m1/L/Nsite*0.5;
    m2 = m2/L/Nsite/Nsite*0.25;
    m4 = m4/L/Nsite/Nsite/Nsite/Nsite*0.0625;
    observable_set_data(obs,0,i_sample%Nsample,m1);
    observable_set_data(obs,1,i_sample%Nsample,m2);
    observable_set_data(obs,2,i_sample%Nsample,m4);
}

void observable_measure_mz_2d(observable* obs, const lattice_struct* las, int i_sample){
    int Nobs = observable_get_Nobs(obs);
    int Nsample = observable_get_Nsample(obs);
    if(Nobs!=3){
        printf("observable_measure_mz_2d : Nobs must be 3\n");
        exit(-1);
    }
    else if(i_sample<0){
        printf("observable_measure_mz_2d : i_sample must be larger or equal to zero!\n");
        exit(-1);
    }

    int i;
    int Nsite = lattice_struct_get_Nsite(las);
    double m=0,m1=0,m2=0,m4=0;
    for(i=0;i<Nsite;++i){
        m += lattice_struct_get_sigma0(las,i);
    }
    m1 = m*0.5/Nsite;
    m2 = m*m*0.25/Nsite/Nsite;
    m4 = m*m*m*m*0.0625/Nsite/Nsite/Nsite/Nsite;
    observable_set_data(obs,0,i_sample%Nsample,m1);
    observable_set_data(obs,1,i_sample%Nsample,m2);
    observable_set_data(obs,2,i_sample%Nsample,m4);
}

void observable_measure_uniform_sus_2d(observable* obs, const lattice_struct* las, double beta, int i_sample){
    int Nobs = observable_get_Nobs(obs);
    int Nsample = observable_get_Nsample(obs);
    if(Nobs!=1){
        printf("observable_measure_mz_2d : Nobs must be 1\n");
        exit(-1);
    }
    else if(i_sample<0){
        printf("observable_measure_mz_2d : i_sample must be larger or equal to zero!\n");
        exit(-1);
    }

    int i;
    int Nsite = lattice_struct_get_Nsite(las);
    double m=0,m2=0;
    for(i=0;i<Nsite;++i){
        m += lattice_struct_get_sigma0(las,i);
    }
    m2 = m*m*0.25/Nsite*beta;

    observable_set_data(obs,0,i_sample%Nsample,m2);
}

#if 0
int main(int argc, char **argv)
{
    int Nx=16;
    int Ny=16;
    double J=0.06;
    double beta=64;
    int L = 1000;
    int Nobs=3,Nsample=2000,Nblock=5,Nther=1000;
    int seed=39829;

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
            observable_measure_ms_2d(obs,las,ops,i);
        }
        observable_result_stdout(obs,beta);
        observable_result_fileout(obs,beta,"test.txt");
        observable_init(obs);
    }

    operator_sequence_free(ops);
    lattice_struct_free(las);
    link_vertex_free(lv);
    observable_free(obs);
    
    return 0;
}
#endif
