#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <gsl/gsl_rng.h>

/*
**  define global variables here
*/

/* The flip rule is the instruction for flipping operator.
**   | 0 | 1 | 2 | 3 | 4 | 5
------------------------------
** a | 1 | 0 | 4 | 5 | 2 | 3
** b | 0 | 1 | 3 | 2 | 5 | 4
** c | 1 | 0 | 4 | 5 | 2 | 3
** d | 0 | 1 | 3 | 2 | 5 | 4
**
**
**  |    |        |    |
** c|   c|       d|   d| 
** --------      --------
** |\\\\\\|------|\\\\\\|
** --------      --------
**  |    |        |    |
** a|   a|       b|   b|
**
*/
static int flip_rule[24]={
    1,0,4,5,2,3,
    0,1,3,2,5,4,
    1,0,4,5,2,3,
    0,1,3,2,5,4};

/*
** L   : total length of the operator sequence
** Noo : number of non-identity operator in the operator sequence
** Sequence : pointer to the operator sequence
** Linkv    : pointer to the link-vertex list (8L)
*/
int L,Noo;
int* Sequence;
int* Linkv;

/*
** Nsite : the total number of the spin
** Sigma0 : the pointer to the current state
** Sigmap : the pointer to the propagate state
** Vfirst : the pointer to the first list
** Vlast  : the pointer to the last list
*/
int Nsite;
int* Sigma0;
int* Sigmap;
int* Vfirst;
int* Vlast;

/*
** Nj : the total number of the J-bonds
** Nq : the totla number of the Q-bond-pairs
** Bond2index : the pointer which mapping the bond on four spin site
** Bondst : the pointer tp the bond strenght
*/
int Nx,Ny,Nz;
int Nj,Nq;
int* Bond2index;
double* Bondst;

/*
** Nobs    : the total number of the observables
** Nsample : the total number of the Monte Carlo sample
** Data : the pointer to the collected data
*/
int Nobs,Nsample;
double* Data;
char Filename[128]="test.txt";

/*
** Nblock  : the total number of the block data
** Beta : inverse temperature
** rng  : gsl_rng
*/
int Nblock,Nther,Seed;
double Beta,Jbond,Qbond;
gsl_rng* rng;


/* -------------------------------------------------- **
** ---------------- SSE algorithm ------------------- **
** -------------------------------------------------- */

void propagate_state(int sp){
    if(sp==-1) return;

    int type   = sp%6;
    if(type==0 || type==2) return;
    
    int i_bond = sp/6;
    if(type==1 || type==4){
        Sigmap[Bond2index[i_bond*4+0]] *= -1;
        Sigmap[Bond2index[i_bond*4+1]] *= -1;
    }
    else if(type==3){
        Sigmap[Bond2index[i_bond*4+2]] *= -1;
        Sigmap[Bond2index[i_bond*4+3]] *= -1;
    }
    else if(type==5){
        Sigmap[Bond2index[i_bond*4+0]] *= -1;
        Sigmap[Bond2index[i_bond*4+1]] *= -1;
        Sigmap[Bond2index[i_bond*4+2]] *= -1;
        Sigmap[Bond2index[i_bond*4+3]] *= -1;
    }
}

void diagonal_update(){
    int i_bond,s1,s2,s3,s4;
    double dis;

    for(int i=0;i<Nsite;++i){
        Sigmap[i] = Sigma0[i];
    }

    for(int p=0;p<L;++p){
        if(Sequence[p]==-1){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*(Nj+Nq));
            s1 = Sigmap[Bond2index[i_bond*4+0]];
            s2 = Sigmap[Bond2index[i_bond*4+1]];
            s3 = Sigmap[Bond2index[i_bond*4+2]];
            s4 = Sigmap[Bond2index[i_bond*4+3]];
            if(i_bond<Nj && s1!=s2){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*2*(L-Noo)<Beta*Bondst[i_bond]*(Nj+Nq)){
                    Sequence[p] = i_bond*6;
                    Noo++;
                }
            }
            else if(s1!=s2 && s3!=s4){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*4*(L-Noo)<Beta*Bondst[i_bond]*(Nj+Nq)){
                    Sequence[p] = i_bond*6+2;
                    Noo++;
                }
            }
        }
        else if(Sequence[p]%6==0){
            i_bond = Sequence[p]/6;
            dis = gsl_rng_uniform_pos(rng);
            if(Beta*(Nj+Nq)*Bondst[i_bond]*dis<2*(L-Noo+1)){
                Sequence[p]=-1;
                Noo--;
            }
        }
        else if(Sequence[p]%6==2){
            i_bond = Sequence[p]/6;
            dis = gsl_rng_uniform_pos(rng);
            if(Beta*(Nj+Nq)*Bondst[i_bond]*dis<4*(L-Noo+1)){
                Sequence[p]=-1;
                Noo--;
            }
        }
        else propagate_state(Sequence[p]);
    }
}

void diagonal_update_exchange(){
    int i_bond,s1,s2,s3,s4;
    double dis;

    for(int i=0;i<Nsite;++i){
        Sigmap[i] = Sigma0[i];
    }

    for(int p=0;p<L;++p){
        if(Sequence[p]%6==0){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*(Nj+Nq));
            s1 = Sigmap[Bond2index[i_bond*4+0]];
            s2 = Sigmap[Bond2index[i_bond*4+1]];
            s3 = Sigmap[Bond2index[i_bond*4+2]];
            s4 = Sigmap[Bond2index[i_bond*4+3]];
            if(i_bond<Nj && s1!=s2){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*Bondst[Sequence[p]/6]<Bondst[i_bond]){
                    Sequence[p] = i_bond*6;
                }
            }
            else if(s1!=s2 && s3!=s4){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*Bondst[Sequence[p]/6]*2<Bondst[i_bond]){
                    Sequence[p] = i_bond*6+2;
                }
            }
        }
        else if(Sequence[p]%6==2){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*(Nj+Nq));
            s1 = Sigmap[Bond2index[i_bond*4+0]];
            s2 = Sigmap[Bond2index[i_bond*4+1]];
            s3 = Sigmap[Bond2index[i_bond*4+2]];
            s4 = Sigmap[Bond2index[i_bond*4+3]];
            if(i_bond<Nj && s1!=s2){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*Bondst[Sequence[p]/6]<Bondst[i_bond]*2){
                    Sequence[p] = i_bond*6;
                }
            }
            else if(s1!=s2 && s3!=s4){
                dis = gsl_rng_uniform_pos(rng);
                if(dis*Bondst[Sequence[p]/6]<Bondst[i_bond]){
                    Sequence[p] = i_bond*6+2;
                }
            }
        }
        else propagate_state(Sequence[p]);
    }
}

void construct_link_vertex_list(){
    for(int i=0;i<(8*L);++i) Linkv[i]=-1;
    for(int i=0;i<Nsite;++i){
        Vfirst[i]=-1;
        Vlast[i] =-1;
    }

    int i_bond,index,nu0,nu1;
    for(int p=0;p<L;++p){
        if(Sequence[p]!=-1){
            i_bond = Sequence[p]/6;
            for(int i=0;i<4;++i){
                index = Bond2index[i_bond*4+i];
                if(index!=-1){
                    nu0 = 8*p+i;
                    nu1 = Vlast[index];
                    if(nu1!=-1){
                        Linkv[nu0] = nu1;
                        Linkv[nu1] = nu0;
                    }
                    else Vfirst[index] = nu0;

                    Vlast[index] = nu0+4;
                }
            }
        }
    }

    for(int i=0;i<Nsite;++i){
        if(Vlast[i]!=-1){
            Linkv[Vlast[i]] = Vfirst[i];
            Linkv[Vfirst[i]] = Vlast[i];
        }
    }
}

void loop_update(){
   int nu0,nup,nun,flip=-1;

    for(nu0=0;nu0<(L*8);nu0+=2){
        if(Linkv[nu0]>=0){
            nun = nu0;
            if(gsl_rng_uniform_pos(rng)<0.5) flip=-1;
            else flip=-2;
            while(Linkv[nun]>=0){
                Linkv[nun] = flip;
                nup = nun^1;
                nun = Linkv[nup];
                Linkv[nup]=flip;
            }
        }
    } 
}

void flip_bit_operator(){
    int nu,i_flip_rule,type,i_bond,index;
    
    for(nu=0;nu<(8*L);nu+=2){
        if(Linkv[nu]==-2){
            type = Sequence[nu/8]%6;
            i_bond = Sequence[nu/8]/6;
            i_flip_rule = ((nu%8)/2)*6+type;
            Sequence[nu/8] = i_bond*6+flip_rule[i_flip_rule];
        }
    }

    for(index=0;index<Nsite;++index){
        nu = Vlast[index];
        if(nu==-1){
            if(gsl_rng_uniform_pos(rng)<0.5){
                Sigma0[index] *= -1;
            }
        }
        else if(nu>=0){
            if(Linkv[nu]==-2){
                Sigma0[index] *= -1;
            }
        }
    }
}



/* --------------------------------------------------------- **
** ----------------------- Estimator ----------------------- **
** --------------------------------------------------------- */

void calc_mean(double* mean){
    for(int i_obs=0;i_obs<Nobs;++i_obs){
        mean[i_obs]=0;
        for(int i=0;i<Nsample;++i) mean[i_obs]+=Data[Nobs*i+i_obs];
        mean[i_obs] = mean[i_obs]/Nsample;
    }
}

void estimator_stdout(){
    double mean[Nobs];
    calc_mean(mean);

    printf("Nsample : %d\n",Nsample);
    printf("beta    : %.5e\n",Beta);
    for(int i=0;i<Nobs;++i) printf("%.5e ",mean[i]);
    printf("\n---------------------------------------\n");
}

void estimator_fileout(char* filename){
    double mean[Nobs];
    calc_mean(mean);

    FILE* ofile = fopen(filename,"a");
    fprintf(ofile,"%.5e ",Beta);
    for(int i=0;i<Nobs;++i) fprintf(ofile,"%.5e ",mean[i]);
    fprintf(ofile,"\n");
    fclose(ofile);
}

void measure_mz(int i_obs, int i_sample){
    double m1,m2,mz=0;
    for(int i=0;i<Nsite;++i) mz+=Sigma0[i];

    m1 = mz*0.5;
    m2 = mz*mz*0.25;
    Data[Nobs*i_sample+i_obs+0] = m1;
    Data[Nobs*i_sample+i_obs+1] = m2*4;
    Data[Nobs*i_sample+i_obs+2] = m2*Beta/Nsite;
    Data[Nobs*i_sample+i_obs+3] = m2*Beta*Beta/Nsite;
}

/* --------------------------------------------------------- **
** ---------------- Setting the model ---------------------- **
** --------------------------------------------------------- */


void set_random_number(int seed){
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
}

void set_lattice_jq_isotropy_2d(int nx, int ny, double jbond){
    int i,j,t,q;
    Nsite = nx*ny;
    Nj = 2*Nsite;
    Nq = 2*Nsite;

    Sigma0 = (int*)malloc(Nsite*sizeof(int));
    Sigmap = (int*)malloc(Nsite*sizeof(int));
    Vfirst = (int*)malloc(Nsite*sizeof(int));
    Vlast  = (int*)malloc(Nsite*sizeof(int));

    Bond2index = (int*)malloc((Nj+Nq)*4*sizeof(int));
    Bondst = (double*)malloc((Nj+Nq)*sizeof(double));

    for(int i_bond=0;i_bond<(Nj+Nq);++i_bond){
        t = i_bond%Nsite;
        q = i_bond/Nsite;
        i = t%nx;
        j = t/nx;

        if(q==0){
            Bond2index[i_bond*4+0] = i+nx*j;
            Bond2index[i_bond*4+1] = ((i+1)%nx)+nx*j;
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;
            Bondst[i_bond] = jbond;
        }
        else if(q==1){
            Bond2index[i_bond*4+0] = i+nx*j;
            Bond2index[i_bond*4+1] = i+nx*((j+1)%ny);
            Bond2index[i_bond*4+2] = -1;
            Bond2index[i_bond*4+3] = -1;
            Bondst[i_bond] = jbond;
        }
        else if(q==2){
            Bond2index[i_bond*4+0] = i+nx*j;
            Bond2index[i_bond*4+1] = ((i+1)%nx)+nx*j;
            Bond2index[i_bond*4+2] = i+nx*((j+1)%ny);
            Bond2index[i_bond*4+3] = ((i+1)%nx)+nx*((j+1)%ny);
            Bondst[i_bond] = 1.0;
        }
        else if(q==3){
            Bond2index[i_bond*4+0] = i+nx*j;
            Bond2index[i_bond*4+1] = i+nx*((j+1)%ny);
            Bond2index[i_bond*4+2] = ((i+1)%nx)+nx*j;
            Bond2index[i_bond*4+3] = ((i+1)%nx)+nx*((j+1)%ny);
            Bondst[i_bond] = 1.0;
        }
    }

    for(i=0;i<Nsite;++i){
        if(gsl_rng_uniform_pos(rng)<0.5) Sigma0[i]=1;
        else Sigma0[i]=-1;
    }
}

void set_sequence_length(int length){
    if(Sequence==NULL){
        Sequence = (int*)malloc(length*sizeof(int));
        Noo=0;
        for(int p=0;p<length;++p){
            Sequence[p]=-1;
        }
    }
    else{
        int* temp = (int*)malloc(length*sizeof(int));
        for(int p=0;p<length;++p){
            if(p<L) temp[p] = Sequence[p];
            else temp[p]=-1;
        }
        free(Sequence);
        Sequence = temp;
    }

    if(Linkv!=NULL) free(Linkv);
    Linkv = (int*)malloc(8*length*sizeof(int));

    L = length;
}

void set_estimator(int n_obs, int n_sample, int n_block){
    Data = (double*)malloc(n_obs*n_sample*sizeof(double));
    Nobs = n_obs;
    Nsample = n_sample;
    Nblock = n_block;
}



/* ----------------------------------------------- **
** ------------------ getopt --------------------- **
** ----------------------------------------------- */ 

int Help,Exc=0;
void set_opt(int argc, char **argv)
{
    int c;
    while((c=getopt(argc,argv,"hx:y:j:b:n:k:t:s:f:e:"))!=-1){
        switch(c){
            case 'h':
                Help=1;
                printf("usage:\n");
                printf("\t-h help\n");
                printf("\t-x <length of x> default 8\n");
                printf("\t-y <length of y> default 8\n");
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
            case 'j':
                Jbond=atof(optarg);
                break;
            case 'b':
                Beta=atof(optarg);
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
                Seed=atoi(optarg);
                break;
            case 'f':
                strcpy(Filename,optarg);
                break;
            case 'm':
                Exc=atoi(optarg);
                break;
        }
    }
}


/* ----------------------------------------------- **
** ------------------ main ----------------------- **
** ----------------------------------------------- */ 

int main(int argc, char** argv){
    int length=1000;
    int n_obs=4;
    double buffer=1.3;

    Beta = 4096;
    Seed = 9237912;
    Nx = 48;
    Ny = 48;
    Jbond = 0.04;
    Nther = 20000;
    Nsample = 2000;
    Nblock = 50;

    set_opt(argc,argv);
    if(Help) return 0;

    set_random_number(Seed);
    set_lattice_jq_isotropy_2d(Nx,Ny,Jbond);
    set_sequence_length(length);
    
    /*---------------Thermalization--------------*/
    for(int i_sample=0;i_sample<Nther;++i_sample){
        diagonal_update();
        if(Exc) diagonal_update_exchange();
        construct_link_vertex_list();
        loop_update();
        flip_bit_operator();
        if(Noo*buffer>length){
            length = (int)(Noo*buffer)+10;
            set_sequence_length(length);
        }
    }

    /*---------------Measurement-----------------*/
    set_estimator(n_obs,Nsample,Nblock);
    for(int k=0;k<Nblock;++k){
        for(int i_sample=0;i_sample<Nsample;++i_sample){
            diagonal_update();
            if(Exc) diagonal_update_exchange();
            construct_link_vertex_list();
            loop_update();
            flip_bit_operator();

            measure_mz(0,i_sample);
        }
        estimator_fileout(Filename);
    }
    
    return 0;
}
