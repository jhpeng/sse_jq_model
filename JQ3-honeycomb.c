#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

/*
**  define global variables here
*/

/* The flip rule is the instruction for flipping operator.
**   | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
** -------------------------------------------
** a | 1 | 0 | 3 | 2 | 6 | 7 | 4 | 5 | 9 | 8 |
** b | 0 | 1 | 4 | 6 | 2 | 8 | 3 | 9 | 5 | 7 |
** c | 0 | 1 | 5 | 7 | 8 | 2 | 9 | 3 | 4 | 6 |
** d | 1 | 0 | 3 | 2 | 6 | 7 | 4 | 5 | 9 | 8 |
** e | 0 | 1 | 4 | 6 | 2 | 8 | 3 | 9 | 5 | 7 |
** f | 0 | 1 | 5 | 7 | 8 | 2 | 9 | 3 | 4 | 6 |
**
**
**
**
**  |    |        |    |        |    |
** d|   d|       e|   e|       f|   f|
** --------      --------      --------
** |\\\\\\|------|\\\\\\|------|\\\\\\|
** --------      --------      --------
**  |    |        |    |        |    |
** a|   a|       b|   b|       c|   c|
**
** label : 0   1   2   3   4   5   6   7   8   9
** type  : d   o   ddd odd dod ddo ood odo doo ooo
*/
static int flip_rule[60]={
    1,0,3,2,6,7,4,5,9,8,
    0,1,4,6,2,8,3,9,5,7,
    0,1,5,7,8,2,9,3,4,6,
    1,0,3,2,6,7,4,5,9,8,
    0,1,4,6,2,8,3,9,5,7,
    0,1,5,7,8,2,9,3,4,6};

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
** Bondst : the pointer to the bond strenght
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
** Beta    : inverse temperature
** rng     : gsl_rng
** Mode    : the mode for calculate observable
**             0 -> normal scheme
**             1 -> beta doubling
*/
int Nblock,Nther,Seed;
double Beta,Jbond,Qbond,P;
gsl_rng* rng;
int Mode,LatticeType;
int Nit;


/* -------------------------------------------------- **
** ---------------- SSE algorithm ------------------- **
** -------------------------------------------------- */

int propagate_rule[60]={ 
     1, 1, 1, 1, 1, 1,
    -1,-1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1,
    -1,-1, 1, 1, 1, 1,
     1, 1,-1,-1, 1, 1,
     1, 1, 1, 1,-1,-1,
    -1,-1,-1,-1, 1, 1,
    -1,-1, 1, 1,-1,-1,
     1, 1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1};

void propagate_state(int sp){
    if(sp==-1) return;

    int type   = sp%10;
    if(type==0 || type==2) return;
    
    int i_bond = sp/10;
    if(type==1){
        Sigmap[Bond2index[i_bond*6+0]] *= -1;
        Sigmap[Bond2index[i_bond*6+1]] *= -1;
    }
    else{
        Sigmap[Bond2index[i_bond*6+0]] *= propagate_rule[type*6+0];
        Sigmap[Bond2index[i_bond*6+1]] *= propagate_rule[type*6+1];
        Sigmap[Bond2index[i_bond*6+2]] *= propagate_rule[type*6+2];
        Sigmap[Bond2index[i_bond*6+3]] *= propagate_rule[type*6+3];
        Sigmap[Bond2index[i_bond*6+4]] *= propagate_rule[type*6+4];
        Sigmap[Bond2index[i_bond*6+5]] *= propagate_rule[type*6+5];
    }
}

void diagonal_update(){
    int i_bond,s1,s2,s3,s4,s5,s6;
    double dis;

    for(int i=0;i<Nsite;++i){
        Sigmap[i] = Sigma0[i];
    }

    for(int p=0;p<L;++p){
        if(Sequence[p]==-1){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*(Nj+Nq));
            if(i_bond<Nj){
                s1 = Sigmap[Bond2index[i_bond*6+0]];
                s2 = Sigmap[Bond2index[i_bond*6+1]];
                if(s1!=s2){
                    dis = gsl_rng_uniform_pos(rng);
                    if(dis*2*(L-Noo)<Beta*Bondst[i_bond]*(Nj+Nq)){
                        Sequence[p] = i_bond*10;
                        Noo++;
                    }
                }
            }
            else {
                s1 = Sigmap[Bond2index[i_bond*6+0]];
                s2 = Sigmap[Bond2index[i_bond*6+1]];
                s3 = Sigmap[Bond2index[i_bond*6+2]];
                s4 = Sigmap[Bond2index[i_bond*6+3]];
                s5 = Sigmap[Bond2index[i_bond*6+4]];
                s6 = Sigmap[Bond2index[i_bond*6+5]];
                if((s1!=s2 && s3!=s4) && s5!=s6){
                    dis = gsl_rng_uniform_pos(rng);
                    if(dis*8*(L-Noo)<Beta*Bondst[i_bond]*(Nj+Nq)){
                        Sequence[p] = i_bond*10+2;
                        Noo++;
                    }
                }
            }
        }
        else if(Sequence[p]%10==0){
            i_bond = Sequence[p]/10;
            dis = gsl_rng_uniform_pos(rng);
            if(Beta*(Nj+Nq)*Bondst[i_bond]*dis<2*(L-Noo+1)){
                Sequence[p]=-1;
                Noo--;
            }
        }
        else if(Sequence[p]%10==2){
            i_bond = Sequence[p]/10;
            dis = gsl_rng_uniform_pos(rng);
            if(Beta*(Nj+Nq)*Bondst[i_bond]*dis<8*(L-Noo+1)){
                Sequence[p]=-1;
                Noo--;
            }
        }
        else propagate_state(Sequence[p]);
    }
}

void construct_link_vertex_list(){
    for(int i=0;i<(12*L);++i) Linkv[i]=-1;
    for(int i=0;i<Nsite;++i){
        Vfirst[i]=-1;
        Vlast[i] =-1;
    }

    int i_bond,index,nu0,nu1;
    for(int p=0;p<L;++p){
        if(Sequence[p]!=-1){
            i_bond = Sequence[p]/10;
            for(int i=0;i<6;++i){
                index = Bond2index[i_bond*6+i];
                if(index!=-1){
                    nu0 = 12*p+i;
                    nu1 = Vlast[index];
                    if(nu1!=-1){
                        Linkv[nu0] = nu1;
                        Linkv[nu1] = nu0;
                    }
                    else Vfirst[index] = nu0;

                    Vlast[index] = nu0+6;
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

    for(nu0=0;nu0<(L*12);nu0+=2){
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
    int nu,p,i,i_flip_rule,type,i_bond,index;
    
    for(nu=0;nu<(12*L);nu+=2){
        if(Linkv[nu]==-2){
            p = nu/12;
            i = nu%12;
            type = Sequence[p]%10;
            i_bond = Sequence[p]/10;
            i_flip_rule = (i/2)*10+type;
            Sequence[p] = i_bond*10+flip_rule[i_flip_rule];
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


void beta_doubling(){
    int* seq = (int*)malloc(2*L*sizeof(int));
    for(int i=0;i<L;++i){
        seq[i] = Sequence[i%L];
    }
    for(int i=0;i<L;++i){
        seq[i+L] = Sequence[L-i-1];
    }
    free(Sequence);
    free(Linkv);

    Sequence = seq;
    Linkv = (int*)malloc(24*L*sizeof(int));

    Noo  *=2;
    L    *=2;
    Beta *=2;
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

    int Ma=0;
    int Mb=0;
    int Dx=0;
    int Dy=0;
    for(int j=0;j<Ny;j++){
        for(int i=0;i<Nx;i++){
            if((i+j)%2){
                Ma += Sigma0[i+j*Nx];
                Dx += Sigma0[i+j*Nx]*Sigma0[(i+1)%Nx+j*Nx];
                Dy += Sigma0[i+j*Nx]*Sigma0[i+((j+1)%Ny)*Nx];
            } else {
                Mb += Sigma0[i+j*Nx];
                Dx -= Sigma0[i+j*Nx]*Sigma0[(i+1)%Nx+j*Nx];
                Dy -= Sigma0[i+j*Nx]*Sigma0[i+((j+1)%Ny)*Nx];
            }
        }
    }

    FILE* ofile = fopen(filename,"a");
    fprintf(ofile,"%.5e ",Beta);
    for(int i=0;i<Nobs;++i) fprintf(ofile,"%.16e ",mean[i]);
    fprintf(ofile,"%d ",Ma);
    fprintf(ofile,"%d ",Mb);
    fprintf(ofile,"%d ",Dx);
    fprintf(ofile,"%d ",Dy);
    fprintf(ofile,"\n");
    fclose(ofile);
}

int* stagger_factor;

void measurement(int i_sample){
    double mz2;
    double mz=0;
    double ms=0;
    double msx=0;
    double ms1=0;
    double ms2=0;
    double ms4=0;
    if(stagger_factor==NULL){
        stagger_factor = (int*)malloc(sizeof(int)*Nx*Ny);
        for(int j=0;j<Ny;j++){
            for(int i=0;i<Nx;i++){
                stagger_factor[i+j*Nx] = (i+j)%2*2-1;
            }
        }
    }
    for(int i=0;i<Nsite;++i){
        mz+=Sigma0[i];
        ms+=Sigma0[i]*stagger_factor[i];
    }

    /*-------------- propagate state ----------------*/
    for(int i=0;i<Nsite;i++) Sigmap[i] = Sigma0[i];
    int sp,i_bond,type,q;
    int s1,s2,s3,s4,s5,s6;
    double wx=0;
    double wy=0;
    for(int p=0;p<L;++p){
        sp = Sequence[p];
        if(sp!=-1){
            type = sp%10;
            i_bond = sp/10;

            q=0;
            if(i_bond<Nsite) q=0;
            else if(i_bond>=Nsite && i_bond<Nj) q=1;
            else if(i_bond>=Nj && i_bond<Nj+Nsite/2)q=2;
            else if(i_bond>=Nj+Nsite/2 && i_bond<Nj+Nq)q=3;
            else exit(-1);

            s1 = Bond2index[i_bond*6+0];
            s2 = Bond2index[i_bond*6+1];
            s3 = Bond2index[i_bond*6+2];
            s4 = Bond2index[i_bond*6+3];
            s5 = Bond2index[i_bond*6+4];
            s6 = Bond2index[i_bond*6+5];

            //if(q==0 && type==1) wx += Sigmap[s1];
            //else if(q==1 && type==1) wy += Sigmap[s1];

            if(type==1){
                ms += -2*Sigmap[s1]*stagger_factor[s1];
                ms += -2*Sigmap[s2]*stagger_factor[s2];

                if(q==0) wx += Sigmap[s1];
                else if(q==1) wy += Sigmap[s1];

                Sigmap[s1] *= -1;
                Sigmap[s2] *= -1;
            } else if(type==3){
                ms += -2*Sigmap[s1]*stagger_factor[s1];
                ms += -2*Sigmap[s2]*stagger_factor[s2];

                wy += Sigmap[s1];

                Sigmap[s1] *= -1;
                Sigmap[s2] *= -1;
            } else if(type==4) {
                ms += -2*Sigmap[s3]*stagger_factor[s3];
                ms += -2*Sigmap[s4]*stagger_factor[s4];

                wx += Sigmap[s3];

                Sigmap[s3] *= -1;
                Sigmap[s4] *= -1;
            } else if(type==5) {
                ms += -2*Sigmap[s5]*stagger_factor[s5];
                ms += -2*Sigmap[s6]*stagger_factor[s6];

                wx += Sigmap[s5];

                Sigmap[s5] *= -1;
                Sigmap[s6] *= -1;
            } else if(type==6) {
                ms += -2*Sigmap[s1]*stagger_factor[s1];
                ms += -2*Sigmap[s2]*stagger_factor[s2];
                ms += -2*Sigmap[s3]*stagger_factor[s3];
                ms += -2*Sigmap[s4]*stagger_factor[s4];

                wy += Sigmap[s1];
                wx += Sigmap[s3];

                Sigmap[s1] *= -1;
                Sigmap[s2] *= -1;
                Sigmap[s3] *= -1;
                Sigmap[s4] *= -1;
            } else if(type==7) {
                ms += -2*Sigmap[s1]*stagger_factor[s1];
                ms += -2*Sigmap[s2]*stagger_factor[s2];
                ms += -2*Sigmap[s5]*stagger_factor[s5];
                ms += -2*Sigmap[s6]*stagger_factor[s6];

                wy += Sigmap[s1];
                wx += Sigmap[s5];

                Sigmap[s1] *= -1;
                Sigmap[s2] *= -1;
                Sigmap[s5] *= -1;
                Sigmap[s6] *= -1;
            } else if(type==8) {
                ms += -2*Sigmap[s3]*stagger_factor[s3];
                ms += -2*Sigmap[s4]*stagger_factor[s4];
                ms += -2*Sigmap[s5]*stagger_factor[s5];
                ms += -2*Sigmap[s6]*stagger_factor[s6];

                wx += Sigmap[s3];
                wx += Sigmap[s5];

                Sigmap[s3] *= -1;
                Sigmap[s4] *= -1;
                Sigmap[s5] *= -1;
                Sigmap[s6] *= -1;
            } else if(type==9) {
                ms += -2*Sigmap[s1]*stagger_factor[s1];
                ms += -2*Sigmap[s2]*stagger_factor[s2];
                ms += -2*Sigmap[s3]*stagger_factor[s3];
                ms += -2*Sigmap[s4]*stagger_factor[s4];
                ms += -2*Sigmap[s5]*stagger_factor[s5];
                ms += -2*Sigmap[s6]*stagger_factor[s6];

                wy += Sigmap[s1];
                wx += Sigmap[s3];
                wx += Sigmap[s5];

                Sigmap[s1] *= -1;
                Sigmap[s2] *= -1;
                Sigmap[s3] *= -1;
                Sigmap[s4] *= -1;
                Sigmap[s5] *= -1;
                Sigmap[s6] *= -1;
            }

            msx += ms;
            ms1 += fabs(ms);
            ms2 += ms*ms;
            ms4 += ms*ms*ms*ms;
        }
    }

    if(Noo!=0) {
        msx = Beta*(msx*msx+ms2)/Noo/(Noo+1)/Nsite/Nsite*0.25;
        ms1 = ms1*0.5/Nsite/Noo;
        ms2 = ms2*0.25/Nsite/Nsite/Noo;
        ms4 = ms4*0.0625/Nsite/Nsite/Nsite/Nsite/Noo;
    } else {
        msx = 0;
        ms1 = 0;
        ms2 = 0;
        ms4 = 0;
    }

    mz2 = mz*mz*0.25;
    Data[Nobs*i_sample+0] = ms1;
    Data[Nobs*i_sample+1] = ms2;
    Data[Nobs*i_sample+2] = ms4;
    Data[Nobs*i_sample+3] = msx;
    Data[Nobs*i_sample+4] = mz2/Nsite;
    Data[Nobs*i_sample+5] = wx*wx/Nx/Nx;
    Data[Nobs*i_sample+6] = wy*wy/Ny/Ny;
    Data[Nobs*i_sample+7] = Noo;
    Data[Nobs*i_sample+8] = (double)Noo*(double)Noo;
}

/* --------------------------------------------------------- **
** ---------------- Setting the model ---------------------- **
** --------------------------------------------------------- */


void set_random_number(int seed){
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
}

void set_lattice_jq3_honeycomb(int nx, int ny, double qbond){
    int i,j,t;
    Nsite = nx*ny;
    Nj = 3*(Nsite/2);
    Nq = Nsite;

    Sigma0 = (int*)malloc(Nsite*sizeof(int));
    Sigmap = (int*)malloc(Nsite*sizeof(int));
    Vfirst = (int*)malloc(Nsite*sizeof(int));
    Vlast  = (int*)malloc(Nsite*sizeof(int));

    Bond2index = (int*)malloc((Nj+Nq)*6*sizeof(int));
    Bondst = (double*)malloc((Nj+Nq)*sizeof(double));

    for(int i_bond=0;i_bond<Nsite;++i_bond){
        t = i_bond;
        i = t%nx;
        j = t/nx;

        Bond2index[i_bond*6+0] = i+nx*j;
        Bond2index[i_bond*6+1] = ((i+1)%nx)+nx*j;
        Bond2index[i_bond*6+2] = -1;
        Bond2index[i_bond*6+3] = -1;
        Bond2index[i_bond*6+4] = -1;
        Bond2index[i_bond*6+5] = -1;
        Bondst[i_bond] = 1.0;
    }
    for(int i_bond=Nsite;i_bond<Nj;++i_bond){
        t = i_bond-Nsite;
        j = t/(nx/2);
        i = (t%(nx/2))*2+j%2;

        Bond2index[i_bond*6+0] = i+nx*j;
        Bond2index[i_bond*6+1] = ((i+0)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+2] = -1;
        Bond2index[i_bond*6+3] = -1;
        Bond2index[i_bond*6+4] = -1;
        Bond2index[i_bond*6+5] = -1;
        Bondst[i_bond] = 1.0;
    }
    for(int i_bond=Nj;i_bond<Nj+Nsite/2;++i_bond){
        t = i_bond-Nj;
        j = t/(nx/2);
        i = (t%(nx/2))*2+j%2;

        Bond2index[i_bond*6+0] = i+nx*j;
        Bond2index[i_bond*6+1] = ((i+0)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+2] = ((i+1)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+3] = ((i+2)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+4] = ((i+1)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+5] = ((i+2)%nx)+nx*((j+1)%ny);
        Bondst[i_bond] = qbond;
    }
    for(int i_bond=Nj+Nsite/2;i_bond<Nj+Nq;++i_bond){
        t = i_bond-Nj-Nsite/2;
        j = t/(nx/2);
        i = (t%(nx/2))*2+j%2;

        Bond2index[i_bond*6+0] = ((i+2)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+1] = ((i+2)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+2] = ((i+0)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+3] = ((i+1)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+4] = ((i+0)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+5] = ((i+1)%nx)+nx*((j+1)%ny);
        Bondst[i_bond] = qbond;
    }

    for(i=0;i<Nsite;++i){
        if(gsl_rng_uniform_pos(rng)<0.5) Sigma0[i]=1;
        else Sigma0[i]=-1;
    }
}

void set_lattice_jq3_honeycomb_spin_half_impurity(int nx, int ny, double qbond){
    int i,j,t;
    Nsite = nx*ny;
    Nj = 3*(Nsite/2);
    Nq = Nsite;

    Sigma0 = (int*)malloc(Nsite*sizeof(int));
    Sigmap = (int*)malloc(Nsite*sizeof(int));
    Vfirst = (int*)malloc(Nsite*sizeof(int));
    Vlast  = (int*)malloc(Nsite*sizeof(int));

    Bond2index = (int*)malloc((Nj+Nq)*6*sizeof(int));
    Bondst = (double*)malloc((Nj+Nq)*sizeof(double));

    for(int i_bond=0;i_bond<Nsite;++i_bond){
        t = i_bond;
        i = t%nx;
        j = t/nx;

        Bond2index[i_bond*6+0] = i+nx*j;
        Bond2index[i_bond*6+1] = ((i+1)%nx)+nx*j;
        Bond2index[i_bond*6+2] = -1;
        Bond2index[i_bond*6+3] = -1;
        Bond2index[i_bond*6+4] = -1;
        Bond2index[i_bond*6+5] = -1;
        Bondst[i_bond] = 1.0;
    }
    for(int i_bond=Nsite;i_bond<Nj;++i_bond){
        t = i_bond-Nsite;
        j = t/(nx/2);
        i = (t%(nx/2))*2+j%2;

        Bond2index[i_bond*6+0] = i+nx*j;
        Bond2index[i_bond*6+1] = ((i+0)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+2] = -1;
        Bond2index[i_bond*6+3] = -1;
        Bond2index[i_bond*6+4] = -1;
        Bond2index[i_bond*6+5] = -1;
        Bondst[i_bond] = 1.0;
    }
    for(int i_bond=Nj;i_bond<Nj+Nsite/2;++i_bond){
        t = i_bond-Nj;
        j = t/(nx/2);
        i = (t%(nx/2))*2+j%2;

        Bond2index[i_bond*6+0] = i+nx*j;
        Bond2index[i_bond*6+1] = ((i+0)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+2] = ((i+1)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+3] = ((i+2)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+4] = ((i+1)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+5] = ((i+2)%nx)+nx*((j+1)%ny);
        Bondst[i_bond] = qbond;
    }
    for(int i_bond=Nj+Nsite/2;i_bond<Nj+Nq;++i_bond){
        t = i_bond-Nj-Nsite/2;
        j = t/(nx/2);
        i = (t%(nx/2))*2+j%2;

        Bond2index[i_bond*6+0] = ((i+2)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+1] = ((i+2)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+2] = ((i+0)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+3] = ((i+1)%nx)+nx*((j+0)%ny);
        Bond2index[i_bond*6+4] = ((i+0)%nx)+nx*((j+1)%ny);
        Bond2index[i_bond*6+5] = ((i+1)%nx)+nx*((j+1)%ny);
        Bondst[i_bond] = qbond;
    }

    for(i=0;i<Nsite;++i){
        if(gsl_rng_uniform_pos(rng)<0.5) Sigma0[i]=1;
        else Sigma0[i]=-1;
    }

    for(int i_bond=0;i_bond<(Nj+Nq);i_bond++){
        for(i=0;i<6;i++){
            if(Bond2index[i_bond*6+i]==Nsite-1){
                Bondst[i_bond] = 0;
            }
        }
    }
    Nsite = Nsite-1;
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
    Linkv = (int*)malloc(12*length*sizeof(int));

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

int Help;
void set_opt(int argc, char **argv)
{
    int c;
    while((c=getopt(argc,argv,"hx:y:q:b:n:k:t:s:f:m:l:p:"))!=-1){
        switch(c){
            case 'h':
                Help=1;
                printf("usage:\n");
                printf("\t-h print this help\n");
                printf("\t-l lattice type for the simulation\n");
                printf("\t\t 0 : 2d JQ3 honeycomb\n");
                printf("\t\t 0 : 2d JQ3 honeycomb (spin-half impurity)\n");
                printf("\t-m mode for calculate observable\n");
                printf("\t\t 0 : normal scheme\n");
                printf("\t\t 1 : beta-doubling scheme\n");
                printf("\t\t 2 : beta increasing scheme\n");
                printf("\t\t 3 : beta increasing scheme without propagate state\n");
                printf("\t-x <length of x> default 8\n");
                printf("\t-y <length of y> default 8\n");
                printf("\t-q <Q3/J ratio> default 1.0\n");
                printf("\t-b <beta> default 4.0\n");
                printf("\t-n <Nsample> default 2000\n");
                printf("\t-k <Nblock> default 50\n");
                printf("\t-t <Nther> default 2000\n");
                printf("\t-e number of iteration for beta-doubing scheme\n");
                printf("\t-s <seed of random number generator> default 1\n");
                printf("\t-f <the file name of output data> default \"test.txt\"\n");
                break;
            case 'l':
                LatticeType=atoi(optarg);
                break;
            case 'm':
                Mode=atoi(optarg);
                break;
            case 'x':
                Nx=atoi(optarg);
                break;
            case 'y':
                Ny=atoi(optarg);
                break;
            case 'q':
                Qbond=atof(optarg);
                break;
            case 'p':
                P=atof(optarg);
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
            case 'e':
                Nit=atoi(optarg);
                break;
            case 's':
                Seed=atoi(optarg);
                break;
            case 'f':
                strcpy(Filename,optarg);
                break;
        }
    }
}


/* ----------------------------------------------- **
** ------------------ main ----------------------- **
** ----------------------------------------------- */ 

int main(int argc, char** argv){
    int length=1000;
    int n_obs=9;
    double buffer=1.3;

    /*--------------default value----------------*/
    Beta = 8;
    Seed = 9237912;
    Nx = 8;
    Ny = 8;
    Qbond = 1.0;
    P = 0;
    Nther = 20000;
    Nsample = 2000;
    Nblock = 50;
    Mode = 0;
    LatticeType = 0;
    Nit = 5;




    /*-------------setting lattice---------------*/
    set_opt(argc,argv);
    if(Help) return 0;

    set_random_number(Seed);

    if(LatticeType==0) set_lattice_jq3_honeycomb(Nx,Ny,Qbond);
    else if(LatticeType==1) set_lattice_jq3_honeycomb_spin_half_impurity(Nx,Ny,Qbond);
    set_sequence_length(length);


    /*---------------lattice debug----------------*/
    if(0){
        int i_bond;
        printf("--------------------------------\n");
        for(i_bond=0;i_bond<Nsite;i_bond++){
            for(int i=0;i<6;i++) printf("%d ",Bond2index[i_bond*6+i]);
            printf("\n");
        }
        printf("--------------------------------\n");
        for(i_bond=Nsite;i_bond<Nj;i_bond++){
            for(int i=0;i<6;i++) printf("%d ",Bond2index[i_bond*6+i]);
            printf("\n");
        }
        printf("--------------------------------\n");
        for(i_bond=Nj;i_bond<Nj+Nsite/2;i_bond++){
            for(int i=0;i<6;i++) printf("%d ",Bond2index[i_bond*6+i]);
            printf("\n");
        }
        printf("--------------------------------\n");
        for(i_bond=Nj+Nsite/2;i_bond<Nj+Nq;i_bond++){
            for(int i=0;i<6;i++) printf("%d ",Bond2index[i_bond*6+i]);
            printf("\n");
        }
        return 0 ;
    }   



/**************************************************************/
/*********************** Normal Scheme ************************/
/**************************************************************/
    if(Mode==0){
        /*---------------Thermalization--------------*/
        for(int i_sample=0;i_sample<Nther;++i_sample){
            diagonal_update();
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
                construct_link_vertex_list();
                loop_update();
                flip_bit_operator();

                measurement(i_sample);
            }
            estimator_fileout(Filename);
        }
    }
/**************************************************************/
/******************** Beta-doubling Scheme ********************/
/**************************************************************/
    else if(Mode==1){
        for(int it=0;it<2*Nit;++it){
            /*---------------Thermalization--------------*/
            for(int i_sample=0;i_sample<Nther;++i_sample){
                diagonal_update();
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
                    construct_link_vertex_list();
                    loop_update();
                    flip_bit_operator();

                    measurement(i_sample);
                }
                estimator_fileout(Filename);
            }
            if(it%2==1) beta_doubling();
        }
    }

    return 0;
}
