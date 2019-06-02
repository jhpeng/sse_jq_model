#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

/*
**  define global variables here
*/

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

/*
** Nblock  : the total number of the block data
** Beta : inverse temperature
** rng  : gsl_rng
*/
int Nblock;
double Beta;
gsl_rng* rng;


/* -------------------------------------------------- **
** ---------------- SSE algorithm ------------------- **
*/ -------------------------------------------------- **

void propagate_state(sp){
    int i_bond = sp/6;
    int type   = sp%6;
    int spin;

    if(sp==-1 || type==0 || type==2) return;
    else if(type==1 || type==4){
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



/* --------------------------------------------------------- **
** ---------------- Setting the model ---------------------- **
*/ --------------------------------------------------------- **


