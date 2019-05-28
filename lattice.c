/* algorithm/lattice.c
*/
#include <gsl/gsl_rng.h>

#include "low_level.h"
#include "mid_level.h"

void lattice_struct_sync_sigmap(lattice_struct* las){
    int i,spin,Nsite = lattice_struct_get_Nsite(las);
    for(i=0;i<Nsite;++i){
        spin = lattice_struct_get_sigma0(las,i);
        lattice_struct_set_sigmap(las,i,spin);
    }
}

void lattice_struct_propagate_state(lattice_struct* las, int sp){
    int i_bond = sp/6;
    int type   = sp%6;
    int spin;

    if(sp==-1) return;
    else if(type==1 || type==4){
        int index[4];
        lattice_struct_get_bond2index(index,las,i_bond);
        spin = lattice_struct_get_sigmap(las,index[0]);
        lattice_struct_set_sigmap(las,index[0],-1*spin);
        spin = lattice_struct_get_sigmap(las,index[1]);
        lattice_struct_set_sigmap(las,index[1],-1*spin);
    }
    else if(type==3){
        int index[4];
        lattice_struct_get_bond2index(index,las,i_bond);
        spin = lattice_struct_get_sigmap(las,index[2]);
        lattice_struct_set_sigmap(las,index[2],-1*spin);
        spin = lattice_struct_get_sigmap(las,index[3]);
        lattice_struct_set_sigmap(las,index[3],-1*spin);
    }
    else if(type==5){
        int index[4];
        lattice_struct_get_bond2index(index,las,i_bond);
        spin = lattice_struct_get_sigmap(las,index[0]);
        lattice_struct_set_sigmap(las,index[0],-1*spin);
        spin = lattice_struct_get_sigmap(las,index[1]);
        lattice_struct_set_sigmap(las,index[1],-1*spin);
        spin = lattice_struct_get_sigmap(las,index[2]);
        lattice_struct_set_sigmap(las,index[2],-1*spin);
        spin = lattice_struct_get_sigmap(las,index[3]);
        lattice_struct_set_sigmap(las,index[3],-1*spin);
    }
}

int lattice_struct_check_propagate_state(lattice_struct* las, operator_sequence* ops){
    int p,sp,L = operator_sequence_get_length(ops);
    int index,Nsite = lattice_struct_get_Nsite(las);
    lattice_struct_sync_sigmap(las);
    
    for(p=0;p<L;++p){
        sp = operator_sequence_get_sequence(ops,p);
        lattice_struct_propagate_state(las,sp);
    }
    int check=0;

    for(index=0;index<Nsite;++index){
        int sigma0 = lattice_struct_get_sigma0(las,index);
        int sigmap = lattice_struct_get_sigmap(las,index);
        if(sigma0!=sigmap){
            check=1;
            return check;
        }
    }
    return check;
}

static void lattice_struct_hot_start(lattice_struct* las, gsl_rng* rng){
    int i,Nsite = lattice_struct_get_Nsite(las);
    for(i=0;i<Nsite;++i){
        if(gsl_rng_uniform_pos(rng)<0.5) lattice_struct_set_sigma0(las,i,1);
        else lattice_struct_set_sigma0(las,i,-1);
    }
}

lattice_struct* lattice_struct_create_model_plaquette_2d(int Nx, int Ny, double J, gsl_rng* rng){
    if(Nx%2!=0 || Ny%2!=0){
        printf("Nx and Ny must be the miltiply of 2!\n");
    }
    int Nsite = Nx*Ny;
    int Nj = 2*Nx*Ny;
    int Nq = Nx*Ny/2;
    int i,j,k,t,q,index[4];
    lattice_struct* las = lattice_struct_alloc(Nsite,Nj,Nq);

    int i_bond;
    for(i_bond=0;i_bond<Nj;++i_bond){
        t = i_bond/Nsite;
        k = i_bond%Nsite;
        i = k%Nx;
        j = k/Nx;
        i+= t^1;
        j+= t;
        i = i%Nx;
        j = j%Ny;
        index[0] = k;
        index[1] = i+j*Nx;
        index[2] = -1;
        index[3] = -1;

        lattice_struct_set_bond2index(las,i_bond,index);
        lattice_struct_set_bond(las,i_bond,J);
    }
    for(i_bond=Nj;i_bond<(Nj+Nq);++i_bond){
        q = i_bond-Nj;
        t = q%2;
        i = ((q/2)%(Nx/2))*2;
        j = ((q/2)/(Nx/2))*2;
        if(t==0){
            index[0] = i     + j*Nx;
            index[1] = (i+1) + j*Nx;
            index[2] = i     + (j+1)*Nx;
            index[3] = (i+1) + (j+1)*Nx;
        }
        else if(t==1){
            index[0] = i     + j*Nx;
            index[1] = i     + (j+1)*Nx;
            index[2] = (i+1) + j*Nx;
            index[3] = (i+1) + (j+1)*Nx;
        }

        lattice_struct_set_bond2index(las,i_bond,index);
        lattice_struct_set_bond(las,i_bond,1);
    }

    lattice_struct_hot_start(las,rng);

    return las;
}

lattice_struct* lattice_struct_create_model_isotropy_2d(int Nx, int Ny, double J, gsl_rng* rng){
    if(Nx%2!=0 || Ny%2!=0){
        printf("Nx and Ny must be the miltiply of 2!\n");
    }
    int Nsite = Nx*Ny;
    int Nj = 2*Nx*Ny;
    int Nq = 2*Nx*Ny;
    int i,j,t,q,index[4];
    lattice_struct* las = lattice_struct_alloc(Nsite,Nj,Nq);

    for(int i_bond=0;i_bond<(4*Nsite);++i_bond){
        t = i_bond%Nsite;
        q = i_bond/Nsite;
        i = t%Nx;
        j = t/Nx;

        /*Initialize the index*/
        index[0]=-1;
        index[1]=-1;
        index[2]=-1;
        index[3]=-1;

        if(q==0){
            index[0] = i+Nx*j;
            index[1] = ((i+1)%Nx)+Nx*j;
            index[2] = -1;
            index[3] = -1;
            lattice_struct_set_bond(las,i_bond,J);
        }
        else if(q==1){
            index[0] = i+Nx*j;
            index[1] = i+Nx*((j+1)%Ny);
            index[2] = -1;
            index[3] = -1;
            lattice_struct_set_bond(las,i_bond,J);
        }
        else if(q==2){
            index[0] = i+Nx*j;
            index[1] = ((i+1)%Nx)+Nx*j;
            index[2] = i+Nx*((j+1)%Ny);
            index[3] = ((i+1)%Nx)+Nx*((j+1)%Ny);
            lattice_struct_set_bond(las,i_bond,1);
        }
        else if(q==3){
            index[0] = i+Nx*j;
            index[1] = i+Nx*((j+1)%Ny);
            index[2] = ((i+1)%Nx)+Nx*j;
            index[3] = ((i+1)%Nx)+Nx*((j+1)%Ny);
            lattice_struct_set_bond(las,i_bond,1);
        }
        lattice_struct_set_bond2index(las,i_bond,index);
    }
    lattice_struct_hot_start(las,rng);

    return las;
}
