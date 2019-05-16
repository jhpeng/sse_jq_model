#include <gsl/gsl_rng>

#include "low_level.h"
#include "mid_level.h"

void lattice_struct_sync_sigmap(lattice_struct* las){
    int i,spin,Nsite = lattice_struct_get_Nsite(las);
    for(i=0;i<Nsite;++i){
        spin = lattice_struct_get_sigma0(las,i);
        lattice_struct_set_sigmap(las,i,spin);
    }
}

void lattice_struct_propagate_state(lattice_strcut* las, int op){
    int i_bond = op/6;
    int type   = op%6;

    lattice_struct_get_bond2index(index,las,i_bond);
    if(type==1 or type==4){
        int index[4];
        lattice_struct_get_bond2index(index,las,i_bond);
        lattice_struct_set_sigmap(las,index[0],-1*spin);
        lattice_struct_set_sigmap(las,index[1],-1*spin);
    }
    else if(type==3){
        int index[4];
        lattice_struct_get_bond2index(index,las,i_bond);
        lattice_struct_set_sigmap(las,index[2],-1*spin);
        lattice_struct_set_sigmap(las,index[3],-1*spin);
    }
    else if(type==5){
        int index[4];
        lattice_struct_get_bond2index(index,las,i_bond);
        lattice_struct_set_sigmap(las,index[0],-1*spin);
        lattice_struct_set_sigmap(las,index[1],-1*spin);
        lattice_struct_set_sigmap(las,index[2],-1*spin);
        lattice_struct_set_sigmap(las,index[3],-1*spin);
    }
}
