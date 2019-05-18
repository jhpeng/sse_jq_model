/*data_struct/mid_level.c
*/
#include <stdio.h>
#include <stdlib.h>

#include "low_level.h"
#include "mid_level.h"

lattice_struct* lattice_struct_alloc(int Nsite, int Nj, int Nq){
    lattice_struct* las = (lattice_struct*)malloc(sizeof(lattice_struct));
    las->Nsite = Nsite;
    las->Nj = Nj;
    las->Nq = Nq;

    las->sigma0 = int_array_alloc(Nsite);
    las->sigmap = int_array_alloc(Nsite);
    las->JQ = double_array_alloc(Nj+Nq);
    las->bond2index = int_array_alloc(4*(Nj+Nq));

    return las;
}

void lattice_struct_free(lattice_struct* las){
    int_array_free(las->sigma0);
    int_array_free(las->sigmap);
    double_array_free(las->JQ);
    int_array_free(las->bond2index);

    free(las);
}

int lattice_struct_get_Nsite(const lattice_struct* las){
    return las->Nsite;
}

int lattice_struct_get_Nj(const lattice_struct* las){
    return las->Nj;
}

int lattice_struct_get_Nq(const lattice_struct* las){
    return las->Nq;
}

int lattice_struct_get_sigma0(const lattice_struct* las, int index){
    return int_array_get(las->sigma0,index);
}

int lattice_struct_get_sigmap(const lattice_struct* las, int index){
    return int_array_get(las->sigma0,index);
}

double lattice_struct_get_bond(const lattice_struct* las, int i_bond){
    return double_array_get(las->JQ,i_bond);
}

void lattice_struct_get_bond2index(int* index, const lattice_struct* las, int i_bond){
    int i;
    for(i=0;i<4;++i){
        index[i] = int_array_get(las->bond2index,i+4*i_bond);
    }
}

void lattice_struct_set_sigma0(lattice_struct* las, int index, int x){
    int_array_set(las->sigma0,index,x);
}

void lattice_struct_set_sigmap(lattice_struct* las, int index, int x){
    int_array_set(las->sigmap,index,x);
}

void lattice_struct_set_bond(lattice_struct* las, int i_bond, double x){
    double_array_set(las->JQ,i_bond,x);
}

void lattice_struct_set_bond2index(lattice_struct* las, int i_bond, const int* index){
    int i;
    for(i=0;i<4;++i){
        int_array_set(las->bond2index,i+4*i_bond,index[i]);
    }
}



operator_sequence* operator_sequence_alloc(int L){
    operator_sequence* ops = (operator_sequence*)malloc(sizeof(operator_sequence));
    ops->L = L;
    ops->n = 0;
    ops->sequence = int_array_alloc(L);

    return ops;
}

void operator_sequence_free(operator_sequence* ops){
    int_array_free(ops->sequence);

    free(ops);
}

int operator_sequence_get_length(const operator_sequence* ops){
    return ops->L;
}

int operator_sequence_get_noo(const operator_sequence* ops){
    return ops->n;
}

int operator_sequence_get_sequence(const operator_sequence* ops, int p){
    return int_array_get(ops->sequence,p);
}

void operator_sequence_set_sequence(operator_sequence* ops, int p, int op){
    int_array_set(ops->sequence,p,op);
}

void operator_sequence_set_noo(operator_sequence* ops, int noo){
    ops->n = noo;
}

void operator_sequence_init(operator_sequence* ops){
    int_array_set_all(ops->sequence,-1);
    operator_sequence_set_noo(ops,0);
}


link_vertex* link_vertex_alloc(int L, int Nsite){
    link_vertex* lv = (link_vertex*)malloc(sizeof(link_vertex));
    lv->L = L;
    lv->Nsite = Nsite;
    lv->lvl = int_array_alloc(8*L);
    lv->v_first = int_array_alloc(Nsite);
    lv->v_last = int_array_alloc(Nsite);

    return lv;
}

void link_vertex_free(link_vertex* lv){
    int_array_free(lv->lvl);
    int_array_free(lv->v_first);
    int_array_free(lv->v_last);
    free(lv);
}

int link_vertex_get_length(const link_vertex* lv){
    return lv->L;
}

int link_vertex_get_Nsite(const link_vertex* lv){
    return lv->Nsite;
}

int link_vertex_get_list(const link_vertex* lv, int nu){
    return int_array_get(lv->lvl,nu);
}

int link_vertex_get_first(const link_vertex* lv, int index){
    return int_array_get(lv->v_first,index);
}

int link_vertex_get_last(const link_vertex* lv, int index){
    return int_array_get(lv->v_last,index);
}

void link_vertex_set_list(link_vertex* lv, int nu, int x){
    int_array_set(lv->lvl,nu,x);
}

void link_vertex_set_first(link_vertex* lv, int index, int x){
    int_array_set(lv->v_first,index,x);
}

void link_vertex_set_last(link_vertex* lv, int index, int x){
    int_array_set(lv->v_last,index,x);
}

void link_vertex_init(link_vertex* lv){
    int_array_set_all(lv->lvl,-1);
    int_array_set_all(lv->v_first,-1);
    int_array_set_all(lv->v_last,-1);
}


observable* observable_alloc(int Nobs, int Nsample){
    observable* obs = (observable*)malloc(sizeof(observable));
    obs->Nobs = Nobs;
    obs->Nsample = Nsample;
    obs->data = double_array_alloc(Nobs*Nsample);
    obs->init = int_array_alloc(Nobs*Nsample);

    return obs;
}

void observable_free(observable* obs){
    double_array_free(obs->data);
    int_array_free(obs->init);
    free(obs);
}

int observable_get_Nobs(const observable* obs){
    return obs->Nobs;
}

int observable_get_Nsample(const observable* obs){
    return obs->Nsample;
}

int observable_get_data(const observable* obs, int i_obs, int i_sample){
    return double_array_get(obs->data,i_obs+obs->Nobs*i_sample);
}

int observable_get_init(const observable* obs, int i_obs, int i_sample){
    return int_array_get(obs->init,i_obs+obs->Nobs*i_sample);
}

void observable_set_data(observable* obs, int i_obs, int i_sample, double x){
    double_array_set(obs->data,i_obs+obs->Nobs*i_sample,x);
    int_array_set(obs->init,i_obs+obs->Nobs*i_sample,1);
}

void observable_init(observable* obs){
    int i;
    for(i=0;i<(obs->Nobs*obs->Nsample);++i){
        int_array_set(obs->init,i,0);
    }
}

int observable_check_finish(const observable* obs){
    int i,check=0;
    for(i=0;i<(obs->Nobs*obs->Nsample);++i){
        if(int_array_get(obs->init,i)^1) return check;
    }

    check=1;
    return check;
}

#if 0
int main(){
    int Nsite = 100;
    int Nj = 200;
    int Nq = 200;
    int L = 10000;
    int Nobs = 8;
    int Nsample = 2000;

    for(int i=0;i<1000;++i){
        lattice_struct* las = lattice_struct_alloc(Nsite,Nj,Nq);
        operator_sequence* ops = operator_sequence_alloc(L);
        link_vertex* lv = link_vertex_alloc(L,Nsite);
        observable* obs = observable_alloc(Nobs,Nsample);

        operator_sequence_init(ops);
        link_vertex_init(lv);
        observable_init(obs);

        lattice_struct_free(las);
        operator_sequence_free(ops);
        link_vertex_free(lv);
        observable_free(obs);
    }
    int check = low_level_check_memleak();
    printf("%d\n",check);

    return 0;
}
#endif
