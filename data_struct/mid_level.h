/*data_struct/mid_level.h
*/
#ifndef MID_LEVEL_H
#define MID_LEVEL_H

/* Note on lattice struct
** Nsite  : total number of spins
** Nj     : total number of J bonds
** Nq     : total number of Q bond-pairs
** sigma0 : storage the current state of the spin
** sigmap : working place of propagate state of spin
** JQ     : Storage the bonds' strength 0~(Nj-1) is for J
**          and Nj~(Nj+Nq-1) is for Q.
** bind2index : the corresponding spin index for bond or 
                bond-pair. leg+4*b
*/
typedef struct lattice_struct{
    int Nsite;
    int Nj;
    int Nq;
    int_array* sigma0;
    int_array* sigmap;
    double_array* JQ;
    int_array* bond2index;
} lattice_struct;

lattice_struct* lattice_struct_alloc(int Nsite, int Nj, int Nq);
void lattice_struct_free(lattice_struct* las);

int lattice_struct_get_Nsite(const lattice_struct* las);
int lattice_struct_get_Nj(const lattice_struct* las);
int lattice_struct_get_Nq(const lattice_struct* las);
int lattice_struct_get_sigma0(const lattice_struct* las, int index);
int lattice_struct_get_sigmap(const lattice_struct* las, int index);
double lattice_struct_get_bond(const lattice_struct* las, int i_bond);
void lattice_struct_get_bond2index(int* index, const lattice_struct* las, int i_bond);

void lattice_struct_set_sigma0(lattice_struct* las, int index, int x);
void lattice_struct_set_sigmap(lattice_struct* las, int index, int x);
void lattice_struct_set_bond(lattice_struct* las, int i_bond, double x);
void lattice_struct_set_bond2index(lattice_struct* las, int i_bond, const int* index);


/* Note on operator sequence
** L        : The length of operator sequence
** n        : Number of non-identity operator
** sequence : storage the operator sequence
*/
typedef struct operator_sequence{
    int L;
    int n;
    int_array* sequence;
} operator_sequence;

operator_sequence* operator_sequence_alloc(int L);
void operator_sequence_free(operator_sequence* ops);

int operator_sequence_get_length(const operator_sequence* ops);
int operator_sequence_get_noo(const operator_sequence* ops);
int operator_sequence_get_sequence(const operator_sequence* ops, int p);
void operator_sequence_set_noo(operator_sequence* ops, int noo);
void operator_sequence_set_sequence(operator_sequence* ops, int p, int op);

void operator_sequence_init(operator_sequence* ops);


/* Note on link vertex list
** L       : The length of operator sequence
** Nsite   : The total number of spin
** lvl     : storage the link vertex list, which has size be 8L, leg+8*p
** v_first : first list 
** v_last  : last list
*/
typedef struct link_vertex{
    int L;
    int Nsite;
    int_array* lvl;
    int_array* v_first;
    int_array* v_last;
} link_vertex;

link_vertex* link_vertex_alloc(int L, int Nsite);
void link_vertex_free(link_vertex* lv);

int link_vertex_get_length(const link_vertex* lv);
int link_vertex_get_Nsite(const link_vertex* lv);
int link_vertex_get_list(const link_vertex* lv, int nu);
int link_vertex_get_first(const link_vertex* lv, int index);
int link_vertex_get_last(const link_vertex* lv, int index);
void link_vertex_set_list(link_vertex* lv, int nu, int x);
void link_vertex_set_first(link_vertex* lv, int index, int x);
void link_vertex_set_last(link_vertex* lv, int index, int x);

void link_vertex_init(link_vertex* lv);


/* Note one observable
** Nobs    : The total number of obsevable
** Nsample : The total number of sample for each block
** data    : storage of each sampling and observable [i_obs+Nobs*i_sample]
** init    : initialization of data. init[i]=1 means data[i] has been changed.
*/
typedef struct observable{
    int Nobs;
    int Nsample;
    double_array* data;
    int_array* init;
} observable;

observable* observable_alloc(int Nobs, int Nsample);
void observable_free(observable* obs);

int observable_get_Nobs(const observable* obs);
int observable_get_Nsample(const observable* obs);
double observable_get_data(const observable* obs, int i_obs, int i_sample);
int observable_get_init(const observable* obs, int i_obs, int i_sample);

void observable_set_data(observable* obs, int i_obs, int i_sample, double x);

void observable_init(observable* obs);
int observable_check_finish(const observable* obs);
#endif
