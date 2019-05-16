/* algorithm/monte_carlo.c
*/
#include <gsl/gsl_rng.h>

#include "low_level.h"
#include "mid_level.h"
#include "lattice.h"

void diagonal_operator_update_J(operator_sequence* ops, const lattice_struct* las, double beta, gsl_rng* rng){
    int p,i_bond;
    int length = operator_sequence_get_length(ops);
    int n      = operator_sequence_get_noo(ops);
    int Nj = lattice_struct_get_Nj(las);
    int index[4],s1,s2,op;
    double dis,J;

    lattice_struct_sync_sigmap(las);

    for(p=0;p<length;++p){
        op = operator_sequence_get_sequence(ops,p);
        if(op==-1){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*Nj);
            lattice_struct_get_bond2index(index,las,i_bond);
            s1 = lattice_struct_get_sigmap(las,index[0]);
            s2 = lattice_struct_get_sigmap(las,index[1]);
            if(s1!=s2){
                dis = gsl_rng_uniform_pos(rng);
                J = lattice_struct_get_bond(las,i_bond);
                if(dis*2*(L-n)<beta*J*Nj){
                    operator_sequence_set_sequecne(ops,p,i_bond*6);
                    n++;
                }
            }
        }
        else if(op%6==0){
            i_bond = op/6;
            dis = gsl_rng_uniform_pos(rng);
            J = lattice_struct_get_bond(las,i_bond);
            if(beta*Nj*J*dis<2*(L-n+1)){
                operator_sequence_bond(ops,p,-1);
                n--;
            }
        }
        else lattice_struct_propagate_state(las,op);
    }
}

void diagonal_operator_update_Q(operator_sequence* ops, const lattice_struct* las, double beta, gsl_rng* rng){
    int p,i_bond;
    int length = operator_sequence_get_length(ops);
    int n      = operator_sequence_get_noo(ops);
    int Nj = lattice_struct_get_Nj(las);
    int Nq = lattice_struct_get_Nq(las);
    int index[4],s1,s2,s3,s4,op;
    double dis,Q;

    lattice_struct_sync_sigmap(las);

    for(p=0;p<length;++p){
        op = operator_sequence_get_sequence(ops,p);
        if(op==-1){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*Nq)+Nj;
            lattice_struct_get_bond2index(index,las,i_bond);
            s1 = lattice_struct_get_sigmap(las,index[0]);
            s2 = lattice_struct_get_sigmap(las,index[1]);
            s3 = lattice_struct_get_sigmap(las,index[2]);
            s4 = lattice_struct_get_sigmap(las,index[3]);
            if(s1!=s2 && s3!=s4){
                dis = gsl_rng_uniform_pos(rng);
                Q = lattice_struct_get_bond(las,i_bond);
                if(dis*4*(L-n)<beta*Q*Nq){
                    operator_sequence_set_sequecne(ops,p,i_bond*6+2);
                    n++;
                }
            }
        }
        else if(op%6==2){
            i_bond = op/6;
            dis = gsl_rng_uniform_pos(rng);
            Q = lattice_struct_get_bond(las,i_bond);
            if(beta*Nq*Q*dis<4*(L-n+1)){
                operator_sequence_bond(ops,p,-1);
                n--;
            }
        }
        else lattice_struct_propagate_state(las,op);
    }
}
