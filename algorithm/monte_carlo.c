/* algorithm/monte_carlo.c
*/
#include <stdio.h>
#include <gsl/gsl_rng.h>

#include "low_level.h"
#include "mid_level.h"
#include "lattice.h"

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

void diagonal_operator_update_J(operator_sequence* ops, lattice_struct* las, double beta, gsl_rng* rng){
    int p,i_bond;
    int L = operator_sequence_get_length(ops);
    int n = operator_sequence_get_noo(ops);
    int Nj = lattice_struct_get_Nj(las);
    int index[4],s1,s2,sp;
    double dis,J;

    lattice_struct_sync_sigmap(las);

    for(p=0;p<L;++p){
        sp = operator_sequence_get_sequence(ops,p);
        if(sp==-1){
            i_bond = (int)(gsl_rng_uniform_pos(rng)*Nj);
            lattice_struct_get_bond2index(index,las,i_bond);
            s1 = lattice_struct_get_sigmap(las,index[0]);
            s2 = lattice_struct_get_sigmap(las,index[1]);
            if(s1!=s2){
                dis = gsl_rng_uniform_pos(rng);
                J = lattice_struct_get_bond(las,i_bond);
                if(dis*2*(L-n)<beta*J*Nj){
                    operator_sequence_set_sequence(ops,p,i_bond*6);
                    n++;
                }
            }
        }
        else if(sp%6==0){
            i_bond = sp/6;
            dis = gsl_rng_uniform_pos(rng);
            J = lattice_struct_get_bond(las,i_bond);
            if(beta*Nj*J*dis<2*(L-n+1)){
                operator_sequence_set_sequence(ops,p,-1);
                n--;
            }
        }
        else lattice_struct_propagate_state(las,sp);
    }
    
    operator_sequence_set_noo(ops,n);
}

void diagonal_operator_update_Q(operator_sequence* ops, lattice_struct* las, double beta, gsl_rng* rng){
    int p,i_bond;
    int L = operator_sequence_get_length(ops);
    int n = operator_sequence_get_noo(ops);
    int Nj = lattice_struct_get_Nj(las);
    int Nq = lattice_struct_get_Nq(las);
    int index[4],s1,s2,s3,s4,sp;
    double dis,Q;

    lattice_struct_sync_sigmap(las);

    for(p=0;p<L;++p){
        sp = operator_sequence_get_sequence(ops,p);
        if(sp==-1){
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
                    operator_sequence_set_sequence(ops,p,i_bond*6+2);
                    n++;
                }
            }
        }
        else if(sp%6==2){
            i_bond = sp/6;
            dis = gsl_rng_uniform_pos(rng);
            Q = lattice_struct_get_bond(las,i_bond);
            if(beta*Nq*Q*dis<4*(L-n+1)){
                operator_sequence_set_sequence(ops,p,-1);
                n--;
            }
        }
        else lattice_struct_propagate_state(las,sp);
    }

    operator_sequence_set_noo(ops,n);
}

void construct_link_vertex_list(link_vertex* lv, const operator_sequence* ops, const lattice_struct* las){
    int i,p,sp,i_bond,type,nu1,nu0;
    int index[4];
    int L = operator_sequence_get_length(ops);
    int Nsite = lattice_struct_get_Nsite(las);

    link_vertex_init(lv);

    for(p=0;p<L;++p){
        sp = operator_sequence_get_sequence(ops,p);
        if(sp!=-1){
            i_bond = sp/6;
            type   = sp%6;
            lattice_struct_get_bond2index(index,las,i_bond);
            for(i=0;i<4;++i){
                if(index[i]!=-1){
                    nu1 = link_vertex_get_last(lv,index[i]);
                    nu0 = 8*p+i;
                    if(nu1!=-1){
                        link_vertex_set_list(lv,nu0,nu1);
                        link_vertex_set_list(lv,nu1,nu0);
                    }
                    else{
                        link_vertex_set_first(lv,index[i],nu0);
                    }
                    link_vertex_set_last(lv,index[i],nu0+4);
                }
            }
        }
    }

    for(i=0;i<Nsite;++i){
        nu0 = link_vertex_get_first(lv,i);
        if(nu0!=-1){
            nu1 = link_vertex_get_last(lv,i);
            link_vertex_set_list(lv,nu0,nu1);
            link_vertex_set_list(lv,nu1,nu0);
        }
    }
}

void loop_update(link_vertex* lv, gsl_rng* rng){
    int nu0,nup,nun,flip=-1;
    int L = link_vertex_get_length(lv);

    for(nu0=0;nu0<(L*8);nu0+=2){
        if(link_vertex_get_list(lv,nu0)>=0){
            nup = nu0;
            nun = nu0;
            if(gsl_rng_uniform_pos(rng)<0.5) flip=-1;
            else flip=-2;
            while(nun>=0){
                link_vertex_set_list(lv,nun,flip);
                nun = nun^1;
                nup = nun;
                nun = link_vertex_get_list(lv,nup);
                link_vertex_set_list(lv,nup,flip);
            }
        }
    }
}

void flip_bit_operator(operator_sequence* ops, lattice_struct* las, const link_vertex* lv, gsl_rng* rng){
    int nu,i_flip_rule,type,i_bond,sp,index,spin;
    int L = link_vertex_get_length(lv);
    int Nsite = lattice_struct_get_Nsite(las);

    for(nu=0;nu<(8*L);nu+=2){
        if(link_vertex_get_list(lv,nu)==-2){
            sp = operator_sequence_get_sequence(ops,nu/8);
            type = sp%6;
            i_bond = sp/6;
            i_flip_rule = ((nu%8)/2)*6+type;
            type = flip_rule[i_flip_rule];
            sp = i_bond*6+type;
            operator_sequence_set_sequence(ops,nu/8,sp);
        }
    }

    for(index=0;index<Nsite;++index){
        nu = link_vertex_get_last(lv,index);
        if(nu==-1){
            if(gsl_rng_uniform_pos(rng)<0.5){
                spin = lattice_struct_get_sigma0(las,index);
                lattice_struct_set_sigma0(las,index,-1*spin);
            }
        }
        else{
            if(link_vertex_get_list(lv,nu)==-2){
                spin = lattice_struct_get_sigma0(las,index);
                lattice_struct_set_sigma0(las,index,-1*spin);
            }
        }
    }
}

void adjust_cutoff(operator_sequence** ops, link_vertex** lv, double buffer){
    int L = operator_sequence_get_length(*ops);
    int n = operator_sequence_get_noo(*ops);

    if(n*buffer>L){
        int sp,L_pre=L;
        int Nsite = link_vertex_get_Nsite(*lv);
        L = (int)(n*buffer)+1;

        operator_sequence* ops_new = operator_sequence_alloc(L);
        link_vertex* lv_new = link_vertex_alloc(L,Nsite);
        operator_sequence_init(ops_new);

        for(int p=0;p<L_pre;++p){
            sp = operator_sequence_get_sequence(*ops,p);
            operator_sequence_set_sequence(ops_new,p,sp);
        }
        operator_sequence_set_noo(ops_new,n);
        operator_sequence_free(*ops);
        link_vertex_free(*lv);
        *ops = ops_new;
        *lv = lv_new;
    }
}

int operator_sequence_count_Q(operator_sequence* ops){
    int sp,n=0;
    int L = operator_sequence_get_length(ops);

    for(int p=0;p<L;++p){
        sp = operator_sequence_get_sequence(ops,p);
        if(sp%6>1) n++;
    }

    return n;
}

#if 1
int main(int argc, char **argv)
{
    int Nx=16;
    int Ny=16;
    double J=atof(argv[1]);
    double beta=32,buffer=1.3;
    int L = 30;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);    
    lattice_struct* las = lattice_struct_create_model_plaquette_2d(Nx,Ny,J,rng);
    operator_sequence* ops = operator_sequence_alloc(L);
    operator_sequence_init(ops);
    
    int Nsite = lattice_struct_get_Nsite(las);

    link_vertex* lv = link_vertex_alloc(L,Nsite);

    for(int i=0;i<2000;++i){
        diagonal_operator_update_Q(ops,las,beta,rng);
        diagonal_operator_update_J(ops,las,beta,rng);
        construct_link_vertex_list(lv,ops,las);
        loop_update(lv,rng);
        flip_bit_operator(ops,las,lv,rng);
        L = operator_sequence_get_length(ops);
        /*for(int j=0;j<L;++j){
            printf("%d ",operator_sequence_get_sequence(ops,j)%6);
        }
        printf("\n");*/
        adjust_cutoff(&ops,&lv,buffer);
        L = operator_sequence_get_length(ops);
    }
    int check_p = lattice_struct_check_propagate_state(las,ops);
    int n_q =  operator_sequence_count_Q(ops);
    printf("%d %d %d %d\n",L,operator_sequence_get_noo(ops),n_q,check_p);

    return 0;
}
#endif
