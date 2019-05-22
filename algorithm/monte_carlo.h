/* algorithm/monte_carlo.h
*/
#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include "low_level.h"
#include "mid_level.h"

void diagonal_operator_update_J(operator_sequence* ops, lattice_struct* las, double beta, gsl_rng* rng);

void diagonal_operator_update_Q(operator_sequence* ops, lattice_struct* las, double beta, gsl_rng* rng);

void construct_link_vertex_list(link_vertex* lv, const operator_sequence* ops, const lattice_struct* las);

void loop_update(link_vertex* lv, gsl_rng* rng);

void flip_bit_operator(operator_sequence* ops, lattice_struct* las, const link_vertex* lv, gsl_rng* rng);

void adjust_cutoff(operator_sequence** ops, link_vertex** lv, double buffer);

int operator_sequence_count_Q(operator_sequence* ops);

#endif
