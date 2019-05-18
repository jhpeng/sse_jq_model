/* algorithm/lattice.h
*/
#include <gsl/gsl_rng.h>

#include "low_level.h"
#include "mid_level.h"

void lattice_struct_sync_sigmap(lattice_struct* las);

void lattice_struct_propagate_state(lattice_struct* las, int sp);

int lattice_struct_check_propagate_state(lattice_struct* las, operator_sequence* ops);

lattice_struct* lattice_struct_create_model_plaquette_2d(int Nx, int Ny, double J, gsl_rng* rng);
