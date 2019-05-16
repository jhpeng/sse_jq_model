/* algorithm/lattice.h
*/
#include <gsl/gsl_rng>

#include "low_level.h"
#include "mid_level.h"

void lattice_struct_sync_sigmap(lattice_struct* las);

void lattice_struct_propagate_state(lattice_strcut* las, int op);
