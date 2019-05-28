/* algorithm/measurement.h
*/

#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "low_level.h"
#include "mid_level.h"

void observable_result_stdout(const observable* obs, double beta);

void observable_result_fileout(const observable* obs, double beta, const char* filename);

void observable_measure_ms_2d(observable* obs, lattice_struct* las, const operator_sequence* ops, int i_sample);

void observable_measure_mz_2d(observable* obs, const lattice_struct* las, int i_sample);

void observable_measure_uniform_sus_2d(observable* obs, const lattice_struct* las, double beta, int i_sample);
#endif
