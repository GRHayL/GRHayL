#ifndef EOS_HYBRID_HEADER_H
#define EOS_HYBRID_HEADER_H

#include "EOS_header.h"

int find_polytropic_K_and_Gamma_index(const eos_parameters *restrict eos, const double rho_in);

void compute_entropy_function( const eos_parameters *restrict eos, const double rho,
                               const double P, double *restrict S );

#endif //EOS_HYBRID_HEADER_H
