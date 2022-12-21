#ifndef EOS_HYBRID_H
#define EOS_HYBRID_H

#include "GRHayL.h"

void compute_P_cold_and_eps_cold(
             const eos_parameters *restrict eos, const double rho_in,
             double *restrict P_cold_ptr, double *restrict eps_cold_ptr);

void compute_P_cold(const eos_parameters *restrict eos, const double rho_in, double *restrict P_cold_ptr);

void compute_entropy_function( const eos_parameters *restrict eos, const double rho,
                               const double P, double *restrict S );

#endif //EOS_HYBRID_H
