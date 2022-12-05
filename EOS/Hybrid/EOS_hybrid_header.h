#ifndef EOS_HYBRID_HEADER_H
#define EOS_HYBRID_HEADER_H

#include "EOS_header.h"

int find_polytropic_index(const eos_parameters *restrict eos, const double rho_in);

void get_K_and_Gamma(const eos_parameters *restrict eos, const double rho_in, double *restrict K, double *restrict Gamma);

void compute_P_cold_and_eps_cold(
             const eos_parameters *restrict eos, const double rho_in,
             double *restrict P_cold_ptr, double *restrict eps_cold_ptr);

void compute_P_cold(const eos_parameters *restrict eos, const double rho_in, double *restrict P_cold_ptr);

void compute_entropy_function( const eos_parameters *restrict eos, const double rho,
                               const double P, double *restrict S );

void setup_K_ppoly_and_eps_integ_consts(eos_parameters *restrict eos);

#endif //EOS_HYBRID_HEADER_H
