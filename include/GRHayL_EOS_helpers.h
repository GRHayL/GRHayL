#ifndef GRHAYL_EOS_HELPERS_H_
#define GRHAYL_EOS_HELPERS_H_

#include "GRHayL.h"

int GRHayL_find_polytropic_index(
      const eos_parameters *restrict eos,
      const double rho_in);

void GRHayL_get_K_and_Gamma(
      const eos_parameters *restrict eos,
      const double rho_in,
      double *restrict K,
      double *restrict Gamma);

void GRHayL_set_K_ppoly_and_eps_integ_consts(eos_parameters *restrict eos);

#endif // GRHAYL_EOS_HELPERS_H_
