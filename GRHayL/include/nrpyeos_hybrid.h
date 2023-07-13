#ifndef NRPYEOS_HYBRID_H_
#define NRPYEOS_HYBRID_H_

#include "ghl.h"

#ifdef __cplusplus
extern "C" {
#endif

int NRPyEOS_find_polytropic_index(
      const ghl_eos_parameters *restrict eos,
      const double rho_in);

void NRPyEOS_get_K_and_Gamma(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict K,
      double *restrict Gamma);

void NRPyEOS_set_K_ppoly_and_eps_integ_consts(ghl_eos_parameters *restrict eos);

void NRPyEOS_compute_P_cold(
      const struct ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr);

void NRPyEOS_compute_P_cold_and_eps_cold(
      const struct ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr,
      double *restrict eps_cold_ptr);

void NRPyEOS_compute_entropy_function(
      const struct ghl_eos_parameters *restrict eos,
      const double rho,
      const double P,
      double *restrict S );

void NRPyEOS_initialize_hybrid_functions(ghl_eos_parameters *restrict eos);

void NRPyEOS_hybrid_compute_enthalpy_and_cs2(
      ghl_eos_parameters const *restrict eos,
      ghl_primitive_quantities const *restrict prims,
      double *restrict enthalpy_ptr,
      double *restrict cs2_ptr );

#ifdef __cplusplus
}
#endif

#endif // NRPYEOS_HYBRID_H_
