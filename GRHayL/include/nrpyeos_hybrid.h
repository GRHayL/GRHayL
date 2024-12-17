#ifndef NRPYEOS_HYBRID_H_
#define NRPYEOS_HYBRID_H_

#include "ghl.h"

#ifdef __cplusplus
extern "C" {
#endif

int NRPyEOS_find_polytropic_index(
      const ghl_eos_parameters *restrict eos,
      const double rho_in);

int NRPyEOS_find_polytropic_index_from_h(
      const ghl_eos_parameters *restrict eos,
      const double h_in);

int NRPyEOS_find_polytropic_index_from_P(
      const ghl_eos_parameters *restrict eos,
      const double P_in);

void NRPyEOS_get_K_and_Gamma(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict K,
      double *restrict Gamma);

void NRPyEOS_set_K_ppoly_and_eps_integ_consts(
      ghl_eos_parameters *restrict eos);

void NRPyEOS_compute_P_cold(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr);

void NRPyEOS_compute_P_cold_and_eps_cold(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr,
      double *restrict eps_cold_ptr);

double NRPyEOS_hybrid_compute_rho_cold_from_h(
      const ghl_eos_parameters *restrict eos,
      const double h_in);

double NRPyEOS_hybrid_compute_rho_cold_from_P_cold(
      const ghl_eos_parameters *restrict eos,
      const double P_in);

double NRPyEOS_hybrid_compute_epsilon(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double P);

double NRPyEOS_compute_entropy_function(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double P);

void NRPyEOS_initialize_hybrid_functions();

ghl_error_codes_t NRPyEOS_hybrid_compute_enthalpy_and_cs2(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      double *restrict enthalpy_ptr,
      double *restrict cs2_ptr);

#ifdef __cplusplus
}
#endif

#endif // NRPYEOS_HYBRID_H_
