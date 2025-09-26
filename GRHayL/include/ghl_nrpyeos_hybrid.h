#ifndef NRPYEOS_HYBRID_H_
#define NRPYEOS_HYBRID_H_

#include "ghl.h"

#ifdef __cplusplus
extern "C" {
#endif

GHL_DEVICE
int NRPyEOS_find_polytropic_index(
      const ghl_eos_parameters *restrict eos,
      const double rho_in);

GHL_DEVICE
int NRPyEOS_find_polytropic_index_from_P(
      const ghl_eos_parameters *restrict eos,
      const double P_in);

GHL_DEVICE
void NRPyEOS_get_K_and_Gamma(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict K,
      double *restrict Gamma);

GHL_DEVICE
void NRPyEOS_set_K_ppoly_and_eps_integ_consts(
      ghl_eos_parameters *restrict eos);

GHL_DEVICE
void NRPyEOS_compute_P_cold(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr);

GHL_DEVICE
void NRPyEOS_compute_P_cold_and_eps_cold(
      const ghl_eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr,
      double *restrict eps_cold_ptr);

GHL_DEVICE
double NRPyEOS_hybrid_compute_rho_cold_from_P_cold(
      const ghl_eos_parameters *restrict eos,
      const double P_in);

GHL_DEVICE
double NRPyEOS_hybrid_compute_epsilon(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double P);

GHL_DEVICE
double NRPyEOS_compute_entropy_function(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double P);

GHL_DEVICE
void NRPyEOS_initialize_hybrid_functions();

GHL_DEVICE
bool NRPyEOS_hybrid_enforce_bounds__rho(const ghl_eos_parameters *restrict eos, double *restrict rho);

GHL_DEVICE
ghl_error_codes_t NRPyEOS_hybrid_compute_enthalpy_and_cs2(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      double *restrict enthalpy_ptr,
      double *restrict cs2_ptr);

#ifdef __cplusplus
}
#endif

#endif // NRPYEOS_HYBRID_H_
