#include "nrpyeos_hybrid.h"

void NRPyEOS_initialize_hybrid_functions() {
  ghl_hybrid_find_polytropic_index            = &NRPyEOS_find_polytropic_index;
  ghl_hybrid_find_polytropic_index_from_P     = &NRPyEOS_find_polytropic_index_from_P;
  ghl_hybrid_get_K_and_Gamma                  = &NRPyEOS_get_K_and_Gamma;
  ghl_hybrid_set_K_ppoly_and_eps_integ_consts = &NRPyEOS_set_K_ppoly_and_eps_integ_consts;
  ghl_hybrid_compute_P_cold                   = &NRPyEOS_compute_P_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold      = &NRPyEOS_compute_P_cold_and_eps_cold;
  ghl_hybrid_compute_entropy_function         = &NRPyEOS_compute_entropy_function;
  ghl_hybrid_compute_epsilon                  = &NRPyEOS_hybrid_compute_epsilon;
  ghl_hybrid_compute_rho_cold_from_P_cold     = &NRPyEOS_hybrid_compute_rho_cold_from_P_cold;
  ghl_compute_h_and_cs2                       = &NRPyEOS_hybrid_compute_enthalpy_and_cs2;
  ghl_hybrid_enforce_bounds__rho              = &NRPyEOS_hybrid_enforce_bounds__rho;
}