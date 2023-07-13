#include "nrpyeos_hybrid.h"

void NRPyEOS_initialize_hybrid_functions(ghl_eos_parameters *restrict eos) {
  ghl_hybrid_find_polytropic_index            = &NRPyEOS_find_polytropic_index;
  ghl_hybrid_get_K_and_Gamma                  = &NRPyEOS_get_K_and_Gamma;
  ghl_hybrid_set_K_ppoly_and_eps_integ_consts = &NRPyEOS_set_K_ppoly_and_eps_integ_consts;
  ghl_hybrid_compute_P_cold                   = &NRPyEOS_compute_P_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold      = &NRPyEOS_compute_P_cold_and_eps_cold;
  ghl_hybrid_compute_entropy_function         = &NRPyEOS_compute_entropy_function;
  ghl_compute_h_and_cs2                       = &NRPyEOS_hybrid_compute_enthalpy_and_cs2;
}
