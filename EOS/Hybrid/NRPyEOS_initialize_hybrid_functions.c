#include "NRPyEOS_Hybrid.h"

void NRPyEOS_initialize_hybrid_functions(eos_parameters *restrict eos) {
  eos->hybrid_find_polytropic_index            = &NRPyEOS_find_polytropic_index;
  eos->hybrid_get_K_and_Gamma                  = &NRPyEOS_get_K_and_Gamma;
  eos->hybrid_set_K_ppoly_and_eps_integ_consts = &NRPyEOS_set_K_ppoly_and_eps_integ_consts;
  eos->hybrid_compute_P_cold                   = &NRPyEOS_compute_P_cold;
  eos->hybrid_compute_P_cold_and_eps_cold      = &NRPyEOS_compute_P_cold_and_eps_cold;
  eos->hybrid_compute_entropy_function         = &NRPyEOS_compute_entropy_function;
  eos->compute_h_and_cs2                       = &NRPyEOS_hybrid_compute_enthalpy_and_cs2;
}
