#include "nrpyeos_hybrid.h"

/*
 * Function     : compute_epsilon()
 * Description  : Computes epsilon from rho and P
 * Documentation:
*/

double NRPyEOS_hybrid_compute_epsilon(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double press) {
  double P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, rho, &P_cold, &eps_cold);
  return eps_cold + (press-P_cold)/(eos->Gamma_th-1.0)/rho;
}
