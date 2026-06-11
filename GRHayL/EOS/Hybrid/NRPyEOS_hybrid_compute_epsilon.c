#include "ghl_nrpyeos_hybrid.h"

/**
 * @ingroup hyb_eos
 * @brief Computes epsilon from rho and P; usually aliased as ghl_hybrid_compute_epsilon
 *
 * @details
 * This function computes the specific internal energy using
 *
 * \f$ \epsilon = \epsilon_\mathrm{cold} + \frac{P - P_\mathrm{cold}}{\rho(\Gamma_\mathrm{th} - 1)} \f$
 *
 * @param[in] eos:   pointer to ghl_eos_parameters struct
 *
 * @param[in] rho:   density value
 *
 * @param[in] press: pressure value
 *
 * @returns the specific internal energy
 */
double NRPyEOS_hybrid_compute_epsilon(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double press) {
  double P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, rho, &P_cold, &eps_cold);
  return eps_cold + (press-P_cold)/(eos->Gamma_th-1.0)/rho;
}
