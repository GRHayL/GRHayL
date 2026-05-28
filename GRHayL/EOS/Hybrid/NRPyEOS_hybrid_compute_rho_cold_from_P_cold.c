#include "ghl_nrpyeos_hybrid.h"

/**
 * @ingroup hyb_eos
 * @brief Computes rho_cold; usually aliased as ghl_hybrid_compute_rho_cold_from_P_cold
 *
 * @details
 * This code handles equations of state of the form defined
 * in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
 * @todo doxygen citation
 *
 * Computes the density using
 *
 * \f$ \rho = \left( \frac{P}{K_i} \right)^\frac{1}{\Gamma_i} \f$
 *
 * @param[in] eos:  pointer to ghl_eos_parameters struct
 *
 * @param[in] P_in: density value
 *
 * @returns the cold density
 */
double NRPyEOS_hybrid_compute_rho_cold_from_P_cold(
      const ghl_eos_parameters *restrict eos,
      const double P_in) {

  const int polytropic_index = ghl_hybrid_find_polytropic_index_from_P(eos, P_in);
  const double K_ppoly       = eos->K_ppoly[polytropic_index];
  const double Gamma_ppoly   = eos->Gamma_ppoly[polytropic_index];
  double rho = pow(P_in / K_ppoly, 1.0 / Gamma_ppoly);
  NRPyEOS_hybrid_enforce_bounds__rho(eos, &rho);
  return rho;
}
