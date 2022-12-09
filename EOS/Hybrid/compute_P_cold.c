#include "EOS_hybrid.h"

/* Function    : compute_P_cold()
 * Authors     : Leo Werneck & Samuel Cupp
 * Description : Computes P_cold.
 *
 * Inputs      : P_cold           - cold pressure
 *             : eps_cold         - cold specific internal energy
 *             : eos              - a struct containing the following
 *                                  relevant quantities:
 *             : neos             - number of polytropic EOSs used
 *             : rho_ppoly_tab         - array of rho values that determine
 *                                  the polytropic EOS to be used.
 *             : Gamma_ppoly_tab       - array of Gamma_cold values to be
 *                                  used in each polytropic EOS.
 *             : K_ppoly_tab           - array of K_ppoly_tab values to be used
 *                                  in each polytropic EOS.
 *             : eps_integ_const  - array of C_{j} values, which are the
 *                                  integration constants that arrise when
 *                                  determining eps_{cold} for a piecewise
 *                                  polytropic EOS.
 *
 * Outputs     : P_cold           - cold pressure (supports SPEOS and PPEOS)
 *
 *             SPEOS: Single-Polytrope Equation of State
 *             PPEOS: Piecewise Polytrope Equation of State
 */

void compute_P_cold(const eos_parameters *restrict eos, const double rho_in, double *restrict P_cold_ptr) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
  if(rho_in==0) {
    *P_cold_ptr   = 0.0;
    return;
  }

  /* In the case of a piecewise polytropic EOS is given by
   *
   *           /   P_{1}      /    K_{1} * rho^(Gamma_{1})       ,      rho_{0} <= rho < rho_{1}
   *           |    ...       |            ...
   * P(rho) = <    P_{j}   = <     K_{j} * rho^(Gamma_{j})       ,    rho_{j-1} <= rho < rho_{j}
   *           |    ...       |            ...
   *           \ P_{neos-2}   \ K_{neos-2} * rho^(Gamma_{neos-2}), rho_{neos-3} <= rho < rho_{neos-2}
   *
   * The index j is determined by the find_polytropic_index() function.
   */

  // Set up useful auxiliary variables
  double K_ppoly;
  double Gamma_ppoly;
  get_K_and_Gamma(eos, rho_in, &K_ppoly, &Gamma_ppoly);

  // Then compute P_{cold}
  *P_cold_ptr = K_ppoly*pow(rho_in,Gamma_ppoly);
}
