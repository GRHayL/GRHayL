#include "EOS_hybrid_header.h"

/* Function    : compute_P_cold_and_eps_cold()
 * Authors     : Leo Werneck
 * Description : Computes P_cold and eps_cold.
 * Dependencies: initialize_igm_eos_parameters_from_input()
 *             : find_polytropic_K_and_Gamma_index()
 *
 * Inputs      : P_cold         - cold pressure
 *             : eps_cold       - cold specific internal energy
 *             : eos            - an initialized eos_parameters struct
 *                                with data for the EOS of the simulation
 *
 * Outputs     : P_cold         - cold pressure (supports SPEOS and PPEOS)
 *             : eps_cold       - cold specific internal energy (supports SPEOS and PPEOS)
 *
 *             SPEOS: Single-Polytrope Equation of State
 *             PPEOS: Piecewise Polytrope Equation of State
 */

void compute_P_cold_and_eps_cold(const eos_parameters *restrict eos, const double rho_in,
                              double *restrict P_cold_ptr, double *restrict eps_cold_ptr) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
  if(rho_in==0) {
    *P_cold_ptr   = 0.0;
    *eps_cold_ptr = 0.0;
    return;
  }

  /*  --------------------------------------------------
   * | Single and Piecewise Polytropic EOS modification |
   *  --------------------------------------------------
   *
   * We now begin our modifications to this function so that
   * it supports both single and piecewise polytropic equations
   * of state.
   *
   * The modifications below currently assume that the user
   * has called the recently added function
   *
   * - initialize_igm_eos_parameters_from_input()
   *
   * *before* this function is called. We can add some feature
   * to check this automatically as well, but we'll keep that as
   * a TODO/FIXME for now.
   */

  /* First, we compute the pressure, which in the case of a
   * piecewise polytropic EOS is given by
   *
   *           /   P_{1}      /    K_{1} * rho^(Gamma_{1})       ,      rho_{0} <= rho < rho_{1}
   *           |    ...       |            ...
   * P(rho) = <    P_{j}   = <     K_{j} * rho^(Gamma_{j})       ,    rho_{j-1} <= rho < rho_{j}
   *           |    ...       |            ...
   *           \ P_{neos-2}   \ K_{neos-2} * rho^(Gamma_{neos-2}), rho_{neos-3} <= rho < rho_{neos-2}
   *
   * The index j is determined by the find_polytropic_K_and_Gamma_index() function.
   */

  // Set up useful auxiliary variables
  int polytropic_index = find_polytropic_index(eos, rho_in);
  double K_ppoly       = eos->K_ppoly[polytropic_index];
  double Gamma_ppoly   = eos->Gamma_ppoly[polytropic_index];
  double eps_integ_const = eos->eps_integ_const[polytropic_index];

  // Then compute P_{cold}
  double P_cold = K_ppoly*pow(rho_in, Gamma_ppoly);

  /* Then we compute the cold component of the specific internal energy,
   * which in the case of a piecewise polytropic EOS is given by (neos -> N)
   *
   *             /   P_{1}/(rho*(Gamma_{1}-1))   + C_{1}  ,   rho_{0} <= rho < rho_{1}
   *             |                     ...
   * eps(rho) = <    P_{j}/(rho*(Gamma_{j}-1))   + C_{j}  , rho_{j-1} <= rho < rho_{j}
   *             |                     ...
   *             \ P_{N-2}/(rho*(Gamma_{N-2}-1)) + C_{N-2}, rho_{N-3} <= rho < rho_{N-2}
   */
  *eps_cold_ptr = P_cold/(rho_in*(Gamma_ppoly-1.0)) + eps_integ_const;

  *P_cold_ptr = P_cold;
}
