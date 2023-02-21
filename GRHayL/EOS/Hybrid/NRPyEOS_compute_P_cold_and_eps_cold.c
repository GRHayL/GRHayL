#include "NRPyEOS_Hybrid.h"

/* Function    : compute_P_cold_and_eps_cold()
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

void NRPyEOS_compute_P_cold_and_eps_cold(
      const eos_parameters *restrict eos,
      const double rho_in,
      double *restrict P_cold_ptr,
      double *restrict eps_cold_ptr) {
  // This code handles equations of state of the form defined
  // in Eqs 13-16 in http://arxiv.org/pdf/0802.0200.pdf
  //
  // Set up useful auxiliary variables
  int polytropic_index   = eos->hybrid_find_polytropic_index(eos, rho_in);
  double K_ppoly         = eos->K_ppoly[polytropic_index];
  double Gamma_ppoly     = eos->Gamma_ppoly[polytropic_index];
  double eps_integ_const = eos->eps_integ_const[polytropic_index];

  // Then compute P_{cold}
  double P_cold = K_ppoly*pow(rho_in, Gamma_ppoly);

  *P_cold_ptr   = P_cold;
  *eps_cold_ptr = P_cold/(rho_in*(Gamma_ppoly-1.0)) + eps_integ_const;
}
