#include "nrpyeos_hybrid.h"

/* Function    : NRPyEOS_hybrid_compute_enthalpy_and_cs2()
 * Description : Computes the enthalpy and the sound speed squared.
 * Dependencies: initialize_eos_parameters_from_input()
 *             : find_polytropic_K_and_Gamma_index()
 *
 * Inputs      : eos            -
 *             : P_cold         - cold pressure
 *             : eps_cold       - cold specific internal energy
 *             : eos            - an initialized ghl_eos_parameters struct
 *                                with data for the EOS of the simulation
 *
 * Outputs     : enthalpy_ptr   - enthalpy
 *             : cs2_ptr        - sound speed squared
 *
 *             SPEOS: Single-Polytrope Equation of State
 *             PPEOS: Piecewise Polytrope Equation of State
 */
void NRPyEOS_hybrid_compute_enthalpy_and_cs2(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      double *restrict enthalpy_ptr,
      double *restrict cs2_ptr) {

  // Step 1: Unpack the prims struct
  const double rho = prims->rho;
  const double P   = prims->press;

  // Step 2: Compute P_cold and eps_cold
  double P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, rho, &P_cold, &eps_cold);

  // Step 3: Set Gamma cold
  const int polytropic_index = ghl_hybrid_find_polytropic_index(eos, rho);
  const double Gamma = eos->Gamma_ppoly[polytropic_index];

  // Step 4: Compute the derivative of cold pressure w.r.t. density,
  //   dP/drho = Gamma K rho^(Gamma-1)
  //           = Gamma/rho (K rho^Gamma)
  //           = Gamma P_cold/rho .
  const double dPcold_drho = Gamma * P_cold / rho;

  // Step 5: Compute eps_thermal
  const double eps_th = (P-P_cold)/( (eos->Gamma_th-1)*rho );

  // Step 6: Compute eps
  const double eps = eps_cold + eps_th;

  // Step 7: Compute the enthalpy
  const double h = 1 + eps + P/rho;

  // Step 8: Compute cs2
  const double cs2 = (dPcold_drho + eos->Gamma_th*(eos->Gamma_th-1)*eps_th)/h;

  // Step 9: Set the output
  *enthalpy_ptr = h;
  *cs2_ptr      = cs2;
}
