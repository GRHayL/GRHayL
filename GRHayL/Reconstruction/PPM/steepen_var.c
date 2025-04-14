#include "ghl_reconstruction.h"

/*
 * Function     : ghl_steepen_var()
 * Description  : applies a steepening algorithm to the
 *                reconstructed variable var
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_steepen_var
*/

void ghl_steepen_var(
      const ghl_parameters *restrict params,
      const double rho[5],
      const double pressure[5],
      const double Gamma_eff,
      double *restrict rhor,
      double *restrict rhol) {

  // Next compute centered differences d RHOB and d^2 RHOB
  // d1rho_b is multiplied by 2 to eliminate common factors in expressions
  const double d1rho_b     = rho[PLUS_1] - rho[MINUS1];
  const double d2rho_b_m1  = rho[PLUS_0] - 2.0*rho[MINUS1] + rho[MINUS2];
  const double d2rho_b_p1  = rho[PLUS_2] - 2.0*rho[PLUS_1] + rho[PLUS_0];
  const double rho_b_min = MIN(rho[PLUS_1], rho[MINUS1]);

  // Gamma_eff = (partial P / partial rho0)_s /(P/rho0)
  const bool contact_discontinuity_check =
                 Gamma_eff * params->ppm_shock_k0 * fabs(d1rho_b) * MIN(pressure[PLUS_1], pressure[MINUS1])
                 >= fabs(pressure[PLUS_1]-pressure[MINUS1])*rho_b_min;
  const bool second_deriv_check = d2rho_b_p1*d2rho_b_m1 <= 0.0;
  const bool relative_change_check = fabs(d1rho_b) >= params->ppm_shock_epsilon*rho_b_min;

  if(contact_discontinuity_check && second_deriv_check && relative_change_check) {

    const double slope_limited_drho_m1 = ghl_slope_limit(rho[MINUS1] - rho[MINUS2], rho[PLUS_0] - rho[MINUS1]);
    const double slope_limited_drho_p1 = ghl_slope_limit(rho[PLUS_1] - rho[PLUS_0], rho[PLUS_2] - rho[PLUS_1]);

    double eta_tilde = 0.0;
    if (fabs(d1rho_b) > 0.0) {
      eta_tilde = -(1.0/6.0)*(d2rho_b_p1 - d2rho_b_m1)/d1rho_b;
    }
    const double eta = MAX(0.0, MIN(params->ppm_shock_eta1*(eta_tilde - params->ppm_shock_eta2), 1.0));

    // Next compute Urp1 and Ul for RHOB, using the MC prescription:
    // Ur_p1 = U_p1   - 0.5*slope_lim_dU_p1
    const double rho_br_mc_p1 = rho[PLUS_1] - 0.5*slope_limited_drho_p1;
    // Ul = U_m1 + 0.5*slope_lim_dU_m1
    // Based on this line of code, Ur[index] = a_j - \delta_m a_j / 2. (cf. Eq. 65 in Marti & Muller's "PPM Method for 1D Relativistic Hydro." paper)
    //    So: Ur[indexp1] = a_{j+1} - \delta_m a_{j+1} / 2. This is why we have rho_br_mc[indexp1]
    const double rho_bl_mc    = rho[MINUS1] + 0.5*slope_limited_drho_m1;

    *rhol = (*rhol)*(1.0 - eta) + rho_bl_mc*eta;
    *rhor = (*rhor)*(1.0 - eta) + rho_br_mc_p1*eta;
  }
}
