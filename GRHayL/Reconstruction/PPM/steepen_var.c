#include "ghl_reconstruction.h"

/**
 * @ingroup ppm_internal
 * @brief Applies a steepening algorithm to a variable.
 *
 * This function takes the right and left face values of a variable and applies
 * a steepening algorithm. This is normally used for the density, as the density
 * profile should be narrowed in the presence of the contact discontinuity
 * (see e.g. Appendix I of \cite Marti_1996. To determine if there is a shock
 * that requires steepening, several conditions must be met. The primary
 * equation is eq. 63 from \cite Marti_1996 :
 *                                            
 * \f[                                        
 * \Gamma K_0 \frac{|\rho_{i+1} - \rho_{i-1}|}{\mathrm{min}\left(\rho_{i+1},
 *                                                       \rho_{i-1}\right)}
 * \ge \frac{|P_{i+1} - P_{i-1}|}{\mathrm{min}\left(P_{i+1}, P_{i-1}\right)}
 * \f]
 *
 * where \f$ \Gamma \f$ is the Gamma\_eff input variable and \f$ K_0 \f$ is
 * a constant. We also check that the second derivatives of the adjacent
 * points are opposite signs
 *
 * \f[
 * \left( \rho_{i+2} + \rho_{i} - 2\rho_{i+1} \right)
 * \left( \rho_{i-2} + \rho_{i} - 2\rho_{i-1} \right) <= 0
 * \f]
 *
 * Finally, we check if the first derivative is a significant fraction of
 * the actual density
 *
 * \f[
 * |\rho_{i+1} - \rho_{i-1}| \ge \epsilon_0\ \mathrm{min}\left(\rho_{i+1},
 *                                                       \rho_{i-1}\right)
 * \f]
 *
 * If all three conditions are met, then we apply the steepening procedure.
 * We start by computing the MC reconstruction with the limited slope.
 * We first compute the slopes using the same prescription as in 
 * ghl_ppm_compute_for_cell:
 *
 * \f[
 * \begin{aligned}
 * s_- &= \mathrm{slope\_limit} \left( U_{i-1} U_{i-2}, U_{i} - U_{i-1} \right) \\
 * s_+ &= \mathrm{slope\_limit} \left( U_{i+1} U_{i}, U_{i+2} - U_{i+1} \right)
 * \end{aligned}
 * \f]
 *
 * Then, we compute the MC-reconstructed face values
 *
 * \f[
 * \begin{aligned}
 * \rho^\mathrm{MC}_R &= \rho_{i+1} - \frac{s_+}{2} \\
 * \rho^\mathrm{MC}_L &= \rho_{i-1} + \frac{s_-}{2}
 * \end{aligned}
 * \f]
 *
 * We also compute \f$ \eta \f$ using
 *
 * \f[
 * \eta = \mathrm{max}\left[0,\ \mathrm{min}\left( 1,\ \eta_1 (\tilde{\eta} - \eta_2) \right) \right]
 * \f]
 *
 * where
 *
 * \f[
 * \tilde{\eta} = - \frac{\partial^2 \rho_+ - \partial^2 \rho_-}{6\partial \rho}
 * \f]
 *
 * and
 *
 * \f[
 * \begin{aligned}
 * \partial \rho     &= \rho_{i+1} - \rho_{i-1} \\
 * \partial^2 \rho_+ &= \rho_{i+2} + \rho_{i} - 2\rho_{i+1} \\
 * \partial^2 \rho_- &= \rho_{i-2} + \rho_{i} - 2\rho_{i-1}
 * \end{aligned}
 * \f]
 *
 * Using these, the steepening procedure is applied to \f$ \rho_R \f$ and \f$ \rho_L \f$:
 *
 * \f[
 * \begin{aligned}
 * \rho_R &= (1 - \eta) \rho_R + \eta\ \rho^\mathrm{MC}_R \\
 * \rho_L &= (1 - \eta) \rho_L + \eta\ \rho^\mathrm{MC}_L \\
 * \end{aligned}
 * \f]
 *
 * The constants \f$ K_0 \f$, \f$ \epsilon_0 \f$, \f$ \eta_1 \f$, and
 * \f$ \eta_2 \f$ are all changeable by the user, and the @ref ghl_initialize_params
 * sets these to reasonable default values based on the original
 * Colella & Woodward paper \cite Colella_1984.
 *
 * @param[in] params:    pointer to a ghl_parameters struct
 *
 * @param[in] rho:       1D array containing stencil of the reconstructed variable
 *
 * @param[in] pressure:  1D array containing stencil of the pressure
 *
 * @param[in] Gamma_eff: value of \f$ \Gamma \f$ to use
 *
 * @param[in,out] rhor:  pointer to reconstructed value of the right face
 *
 * @param[in,out] rhol:  pointer to reconstructed value of the left face
 *
 * @returns void
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
    const double eta = ghl_clamp(params->ppm_shock_eta1*(eta_tilde - params->ppm_shock_eta2), 0.0, 1.0);

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
