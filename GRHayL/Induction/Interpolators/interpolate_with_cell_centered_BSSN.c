#include "ghl_induction.h"
#include "ghl_induction_helpers.h"

/**
 * @ingroup mag_gauge
 * @brief Interpolate induction gauge quantities using a
 *        cell-centered BSSN metric.
 *
 * @details
 * This function computes the elements of ghl_induction_interp_vars
 * using a cell-centered BSSN metric input. The `metric_stencil` elements
 * requires some auxiliary quantities to be filled, so it is recommended
 * to use the @ref ghl_initialize_metric function to fill each element
 * of this struct. Interpolations are handled by the internal
 * functions @ref ghl_BSSN_cell_interp and @ref ghl_A_i_avg .
 *
 * These two averaging loops are split because the stencils are of
 * different sizes. The stencils are centered around the staggered
 * point. This means that the metric quantities have an even stencil,
 * and the \f$ A_i \f$ have an odd stencil.
 *
 * @param[in] metric_stencil: 3D stencil array of ghl_metric_quantities from
 *                            \f$ (i, j, k) \f$ to \f$ (i+1, j+1, k+1) \f$
 *
 * @param[in] psi_stencil:    3D stencil array of \f$ \psi^6 \f$ from
 *                            \f$ (i, j, k) \f$ to \f$ (i+1, j+1, k+1) \f$
 *
 * @param[in] Ax_stencil:     3D stencil array containing \f$ A_x \f$ from
 *                            \f$ (i-1, j-\frac{1}{2}, k-\frac{1}{2}) \f$ to
 *                            \f$ (i+1, j+\frac{3}{2}, k+\frac{3}{2}) \f$
 *
 * @param[in] Ay_stencil:     3D stencil array containing \f$ A_y \f$ from
 *                            \f$ (i-\frac{1}{2}, j-1, k-\frac{1}{2}) \f$ to
 *                            \f$ (i+\frac{3}{2}, j+1, k+\frac{3}{2}) \f$
 *
 * @param[in] Az_stencil:     3D stencil array containing \f$ A_z \f$ from
 *                            \f$ (i-\frac{1}{2}, j-\frac{1}{2}, k-1) \f$ to
 *                            \f$ (i+\frac{3}{2}, j+\frac{3}{2}, k+1) \f$
 *
 * @param[in] phitilde:       value of \f$ \tilde{\Phi} \f$ at the staggered point
 *
 * @param[out] interp_vars:   ghl_induction_interp_vars with interpolated values
 *                            needed by @ref ghl_calculate_phitilde_rhs
 *
 * @returns void
 */
void ghl_interpolate_with_cell_centered_BSSN(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double psi_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      ghl_induction_interp_vars *restrict interp_vars) {

  ghl_metric_quantities metric_interp;
  double lapse_psi2_interp[3], lapse_over_psi6_interp;
  ghl_BSSN_cell_interp(metric_stencil, psi_stencil, &metric_interp, lapse_psi2_interp, &lapse_over_psi6_interp);

  double A_to_phitilde[3], A_to_Ax[3], A_to_Ay[3], A_to_Az[3];
  ghl_A_i_avg(Ax_stencil, Ay_stencil, Az_stencil, A_to_phitilde, A_to_Ax, A_to_Ay, A_to_Az);

  // Interpolating lapse and shift to (i+1/2, j+1/2, k+1/2) is needed for computing phitilde_rhs.
  interp_vars->alpha    = metric_interp.lapse;
  interp_vars->betai[0] = metric_interp.betaU[0];
  interp_vars->betai[1] = metric_interp.betaU[1];
  interp_vars->betai[2] = metric_interp.betaU[2];

  // Compute \partial_t phitilde = -\partial_i (  \alpha psi^6 A^i - phitilde \beta^i)
  // A^x term (interpolated to (i, j+1/2, k+1/2) )
  // \alpha \sqrt{\gamma} A^x = \alpha psi^6 A^x (RHS of \partial_i psi6phi)
  // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
  interp_vars->sqrtg_Ai[0] = lapse_psi2_interp[0]*
                                     ( metric_interp.gammaUU[0][0]*A_to_Ax[0]
                                     + metric_interp.gammaUU[0][1]*A_to_Ax[1]
                                     + metric_interp.gammaUU[0][2]*A_to_Ax[2] );

  // A^y term (interpolated to (i+1/2, j, k+1/2) )
  // \alpha \sqrt{\gamma} A^y = \alpha psi^6 A^y (RHS of \partial_i psi6phi)
  // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
  interp_vars->sqrtg_Ai[1] = lapse_psi2_interp[1]*
                                     ( metric_interp.gammaUU[1][0]*A_to_Ay[0]
                                     + metric_interp.gammaUU[1][1]*A_to_Ay[1]
                                     + metric_interp.gammaUU[1][2]*A_to_Ay[2] );

  // A^z term (interpolated to (i+1/2, j+1/2, k) )
  // \alpha \sqrt{\gamma} A^z = \alpha psi^6 A^z (RHS of \partial_i psi6phi)
  // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
  interp_vars->sqrtg_Ai[2] = lapse_psi2_interp[2]*
                                     ( metric_interp.gammaUU[2][0]*A_to_Az[0]
                                     + metric_interp.gammaUU[2][1]*A_to_Az[1]
                                     + metric_interp.gammaUU[2][2]*A_to_Az[2] );

  // Next set \alpha \Phi - \beta^j A_j at (i+1/2,j+1/2,k+1/2):
  interp_vars->alpha_Phi_minus_betaj_A_j = phitilde*lapse_over_psi6_interp
                                            - ( interp_vars->betai[0]*A_to_phitilde[0]
                                              + interp_vars->betai[1]*A_to_phitilde[1]
                                              + interp_vars->betai[2]*A_to_phitilde[2] );
}
