#include "ghl_induction.h"
#include "ghl_induction_helpers.h"

/**
 * @ingroup mag_gauge
 * @brief Interpolate induction gauge quantities using a
 *        vertex-centered ADM metric.
 *
 * @details
 * This function computes the elements of ghl_induction_interp_vars
 * using a vertex-centered ADM metric input. The `metric_stencil` elements
 * requires some auxiliary quantities to be filled, so it is recommended
 * to use the @ref ghl_initialize_metric function to fill each element
 * of this struct. Interpolations are handled by the internal
 * functions @ref ghl_ADM_vertex_interp and @ref ghl_A_i_avg .
 *
 * These two averaging loops are split because the stencils are of
 * different sizes. The stencils are centered around the staggered
 * point. This means that the metric quantities have an even stencil,
 * and the \f$ A_i \f$ have an odd stencil.
 *
 * We also note for completion that the metric interpolation returns
 * \f$ \sqrt{\gamma}\gamma^{ij} \f$ instead of \f$ \gamma^{ij} \f$
 * to condense some of the mathematical operations needed.
 *
 * @param[in] metric_stencil: 3D stencil array of ghl_metric_quantities from
 *                            \f$ (i-\frac{1}{2}, j-\frac{1}{2}, k-\frac{1}{2}) \f$ to
 *                            \f$ (i+\frac{3}{2}, j+\frac{3}{2}, k+\frac{3}{2}) \f$
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
void ghl_interpolate_with_vertex_centered_ADM(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      ghl_induction_interp_vars *restrict interp_vars) {

  double gammaUU_interp[3][3];
  ghl_ADM_vertex_interp(metric_stencil, gammaUU_interp);

  double A_to_phitilde[3], A_to_Ax[3], A_to_Ay[3], A_to_Az[3];
  ghl_A_i_avg(Ax_stencil, Ay_stencil, Az_stencil, A_to_phitilde, A_to_Ax, A_to_Ay, A_to_Az);

  // Compute \partial_t phitilde = -\partial_i (  \alpha psi^6 A^i - phitilde \beta^i)
  // A^x term (interpolated to (i, j+1/2, k+1/2) )
  // \sqrt{-g} A^x = \alpha \sqrt{\gamma} A^x (RHS of \partial_i psi6phi)
  interp_vars->sqrtg_Ai[0] = gammaUU_interp[0][0]*A_to_Ax[0]
                           + gammaUU_interp[0][1]*A_to_Ax[1]
                           + gammaUU_interp[0][2]*A_to_Ax[2];

  // A^y term (interpolated to (i+1/2, j, k+1/2) )
  // \sqrt{-g} A^y = \alpha \sqrt{\gamma} A^y (RHS of \partial_i psi6phi)
  interp_vars->sqrtg_Ai[1] = gammaUU_interp[1][0]*A_to_Ay[0]
                           + gammaUU_interp[1][1]*A_to_Ay[1]
                           + gammaUU_interp[1][2]*A_to_Ay[2];

  // A^z term (interpolated to (i+1/2, j+1/2, k) )
  // \sqrt{-g} A^z = \alpha \sqrt{\gamma} A^z (RHS of \partial_i psi6phi)
  interp_vars->sqrtg_Ai[2] = gammaUU_interp[2][0]*A_to_Az[0]
                           + gammaUU_interp[2][1]*A_to_Az[1]
                           + gammaUU_interp[2][2]*A_to_Az[2];

  // Next set \alpha \Phi - \beta^j A_j at (i+1/2,j+1/2,k+1/2)
  // \alpha \Phi = \alpha \tilde{\Phi} / psi^6
  //             = \alpha \tilde{\Phi} / \sqrt{\gamma}
  interp_vars->alpha_Phi_minus_betaj_A_j = phitilde*metric_stencil[0][0][0].lapse/metric_stencil[0][0][0].sqrt_detgamma
                                            - ( metric_stencil[0][0][0].betaU[0]*A_to_phitilde[0]
                                              + metric_stencil[0][0][0].betaU[1]*A_to_phitilde[1]
                                              + metric_stencil[0][0][0].betaU[2]*A_to_phitilde[2] );
}
