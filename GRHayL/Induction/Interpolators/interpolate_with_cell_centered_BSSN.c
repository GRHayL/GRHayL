#include "ghl_induction.h"
#include "ghl_induction_helpers.h"

/*
 * Function      : ghl_interpolate_with_cell_centered_BSSN()
 * Description   : computes several interpolated quantities for computing the RHS
 *                 for tilde{phi} and the gauge contributions to A_i; these are
 *                 used in ghl_calculate_phitilde_rhs() function
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_interpolate_with_cell_centered_BSSN
*/

void ghl_interpolate_with_cell_centered_BSSN(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double psi_stencil[2][2][2],
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
      const double phitilde,
      induction_interp_vars *restrict interp_vars) {
  /* Compute \partial_t psi6phi = -\partial_i (  \alpha psi^6 A^i - psi6phi \beta^i)
   *    (Eq 13 of http://arxiv.org/pdf/1110.4633.pdf), using Lorenz gauge.
   * Note that the RHS consists of a shift advection term on psi6phi and
   *    a term depending on the vector potential.
   * psi6phi is defined at (i+1/2,j+1/2,k+1/2), but instead of reconstructing
   *    to compute the RHS of \partial_t psi6phi, we instead use standard
   *    interpolations.
   */

  // First compute \partial_j \alpha \sqrt{\gamma} A^j (RHS of \partial_i psi6phi)
  // FIXME: Would be much cheaper & easier to unstagger A_i, raise, then interpolate A^i.
  //        However, we keep it this way to be completely compatible with the original
  //        Illinois GRMHD thorn, called mhd_evolve.

  /*
     We need to interpolate several quantities to several different points depending on the quantities
     we're computing. The staggered gridpoints for these variables are
       phitilde: (i+1/2, j+1/2, k+1/2)
       A_x:      (i,     j+1/2, k+1/2)
       A_y:      (i+1/2, j,     k+1/2)
       A_z:      (i+1/2, j+1/2, k    )
     For metric quantities, we use ghl_BSSN_cell_interp(), which computes most of the needed quantities.
     It interpolates (via averaging) the lapse and shift to phitilde's location. The metric
     is interpolated to 3 different points:
       gammaUU[0][i] is at A_x's location
       gammaUU[1][i] is at A_y's location
       gammaUU[2][i] is at A_z's location

     Similarly, the function ghl_A_i_avg() interpolates A_i to these points, storing the interpolated
     data in the 4 arrays A_to_phitilde, A_to_Ax, A_to_Ay, and A_to_Az.

     These two averaging loops are split because the stencils are of different sizes. The stencils are
     centered around the staggered point. This means that the metric have an even stencil, and the A_i
     have an odd stencil.
  */
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
