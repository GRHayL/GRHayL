#include "induction.h"

// Helper function for interpolating quantities via averaging
void grhayl_A_i_avg(
      const int size,
      const double Ax_stencil[size][size][size],
      const double Ay_stencil[size][size][size],
      const double Az_stencil[size][size][size],
      double A_to_phitilde[3],
      double A_to_Ax[3],
      double A_to_Ay[3],
      double A_to_Az[3]);

void grhayl_metric_avg(
      const int size,
      const metric_quantities metric_stencil[size][size][size],
      const double psi_stencil[size][size][size],
      metric_quantities *restrict metric_interp,
      double lapse_psi2_interp[3],
      double *restrict lapse_over_psi6_interp);

/* Function    : grhayl_interpolate_for_A_gauge_rhs()
 * Description : computes several interpolated quantities for computing the RHS
 *               for tilde{phi} and the gauge contributions to A_i; these are
 *               used in grhayl_calculate_phitilde_rhs() function
 *
 * Inputs      : gauge_vars      - A_gauge_vars struct containing stencils for
 *                                 variables to compute interpolated values
 *
 * Outputs     : interp_vars  - induction_interp_vars struct containing interpolated
 *                                 values for use in
 *
 */
void grhayl_interpolate_for_induction_rhs(
      const metric_quantities metric_stencil[2][2][2],
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
  //
  //Step 1) j=x: Need to raise A_i, but to do that, we must have all variables at the same gridpoints:
  // The goal is to compute \partial_j (\alpha \sqrt{\gamma} A^j) at (i+1/2,j+1/2,k+1/2)
  //    We do this by first interpolating (RHS1x) = (\alpha \sqrt{\gamma} A^x) at
  //    (i,j+1/2,k+1/2)and (i+1,j+1/2,k+1/2), then taking \partial_x (RHS1x) =
  //    [ RHS1x(i+1,j+1/2,k+1/2) - RHS1x(i,j+1/2,k+1/2) ]/dX.


  /*
     We need to interpolate several quantities to several different points depending on the quantities
     we're computing. The staggered gridpoints for these variables are
       phitilde: (i+1/2, j+1/2, k+1/2)
       A_x:      (i,     j+1/2, k+1/2)
       A_y:      (i+1/2, j,     k+1/2)
       A_z:      (i+1/2, j+1/2, k    )
     For metric quantities, we use grhayl_metric_avg(), which computes most of the needed quantities.
     It interpolates (via averaging) the lapse and shift to phitilde's location. The metric
     is interpolated to 3 different points:
       gammaUU[0][i] is a A_x's location
       gammaUU[1][i] is a A_y's location
       gammaUU[2][i] is a A_z's location

     Similarly, the function grhayl_A_i_avg() interpolates A_i to these points, storing the interpolated
     data in the 4 arrays A_to_phitilde, A_to_Ax, A_to_Ay, and A_to_Az.

     These two averaging loops are split because the stencils are of different sizes. The stencils are
     centered around the staggered point. This means that the metric have an even stencil, and the A_i
     have an odd stencil.
  */
  metric_quantities metric_interp;
  double lapse_psi2_interp[3], lapse_over_psi6_interp;
  grhayl_metric_avg(2, metric_stencil, psi_stencil, &metric_interp, lapse_psi2_interp, &lapse_over_psi6_interp);

  double A_to_phitilde[3], A_to_Ax[3], A_to_Ay[3], A_to_Az[3];
  grhayl_A_i_avg(3, Ax_stencil, Ay_stencil, Az_stencil, A_to_phitilde, A_to_Ax, A_to_Ay, A_to_Az);

  // Interpolating lapse and shift to (i+1/2, j+1/2, k+1/2) is needed for computing phitilde_rhs.
  interp_vars->alpha_interp  = metric_interp.lapse;
  interp_vars->shifti_interp[0] = metric_interp.betaU[0];
  interp_vars->shifti_interp[1] = metric_interp.betaU[1];
  interp_vars->shifti_interp[2] = metric_interp.betaU[2];

  // A^x term (interpolated to (i, j+1/2, k+1/2) )
  // \alpha \sqrt{\gamma} A^x = \alpha psi^6 A^x (RHS of \partial_i psi6phi)
  // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
  interp_vars->alpha_sqrtg_Ai_interp[0] = lapse_psi2_interp[0]*
                                            ( metric_interp.gammaUU[0][0]*A_to_Ax[0]
                                            + metric_interp.gammaUU[0][1]*A_to_Ax[1]
                                            + metric_interp.gammaUU[0][2]*A_to_Ax[2] );

  // A^y term (interpolated to (i+1/2, j, k+1/2) )
  // \alpha \sqrt{\gamma} A^y = \alpha psi^6 A^y (RHS of \partial_i psi6phi)
  // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
  interp_vars->alpha_sqrtg_Ai_interp[1] = lapse_psi2_interp[1]*
                                            ( metric_interp.gammaUU[1][0]*A_to_Ay[0]
                                            + metric_interp.gammaUU[1][1]*A_to_Ay[1]
                                            + metric_interp.gammaUU[1][2]*A_to_Ay[2] );

  // A^z term (interpolated to (i+1/2, j+1/2, k) )
  // \alpha \sqrt{\gamma} A^z = \alpha psi^6 A^z (RHS of \partial_i psi6phi)
  // Note that gupij is \tilde{\gamma}^{ij}, so we need to multiply by \psi^{-4}.
  interp_vars->alpha_sqrtg_Ai_interp[2] = lapse_psi2_interp[2]*
                                            ( metric_interp.gammaUU[2][0]*A_to_Az[0]
                                            + metric_interp.gammaUU[2][1]*A_to_Az[1]
                                            + metric_interp.gammaUU[2][2]*A_to_Az[2] );

  // Next set \alpha \Phi - \beta^j A_j at (i+1/2,j+1/2,k+1/2):
  interp_vars->alpha_Phi_minus_betaj_A_j_interp = phitilde*lapse_over_psi6_interp
                                                   - ( interp_vars->shifti_interp[0]*A_to_phitilde[0]
                                                     + interp_vars->shifti_interp[1]*A_to_phitilde[1]
                                                     + interp_vars->shifti_interp[2]*A_to_phitilde[2] );
}

// Note: the following functions, while having a size argument, are specialized for array size of 3 and 2, respectively.
//       They will need to be extended to support arbitrary stencil averaging. This involves determining the element
//       of the array which corresponds to the current location. grhayl_A_i_avg() assumes 1 (centered 3-element stencil),
//       and grhayl_metric_avg() assumes 0 (2-element stencil centered around staggered gridpoint)
void grhayl_A_i_avg(
      const int size,
      const double Ax_stencil[size][size][size],
      const double Ay_stencil[size][size][size],
      const double Az_stencil[size][size][size],
      double A_to_phitilde[3],
      double A_to_Ax[3],
      double A_to_Ay[3],
      double A_to_Az[3]) {
  
  A_to_phitilde[0] = 0.0;
  A_to_phitilde[1] = 0.0;
  A_to_phitilde[2] = 0.0;
  A_to_Ax[1]       = 0.0;
  A_to_Ax[2]       = 0.0;
  A_to_Ay[0]       = 0.0;
  A_to_Ay[2]       = 0.0;
  A_to_Az[0]       = 0.0;
  A_to_Az[1]       = 0.0;

  // A_i is already at the location of A_i :)
  A_to_Ax[0] = Ax_stencil[1][1][1];
  A_to_Ay[1] = Ay_stencil[1][1][1];
  A_to_Az[2] = Az_stencil[1][1][1];

  int sum_1D = 0;
  int sum_2D = 0;
  // The stenciling sum here is reasonably simple. The index corresponding to the component needs to be
  // interpolated to +1/2 (the kk loop). Thus, the kk index is always in the same index as the component
  // of A being interpolated. For the jj loop, the index corresponding to the interpolated A
  // (i.e. x index for A_to_Ax) needs to be interpolated back by -1/2. Hence, the kk loop is in the range
  // [1,size) and the jj loop is in the range [0,size-1).
  for(int kk=1; kk<size; kk++) {
    A_to_phitilde[0] += Ax_stencil[1][1][kk];
    A_to_phitilde[1] += Ay_stencil[1][kk][1];
    A_to_phitilde[2] += Az_stencil[kk][1][1];
    sum_1D++;
    for(int jj=0; jj<size-1; jj++) {
      A_to_Ax[1] += Ay_stencil[1][kk][jj];
      A_to_Ax[2] += Az_stencil[kk][1][jj];

      A_to_Ay[0] += Ax_stencil[1][jj][kk];
      A_to_Ay[2] += Az_stencil[kk][jj][1];

      A_to_Az[0] += Ax_stencil[jj][1][kk];
      A_to_Az[1] += Ay_stencil[jj][kk][1];
      sum_2D++;
    }
  }
  A_to_phitilde[0] /= sum_1D;
  A_to_phitilde[1] /= sum_1D;
  A_to_phitilde[2] /= sum_1D;
  A_to_Ax[1]       /= sum_2D;
  A_to_Ax[2]       /= sum_2D;
  A_to_Ay[0]       /= sum_2D;
  A_to_Ay[2]       /= sum_2D;
  A_to_Az[0]       /= sum_2D;
  A_to_Az[1]       /= sum_2D;
}

void grhayl_metric_avg(
      const int size,
      const metric_quantities metric_stencil[size][size][size],
      const double psi_stencil[size][size][size],
      metric_quantities *restrict metric_interp,
      double lapse_psi2_interp[3],
      double *restrict lapse_over_psi6_interp) {

  int sum_2D = 0;
  int sum_3D = 0;
  metric_interp->lapse         = 0.0;
  metric_interp->betaU[0]      = 0.0;
  metric_interp->betaU[1]      = 0.0;
  metric_interp->betaU[2]      = 0.0;
  metric_interp->gammaUU[0][0] = 0.0;
  metric_interp->gammaUU[0][1] = 0.0;
  metric_interp->gammaUU[0][2] = 0.0;
  metric_interp->gammaUU[1][0] = 0.0;
  metric_interp->gammaUU[1][1] = 0.0;
  metric_interp->gammaUU[1][2] = 0.0;
  metric_interp->gammaUU[2][0] = 0.0;
  metric_interp->gammaUU[2][1] = 0.0;
  metric_interp->gammaUU[2][2] = 0.0;
  lapse_psi2_interp[0]         = 0.0;
  lapse_psi2_interp[1]         = 0.0;
  lapse_psi2_interp[2]         = 0.0;
  *lapse_over_psi6_interp      = 0.0;

  double lapse_psi2[size][size][size];
  double lapse_over_psi6[size][size][size];
  for(int kk=0; kk<size; kk++) {
    for(int jj=0; jj<size; jj++) {
      for(int ii=0; ii<size; ii++) {
        const double psi2 = psi_stencil[kk][jj][ii]*psi_stencil[kk][jj][ii];
        lapse_psi2[kk][jj][ii] = metric_stencil[kk][jj][ii].lapse*psi2;
        lapse_over_psi6[kk][jj][ii] = metric_stencil[kk][jj][ii].lapse/(psi2*psi2*psi2);
      }
    }
  }

  for(int kk=0; kk<size; kk++) {
    for(int jj=0; jj<size; jj++) {

      // Interpolate xx, xy, xz to (i, j+1/2, k+1/2) for A_x
      metric_interp->gammaUU[0][0] += metric_stencil[kk][jj][0].gammaUU[0][0];
      metric_interp->gammaUU[0][1] += metric_stencil[kk][jj][0].gammaUU[0][1];
      metric_interp->gammaUU[0][2] += metric_stencil[kk][jj][0].gammaUU[0][2];
      lapse_psi2_interp[0]         += lapse_psi2[kk][jj][0];

      // Interpolate yx, yy, yz to (i+1/2, j, k+1/2) for A_y
      metric_interp->gammaUU[1][0] += metric_stencil[kk][0][jj].gammaUU[0][1];
      metric_interp->gammaUU[1][1] += metric_stencil[kk][0][jj].gammaUU[1][1];
      metric_interp->gammaUU[1][2] += metric_stencil[kk][0][jj].gammaUU[1][2];
      lapse_psi2_interp[1]         += lapse_psi2[kk][0][jj];

      // Interpolate zx, zy, zz to (i+1/2, j+1/2, k) for A_z
      metric_interp->gammaUU[2][0] += metric_stencil[0][kk][jj].gammaUU[0][2];
      metric_interp->gammaUU[2][1] += metric_stencil[0][kk][jj].gammaUU[1][2];
      metric_interp->gammaUU[2][2] += metric_stencil[0][kk][jj].gammaUU[2][2];
      lapse_psi2_interp[2]         += lapse_psi2[0][kk][jj];

      for(int ii=0; ii<size; ii++) {
        metric_interp->lapse    += metric_stencil[kk][jj][ii].lapse;
        metric_interp->betaU[0] += metric_stencil[kk][jj][ii].betaU[0];
        metric_interp->betaU[1] += metric_stencil[kk][jj][ii].betaU[1];
        metric_interp->betaU[2] += metric_stencil[kk][jj][ii].betaU[2];
        *lapse_over_psi6_interp += lapse_over_psi6[kk][jj][ii];
        sum_3D++;
      }
      sum_2D++;
    }
  }
  metric_interp->lapse    /= sum_3D;
  metric_interp->betaU[0] /= sum_3D;
  metric_interp->betaU[1] /= sum_3D;
  metric_interp->betaU[2] /= sum_3D;
  *lapse_over_psi6_interp /= sum_3D;

  metric_interp->gammaUU[0][0] /= sum_2D;
  metric_interp->gammaUU[0][1] /= sum_2D;
  metric_interp->gammaUU[0][2] /= sum_2D;
  metric_interp->gammaUU[1][0] /= sum_2D;
  metric_interp->gammaUU[1][1] /= sum_2D;
  metric_interp->gammaUU[1][2] /= sum_2D;
  metric_interp->gammaUU[2][0] /= sum_2D;
  metric_interp->gammaUU[2][1] /= sum_2D;
  metric_interp->gammaUU[2][2] /= sum_2D;
  lapse_psi2_interp[0]         /= sum_2D;
  lapse_psi2_interp[1]         /= sum_2D;
  lapse_psi2_interp[2]         /= sum_2D;
}
