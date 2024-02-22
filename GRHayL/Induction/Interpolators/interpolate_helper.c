#include "induction.h"
#include "ghl_induction_helpers.h"

/*
   Recall that
       Ax is staggered to (i, j+1/2, k+1/2)
       Ay is staggered to (i+1/2, j, k+1/2)
       Az is staggered to (i+1/2, j, k+1/2)
       \tilde{\Phi} is staggered to (i+1/2, j+1/2, k+1/2)
   We need the A_i values at all four of these gridpoints, so we interpolate the three
   components and store them in A_to_location variables.
*/
void ghl_A_i_avg(
      const double Ax_stencil[3][3][3],
      const double Ay_stencil[3][3][3],
      const double Az_stencil[3][3][3],
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

  /*
     The stenciling sum here is reasonably simple. The index corresponding to the component needs to be
     interpolated to +1/2 (the kk loop). Thus, the kk index is always in the same index as the component
     of A being interpolated. For the jj loop, the index corresponding to the interpolated A
     (i.e. x index for A_to_Ax) needs to be interpolated back by -1/2. Hence, the kk loop is in the range
     [1,size) and the jj loop is in the range [0,size-1).
  */
  for(int kk=1; kk<3; kk++) {
    A_to_phitilde[0] += Ax_stencil[1][1][kk];
    A_to_phitilde[1] += Ay_stencil[1][kk][1];
    A_to_phitilde[2] += Az_stencil[kk][1][1];
    for(int jj=0; jj<2; jj++) {
      A_to_Ax[1] += Ay_stencil[1][kk][jj];
      A_to_Ax[2] += Az_stencil[kk][1][jj];

      A_to_Ay[0] += Ax_stencil[1][jj][kk];
      A_to_Ay[2] += Az_stencil[kk][jj][1];

      A_to_Az[0] += Ax_stencil[jj][1][kk];
      A_to_Az[1] += Ay_stencil[jj][kk][1];
    }
  }
  A_to_phitilde[0] /= 2.0;
  A_to_phitilde[1] /= 2.0;
  A_to_phitilde[2] /= 2.0;
  A_to_Ax[1]       /= 4.0;
  A_to_Ax[2]       /= 4.0;
  A_to_Ay[0]       /= 4.0;
  A_to_Ay[2]       /= 4.0;
  A_to_Az[0]       /= 4.0;
  A_to_Az[1]       /= 4.0;
}

/*
   All of the following code aims to compute metric quantities
   at several points:
     1) lapse psi^2 (BSSN) or sqrt{-g} (ADM) at each A_i gridpoint
     2) lapse, shift, and either lapse / psi^6 (BSSN) or lapse^2 / sqrt{gamma}
        at the location of \tilde{\Phi}
   The cell-centered metric code requires more interpolations since \tilde{\Phi}
   is vertex-centered.
*/
void ghl_BSSN_cell_interp(
      const ghl_metric_quantities metric_stencil[2][2][2],
      const double psi_stencil[2][2][2],
      ghl_metric_quantities *restrict metric_interp,
      double lapse_psi2_interp[3],
      double *restrict lapse_over_psi6_interp) {

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

  double lapse_psi2[2][2][2];
  double lapse_over_psi6[2][2][2];
  for(int kk=0; kk<2; kk++) {
    for(int jj=0; jj<2; jj++) {
      for(int ii=0; ii<2; ii++) {
        const double psi2 = psi_stencil[kk][jj][ii]*psi_stencil[kk][jj][ii];
        lapse_psi2[kk][jj][ii] = metric_stencil[kk][jj][ii].lapse*psi2;
        lapse_over_psi6[kk][jj][ii] = metric_stencil[kk][jj][ii].lapse/(psi2*psi2*psi2);
      }
    }
  }

  for(int kk=0; kk<2; kk++) {
    for(int jj=0; jj<2; jj++) {
      // Interpolate xx, xy, xz from ccc to cvv centering for A_x
      metric_interp->gammaUU[0][0] += metric_stencil[kk][jj][0].gammaUU[0][0];
      metric_interp->gammaUU[0][1] += metric_stencil[kk][jj][0].gammaUU[0][1];
      metric_interp->gammaUU[0][2] += metric_stencil[kk][jj][0].gammaUU[0][2];
      lapse_psi2_interp[0]         += lapse_psi2[kk][jj][0];

      // Interpolate yx, yy, yz from ccc to vcv centering for A_y
      metric_interp->gammaUU[1][0] += metric_stencil[kk][0][jj].gammaUU[0][1];
      metric_interp->gammaUU[1][1] += metric_stencil[kk][0][jj].gammaUU[1][1];
      metric_interp->gammaUU[1][2] += metric_stencil[kk][0][jj].gammaUU[1][2];
      lapse_psi2_interp[1]         += lapse_psi2[kk][0][jj];

      // Interpolate zx, zy, zz from ccc to vvc centering for A_z
      metric_interp->gammaUU[2][0] += metric_stencil[0][kk][jj].gammaUU[0][2];
      metric_interp->gammaUU[2][1] += metric_stencil[0][kk][jj].gammaUU[1][2];
      metric_interp->gammaUU[2][2] += metric_stencil[0][kk][jj].gammaUU[2][2];
      lapse_psi2_interp[2]         += lapse_psi2[0][kk][jj];

      for(int ii=0; ii<2; ii++) {
        metric_interp->lapse    += metric_stencil[kk][jj][ii].lapse;
        metric_interp->betaU[0] += metric_stencil[kk][jj][ii].betaU[0];
        metric_interp->betaU[1] += metric_stencil[kk][jj][ii].betaU[1];
        metric_interp->betaU[2] += metric_stencil[kk][jj][ii].betaU[2];
        *lapse_over_psi6_interp += lapse_over_psi6[kk][jj][ii];
      }
    }
  }
  metric_interp->lapse    /= 8.0;
  metric_interp->betaU[0] /= 8.0;
  metric_interp->betaU[1] /= 8.0;
  metric_interp->betaU[2] /= 8.0;
  *lapse_over_psi6_interp /= 8.0;

  metric_interp->gammaUU[0][0] /= 4.0;
  metric_interp->gammaUU[0][1] /= 4.0;
  metric_interp->gammaUU[0][2] /= 4.0;
  metric_interp->gammaUU[1][0] /= 4.0;
  metric_interp->gammaUU[1][1] /= 4.0;
  metric_interp->gammaUU[1][2] /= 4.0;
  metric_interp->gammaUU[2][0] /= 4.0;
  metric_interp->gammaUU[2][1] /= 4.0;
  metric_interp->gammaUU[2][2] /= 4.0;
  lapse_psi2_interp[0]         /= 4.0;
  lapse_psi2_interp[1]         /= 4.0;
  lapse_psi2_interp[2]         /= 4.0;
}

void ghl_ADM_cell_interp(
      const ghl_metric_quantities metric_stencil[2][2][2],
      ghl_metric_quantities *restrict metric_interp,
      double *restrict lapse_over_psi6_interp) {

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
  *lapse_over_psi6_interp      = 0.0;

  for(int kk=0; kk<2; kk++) {
    for(int jj=0; jj<2; jj++) {
      // Interpolate xx, xy, xz from ccc to cvv centering for A_x
      const double detgx = metric_stencil[kk][jj][0].lapse*metric_stencil[kk][jj][0].sqrt_detgamma;
      metric_interp->gammaUU[0][0] += detgx*metric_stencil[kk][jj][0].gammaUU[0][0];
      metric_interp->gammaUU[0][1] += detgx*metric_stencil[kk][jj][0].gammaUU[0][1];
      metric_interp->gammaUU[0][2] += detgx*metric_stencil[kk][jj][0].gammaUU[0][2];

      // Interpolate yx, yy, yz from ccc to vcv centering for A_y
      const double detgy = metric_stencil[kk][0][jj].lapse*metric_stencil[kk][0][jj].sqrt_detgamma;
      metric_interp->gammaUU[1][0] += detgy*metric_stencil[kk][0][jj].gammaUU[0][1];
      metric_interp->gammaUU[1][1] += detgy*metric_stencil[kk][0][jj].gammaUU[1][1];
      metric_interp->gammaUU[1][2] += detgy*metric_stencil[kk][0][jj].gammaUU[1][2];

      // Interpolate zx, zy, zz from ccc to vvc centering for A_z
      const double detgz = metric_stencil[0][kk][jj].lapse*metric_stencil[0][kk][jj].sqrt_detgamma;
      metric_interp->gammaUU[2][0] += detgz*metric_stencil[0][kk][jj].gammaUU[0][2];
      metric_interp->gammaUU[2][1] += detgz*metric_stencil[0][kk][jj].gammaUU[1][2];
      metric_interp->gammaUU[2][2] += detgz*metric_stencil[0][kk][jj].gammaUU[2][2];

      for(int ii=0; ii<2; ii++) {
        metric_interp->lapse    += metric_stencil[kk][jj][ii].lapse;
        metric_interp->betaU[0] += metric_stencil[kk][jj][ii].betaU[0];
        metric_interp->betaU[1] += metric_stencil[kk][jj][ii].betaU[1];
        metric_interp->betaU[2] += metric_stencil[kk][jj][ii].betaU[2];
        *lapse_over_psi6_interp += metric_stencil[kk][jj][ii].lapse/metric_stencil[kk][jj][ii].sqrt_detgamma;
      }
    }
  }
  metric_interp->lapse    /= 8.0;
  metric_interp->betaU[0] /= 8.0;
  metric_interp->betaU[1] /= 8.0;
  metric_interp->betaU[2] /= 8.0;
  *lapse_over_psi6_interp /= 8.0;

  metric_interp->gammaUU[0][0] /= 4.0;
  metric_interp->gammaUU[0][1] /= 4.0;
  metric_interp->gammaUU[0][2] /= 4.0;
  metric_interp->gammaUU[1][0] /= 4.0;
  metric_interp->gammaUU[1][1] /= 4.0;
  metric_interp->gammaUU[1][2] /= 4.0;
  metric_interp->gammaUU[2][0] /= 4.0;
  metric_interp->gammaUU[2][1] /= 4.0;
  metric_interp->gammaUU[2][2] /= 4.0;
}

void ghl_ADM_vertex_interp(
      const ghl_metric_quantities metric_stencil[2][2][2],
      double gammaUU_interp[3][3]) {

  gammaUU_interp[0][0] = 0.0;
  gammaUU_interp[0][1] = 0.0;
  gammaUU_interp[0][2] = 0.0;
  gammaUU_interp[1][0] = 0.0;
  gammaUU_interp[1][1] = 0.0;
  gammaUU_interp[1][2] = 0.0;
  gammaUU_interp[2][0] = 0.0;
  gammaUU_interp[2][1] = 0.0;
  gammaUU_interp[2][2] = 0.0;

  for(int ii=0; ii<2; ii++) {
    // Interpolate xx, xy, xz from vvv to cvv centering for A_x
    const double detgx = metric_stencil[0][0][ii].lapse*metric_stencil[0][0][ii].sqrt_detgamma;
    gammaUU_interp[0][0] += detgx*metric_stencil[0][0][ii].gammaUU[0][0];
    gammaUU_interp[0][1] += detgx*metric_stencil[0][0][ii].gammaUU[0][1];
    gammaUU_interp[0][2] += detgx*metric_stencil[0][0][ii].gammaUU[0][2];

    // Interpolate yx, yy, yz from vvv to vcv centering for A_y
    const double detgy = metric_stencil[0][ii][0].lapse*metric_stencil[0][ii][0].sqrt_detgamma;
    gammaUU_interp[1][0] += detgy*metric_stencil[0][ii][0].gammaUU[0][1];
    gammaUU_interp[1][1] += detgy*metric_stencil[0][ii][0].gammaUU[1][1];
    gammaUU_interp[1][2] += detgy*metric_stencil[0][ii][0].gammaUU[1][2];

    // Interpolate zx, zy, zz from vvv to vvc centering for A_z
    const double detgz = metric_stencil[ii][0][0].lapse*metric_stencil[ii][0][0].sqrt_detgamma;
    gammaUU_interp[2][0] += detgz*metric_stencil[ii][0][0].gammaUU[0][2];
    gammaUU_interp[2][1] += detgz*metric_stencil[ii][0][0].gammaUU[1][2];
    gammaUU_interp[2][2] += detgz*metric_stencil[ii][0][0].gammaUU[2][2];
  }
  gammaUU_interp[0][0] /= 2.0;
  gammaUU_interp[0][1] /= 2.0;
  gammaUU_interp[0][2] /= 2.0;
  gammaUU_interp[1][0] /= 2.0;
  gammaUU_interp[1][1] /= 2.0;
  gammaUU_interp[1][2] /= 2.0;
  gammaUU_interp[2][0] /= 2.0;
  gammaUU_interp[2][1] /= 2.0;
  gammaUU_interp[2][2] /= 2.0;
}
