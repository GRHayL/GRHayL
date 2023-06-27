#include "induction.h"

// Note: the following functions, while having a size argument, are specialized for array size of 3 and 2, respectively.
//       They will need to be extended to support arbitrary stencil averaging. This involves determining the element
//       of the array which corresponds to the current location. ghl_A_i_avg() assumes 1 (centered 3-element stencil),
//       and ghl_metric_avg() assumes 0 (2-element stencil centered around staggered gridpoint)

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
  /*
     The stenciling sum here is reasonably simple. The index corresponding to the component needs to be
     interpolated to +1/2 (the kk loop). Thus, the kk index is always in the same index as the component
     of A being interpolated. For the jj loop, the index corresponding to the interpolated A
     (i.e. x index for A_to_Ax) needs to be interpolated back by -1/2. Hence, the kk loop is in the range
     [1,size) and the jj loop is in the range [0,size-1).
  */
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


void ghl_ADM_cell_interp(
      const int size,
      const metric_quantities metric_stencil[size][size][size],
      metric_quantities *restrict metric_interp,
      double detg_interp[3],
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
  detg_interp[0]               = 0.0;
  detg_interp[1]               = 0.0;
  detg_interp[2]               = 0.0;
  *lapse_over_psi6_interp      = 0.0;

  for(int kk=0; kk<size; kk++) {
    for(int jj=0; jj<size; jj++) {
      // Interpolate xx, xy, xz to (i, j+1/2, k+1/2) for A_x
      metric_interp->gammaUU[0][0] += metric_stencil[kk][jj][0].gammaUU[0][0];
      metric_interp->gammaUU[0][1] += metric_stencil[kk][jj][0].gammaUU[0][1];
      metric_interp->gammaUU[0][2] += metric_stencil[kk][jj][0].gammaUU[0][2];
      detg_interp[0]               += metric_stencil[kk][jj][0].lapse*metric_stencil[kk][jj][0].sqrt_detgamma;

      // Interpolate yx, yy, yz to (i+1/2, j, k+1/2) for A_y
      metric_interp->gammaUU[1][0] += metric_stencil[kk][0][jj].gammaUU[0][1];
      metric_interp->gammaUU[1][1] += metric_stencil[kk][0][jj].gammaUU[1][1];
      metric_interp->gammaUU[1][2] += metric_stencil[kk][0][jj].gammaUU[1][2];
      detg_interp[1]               += metric_stencil[kk][0][jj].lapse*metric_stencil[kk][0][jj].sqrt_detgamma;

      // Interpolate zx, zy, zz to (i+1/2, j+1/2, k) for A_z
      metric_interp->gammaUU[2][0] += metric_stencil[0][kk][jj].gammaUU[0][2];
      metric_interp->gammaUU[2][1] += metric_stencil[0][kk][jj].gammaUU[1][2];
      metric_interp->gammaUU[2][2] += metric_stencil[0][kk][jj].gammaUU[2][2];
      detg_interp[2]               += metric_stencil[0][kk][jj].lapse*metric_stencil[0][kk][jj].sqrt_detgamma;

      for(int ii=0; ii<size; ii++) {
        metric_interp->lapse    += metric_stencil[kk][jj][ii].lapse;
        metric_interp->betaU[0] += metric_stencil[kk][jj][ii].betaU[0];
        metric_interp->betaU[1] += metric_stencil[kk][jj][ii].betaU[1];
        metric_interp->betaU[2] += metric_stencil[kk][jj][ii].betaU[2];
        *lapse_over_psi6_interp += pow(metric_stencil[kk][jj][ii].lapse,2)/metric_stencil[kk][jj][ii].sqrt_detgamma;
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
  detg_interp[0]               /= sum_2D;
  detg_interp[1]               /= sum_2D;
  detg_interp[2]               /= sum_2D;
}
