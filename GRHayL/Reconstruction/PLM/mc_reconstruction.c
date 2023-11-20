#include "reconstruction.h"

/*
 * Function     : ghl_mc_reconstruction()
 * Description  : reconstructs variables at the points
 *                    Ur(i) = U(i-1/2+epsilon)
 *                    Ul(i) = U(i-1/2-epsilon)
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_mc_reconstruction
*/

void ghl_mc_reconstruction(
      const double U[4],
      double *restrict Ur,
      double *restrict Ul) {

  const double DeltaU_0_m1 = U[2] - U[1];

  const double tmp_sigma_i   = ghl_minmod(DeltaU_0_m1, U[3] - U[2]);
  const double tmp_sigma_im1 = ghl_minmod(U[1] - U[0], DeltaU_0_m1);

  const double sigma_i   = ghl_minmod(0.5*(U[3] - U[1]), 2*tmp_sigma_i);
  const double sigma_im1 = ghl_minmod(0.5*(U[2] - U[0]), 2*tmp_sigma_im1);

  *Ur = U[2] - 0.5*sigma_i;
  *Ul = U[1] + 0.5*sigma_im1;
 }
