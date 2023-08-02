#include "reconstruction.h"

/*
 * Function     : ghl_mc_reconstruction()
 * Description  : reconstructs variables at the points
 *                    Ur(i) = U(i-1/2+epsilon)
 *                    Ul(i) = U(i-1/2-epsilon)
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_mc_reconstruction
*/

void ghl_mc_reconstruction(
      const double U_m2,
      const double U_m1,
      const double U,
      const double U_p1,
      double *restrict Ur,
      double *restrict Ul) {

  const double DeltaU_0_m1 = U - U_m1;

  const double tmp_sigma_i   = ghl_minmod(DeltaU_0_m1, U_p1 - U   );
  const double tmp_sigma_im1 = ghl_minmod(U_m1 - U_m2, DeltaU_0_m1);

  const double sigma_i   = ghl_minmod(0.5*(U_p1 - U_m1), 2*tmp_sigma_i);
  const double sigma_im1 = ghl_minmod(0.5*(U    - U_m2), 2*tmp_sigma_im1);

  *Ur = U    - 0.5*sigma_i;
  *Ul = U_m1 + 0.5*sigma_im1;
 }
