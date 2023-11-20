#include "reconstruction.h"

/*
 * Function     : ghl_superbee_reconstruction()
 * Description  : reconstructs variables at the points
 *                    Ur(i) = U(i-1/2+epsilon)
 *                    Ul(i) = U(i-1/2-epsilon)
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_superbee_reconstruction
*/

void ghl_superbee_reconstruction(
      const double U[4],
      double *restrict Ur,
      double *restrict Ul) {

  const double Um1m2 = U[1] - U[0];
  const double U0m1  = U[2] - U[1];
  const double Up10  = U[3] - U[2];

  const double sigma_i   = ghl_maxmod(ghl_minmod(U0m1, 2.0*Up10),
                                      ghl_minmod(2.0*U0m1, Up10));

  const double sigma_im1 = ghl_maxmod(ghl_minmod(Um1m2, 2.0*U0m1),
                                      ghl_minmod(2.0*Um1m2, U0m1));

  *Ur = U[2] - 0.5*sigma_i;
  *Ul = U[1] + 0.5*sigma_im1;
 }
