#include "reconstruction.h"

/* Function    : ghl_superbee_reconstruction()
 * Description : reconstructs variables at the points
 *                   Ur(i) = U(i-1/2+epsilon)
 *                   Ul(i) = U(i-1/2-epsilon)
 *               input stencils require primitives in stencil (i-2,i-1,i,i+1)
 *               along reconstruction direction
 *
 * Inputs      : U_m2 - variable at index i-2
 *             : U_m1 - variable at index i-1
 *             : U    - variable at index i
 *             : U_p1 - variable at index i+1
 *
 * Outputs     : Ur - reconstructed U on left face
 *             : Ul - reconstructed U on left face
 */

void ghl_superbee_reconstruction(
      const double U_m2,
      const double U_m1,
      const double U,
      const double U_p1,
      double *restrict Ur,
      double *restrict Ul) {

  const double Um1m2 = U_m1 - U_m2;
  const double U0m1  = U    - U_m1;
  const double Up10  = U_p1 - U;

  const double sigma_i   = ghl_maxmod(ghl_minmod(U0m1, 2.0*Up10),
                                      ghl_minmod(2.0*U0m1, Up10));

  const double sigma_im1 = ghl_maxmod(ghl_minmod(Um1m2, 2.0*U0m1),
                                      ghl_minmod(2.0*Um1m2, U0m1));

  *Ur = U    - 0.5*sigma_i;
  *Ul = U_m1 + 0.5*sigma_im1;
 }
