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

  double sigma_i_1, sigma_i_2, sigma_im1_1, sigma_im1_2;

  sigma_i_1   = ghl_minmod( (U    - U_m1), 2.0*(U_p1 - U   ) );
  sigma_im1_1 = ghl_minmod( (U_m1 - U_m2), 2.0*(U    - U_m1) );

  sigma_i_2   = ghl_minmod( 2.0*(U    - U_m1), (U_p1 - U   ) );
  sigma_im1_2 = ghl_minmod( 2.0*(U_m1 - U_m2), (U    - U_m1) );

  double sigma_i, sigma_im1;

  sigma_i   = ghl_maxmod(sigma_i_1,   sigma_i_2);
  sigma_im1 = ghl_maxmod(sigma_im1_1, sigma_im1_2);

  *Ur = U    - 0.5*sigma_i;
  *Ul = U_m1 + 0.5*sigma_im1;
 }
