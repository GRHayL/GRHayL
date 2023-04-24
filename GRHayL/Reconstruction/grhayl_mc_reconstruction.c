#include "reconstruction.h"

/* Function    : grhayl_mc_reconstruction()
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

void grhayl_mc_reconstruction(const double U_m2,
                              const double U_m1,
                              const double U,
                              const double U_p1,
                              double *restrict Ur,
                              double *restrict Ul) {

  double tmp_sigma_i, tmp_sigma_im1, sigma_i, sigma_im1;

  tmp_sigma_i   = grhayl_minmod(U    - U_m1, U_p1 - U   );
  tmp_sigma_im1 = grhayl_minmod(U_m1 - U_m2, U    - U_m1);

  sigma_i   = grhayl_minmod(0.5*(U_p1 - U_m1), 2*tmp_sigma_i);
  sigma_im1 = grhayl_minmod(0.5*(U    - U_m2), 2*tmp_sigma_im1);

  *Ul = U    - 0.5*sigma_i;
  *Ur = U_m1 + 0.5*sigma_im1;
 }