#include "ghl_reconstruction.h"

/*
 * Function     : ghl_wenoz_reconstruction()
 * Description  : reconstructs variables at the points
 *                    Ur(i) = U(i-1/2+epsilon)
 *                    Ul(i) = U(i-1/2-epsilon)
 *                using the WENO-z reconstruction algorithm
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/wenoz
 *
*/

void ghl_wenoz_reconstruction(
      const double U[6],
      double *restrict Ur,
      double *restrict Ul) {

  double tmpr, tmpl;

  ghl_wenoz_reconstruction_right_left_faces(&U[0], &tmpr, &tmpl);
  *Ul = tmpr;

  ghl_wenoz_reconstruction_right_left_faces(&U[1], &tmpr, &tmpl);
  *Ur = tmpl;
 }