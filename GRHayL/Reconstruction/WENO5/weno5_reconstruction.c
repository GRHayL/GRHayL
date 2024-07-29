#include "reconstruction.h"

/*
 * Function     : ghl_weno5_reconstruction()
 * Description  : reconstructs variables at the points
 *                    Ur(i) = U(i-1/2+epsilon)
 *                    Ul(i) = U(i-1/2-epsilon)
 *                using the WENO-5 reconstruction algorithm
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/weno5
 * 
*/

void ghl_weno5_reconstruction(
      const double U[6],
      double *restrict Ul,
      double *restrict Ur) {

  double tmpr, tmpl;        

  ghl_weno5_reconstruction_right_left_faces( U,    &tmpr, &tmpl);
  *Ul = tmpr;
  
  ghl_weno5_reconstruction_right_left_faces(&U[1], &tmpr, &tmpl);
  *Ur = tmpl;
 }