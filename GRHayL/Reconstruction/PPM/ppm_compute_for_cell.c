#include "ghl_reconstruction.h"

/*
 * Function     : ghl_ppm_compute_for_cell()
 * Description  : reconstructs variables at the points
 *                    Ur(i) = U(i-1/2+epsilon)
 *                    Ul(i) = U(i+1/2-epsilon)
 * Documentation:  https://github.com/GRHayL/GRHayL/wiki/ghl_ppm_compute_for_cell
*/

#define SLOPE_LIMITER_COEFF 2.0
void ghl_ppm_compute_for_cell(
      const double ftilde,
      const double U[5],
      double *restrict Ur_ptr,
      double *restrict Ul_ptr) {

  const double U0 = U[PLUS_0];

  const double slope_limited_dU_m1 = ghl_slope_limit(U[MINUS1] - U[MINUS2], U0        - U[MINUS1]);
  const double slope_limited_dU_p0 = ghl_slope_limit(U0        - U[MINUS1], U[PLUS_1] - U0);
  const double slope_limited_dU_p1 = ghl_slope_limit(U[PLUS_1] - U0,        U[PLUS_2] - U[PLUS_1]);

  double Ur = 0.5*(U[PLUS_1] + U0) + (1.0/6.0)*(slope_limited_dU_p0 - slope_limited_dU_p1);
  double Ul = 0.5*(U0 + U[MINUS1]) + (1.0/6.0)*(slope_limited_dU_m1 - slope_limited_dU_p0);

  // First detect shocks / steep gradients
  // and flatten variables
  Ur = U0*ftilde + Ur*(1.0 - ftilde);
  Ul = U0*ftilde + Ul*(1.0 - ftilde);

  // Then monotonize all variables
  if ( (Ur - U0)*(U0 - Ul) <= 0.0) {
    *Ur_ptr = U0;
    *Ul_ptr = U0;
    return;
  }

  const double dU = Ur - Ul;
  const double Utmp = dU*( U0 - 0.5*(Ur + Ul) );

  if ( Utmp > (1.0/6.0)*(dU*dU)) {
    *Ur_ptr = Ur;
    *Ul_ptr = 3.0*U0 - 2.0*Ur;
  } else if ( Utmp < -(1.0/6.0)*(dU*dU)) {
    *Ur_ptr = 3.0*U0 - 2.0*Ul;
    *Ul_ptr = Ul;
  } else {
    *Ur_ptr = Ur;
    *Ul_ptr = Ul;
  }
}
