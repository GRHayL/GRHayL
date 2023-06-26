#include "reconstruction.h"

void ghl_compute_UrUl_onevar(
      const double U[5],
      double *restrict Ur,
      double *restrict Ul) {

  const double slope_limited_dU_m1 = ghl_slope_limit(U[MINUS1] - U[MINUS2], U[PLUS_0] - U[MINUS1]);
  const double slope_limited_dU_p0 = ghl_slope_limit(U[PLUS_0] - U[MINUS1], U[PLUS_1] - U[PLUS_0]);
  const double slope_limited_dU_p1 = ghl_slope_limit(U[PLUS_1] - U[PLUS_0], U[PLUS_2] - U[PLUS_1]);

  *Ur = 0.5*(U[PLUS_1] + U[PLUS_0]) + (1.0/6.0)*(slope_limited_dU_p0 - slope_limited_dU_p1);
  *Ul = 0.5*(U[PLUS_0] + U[MINUS1]) + (1.0/6.0)*(slope_limited_dU_m1 - slope_limited_dU_p0);
}
