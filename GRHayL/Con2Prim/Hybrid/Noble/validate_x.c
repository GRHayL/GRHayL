#include "../../utils_Noble.h"

void ghl_validate_1D(
      const harm_aux_vars_struct *restrict harm_aux,
      const double dx[1],
      const double x_old[1],
      double x[1]) {

  // Try the full step, then fixed fractions back to the old iterate. Requiring
  // positive Z and v^2 < 1 keeps every accepted trial physical. Never cross
  // x_old; if even the old iterate is invalid, report failure to the driver.
  for(int backstep = 0; backstep <= 10; backstep++) {
    if(backstep > 0) {
      x[0] = x_old[0] + (1.0 - 0.1 * backstep) * dx[0];
    }
    if(isfinite(x[0]) && x[0] > 0.0) {
      const double vsq = ghl_vsq_calc(harm_aux, x[0]);
      if(isfinite(vsq) && vsq < 1.0) {
        return;
      }
    }
  }
  x[0] = NAN;
}

void ghl_validate_1D_entropy(
      const harm_aux_vars_struct *restrict harm_aux,
      const double dx[1],
      const double x_old[1],
      double x[1]) {

  x[0] = fabs(x[0]);
}

void ghl_validate_2D(
      const harm_aux_vars_struct *restrict harm_aux,
      const double dx[2],
      const double x_old[2],
      double x[2]) {
  const double dv = 1.e-15;

  x[0] = fabs(x[0]);
  x[0] = (x[0] > 1e20) ? x_old[0] : x[0];

  x[1] = (x[1] < 0.) ? 0. : x[1];        /* if it's too small */
  x[1] = (x[1] > 1.) ? (1. - dv) : x[1]; /* if it's too big   */
}
