#include "../../utils_Noble.h"

void ghl_validate_1D(
      const double x0[1],
      double x[1]) {
  x[0] = fabs(x[0]);
}

void ghl_validate_2D(
      const double x0[2],
      double x[2]) {
  const double dv = 1.e-15;

  x[0] = fabs(x[0]);
  x[0] = (x[0] > Z_TOO_BIG) ?  x0[0] : x[0];

  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] > 1.) ?  (1. - dv) : x[1];  /* if it's too big   */
}

