#ifndef GHL_RUSANOV_PRIVATE_H_
#define GHL_RUSANOV_PRIVATE_H_

#include "ghl.h"
#include "ghl_m1.h"

/* ghl_calculate_Rusanov_flux() is the public component-wise helper; its
 * declaration lives with the installed M1 transport surface in ghl_m1.h. */

/* Scalar arithmetic shared by the undensitized wrapper and the
 * prepared-face M1 wrapper. Operands must be finite; callers own
 * validation and output publication. */
double ghl_rusanov_average_finite(double left, double right);

double ghl_rusanov_candidate_finite(
      double state_L,
      double state_R,
      double physical_flux_L,
      double physical_flux_R,
      double speed);

#endif // GHL_RUSANOV_PRIVATE_H_
