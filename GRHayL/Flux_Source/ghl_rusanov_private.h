#ifndef GHL_RUSANOV_PRIVATE_H_
#define GHL_RUSANOV_PRIVATE_H_

/* Scalar arithmetic shared by the public undensitized wrapper and the
 * private prepared-face M1 wrapper. Operands must be finite; callers own
 * validation and output publication. */
double ghl_rusanov_average_finite(double left, double right);

double ghl_rusanov_candidate_finite(
      double state_L,
      double state_R,
      double physical_flux_L,
      double physical_flux_R,
      double speed);

#endif
