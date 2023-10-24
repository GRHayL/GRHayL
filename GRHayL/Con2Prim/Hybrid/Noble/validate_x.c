#include "../../utils_Noble.h"

void ghl_validate_1D(
      const harm_aux_vars_struct *restrict harm_aux,
      const double dx[1],
      const double x_old[1],
      double x[1]) {

  // This is only in the Noble1D version of HARM's Newton-Raphson solver.
  // It is unclear what exactly this does, but it is necessary for convergence.
  double i_increase = 0;
  while( (( x[0]*x[0]*x[0] * ( x[0] + 2.*harm_aux->Bsq ) -
            harm_aux->QdotBsq*(2.*x[0] + harm_aux->Bsq) ) <= x[0]*x[0]*(harm_aux->Qtsq-harm_aux->Bsq*harm_aux->Bsq))
         && (i_increase < 10) ) {
    x[0] -= (1.*i_increase) * dx[0] / 10. ;
    i_increase++;
  }

  x[0] = fabs(x[0]);
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
  x[0] = (x[0] > Z_TOO_BIG) ?  x_old[0] : x[0];

  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] > 1.) ?  (1. - dv) : x[1];  /* if it's too big   */
}

