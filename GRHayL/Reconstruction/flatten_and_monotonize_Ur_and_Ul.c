#include "reconstruction.h"

void ghl_flatten_and_monotonize_Ur_and_Ul(
      const double U,
      const double ftilde,
      double *restrict Ur,
      double *restrict Ul) {

  // First detect shocks / steep gradients:
  *Ur = U*ftilde + (*Ur)*(1.0 - ftilde);
  *Ul = U*ftilde + (*Ul)*(1.0 - ftilde);

  // Then monotonize all variables
  const double dU = (*Ur) - (*Ul);
  const double Utmp = dU*( U - 0.5*((*Ur) + (*Ul)) );

  if ( ((*Ur) - U)*(U - (*Ul)) <= 0.0) {
    (*Ur) = U;
    (*Ul) = U;
    return;
  }
  if ( Utmp > (1.0/6.0)*(dU*dU)) {
    (*Ul) = 3.0*U - 2.0*(*Ur);
    return;
  }
  if ( Utmp < -(1.0/6.0)*(dU*dU)) {
    (*Ur) = 3.0*U - 2.0*(*Ul);
    return;
  }
}
