// Thorn      : IllinoisGRMHD
// File       : con2prim_set_cons_and_prim_from_CONSERVS_and_PRIMS.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This provides functions which 1. convert IllinoisGRMHD's set
//              of conservative variables into the appropriate variables
//              required by the C2P routines and 2. set appropriate primitive
//              guesses.

#include "con2prim_header.h"

void undensitize( const eos_parameters *restrict eos,
                  const int c2p_key,
                  const metric_quantities *restrict metric,
                  const primitive_quantities *restrict prims,
                  const conservative_quantities *restrict cons,
                  conservative_quantities *restrict cons_undens ) {

  // IllinoisGRMHD's variable is the "densitised" version of
  // the standard conservative variables (D,tau,S_{i}). In
  // other words, we have the relationships:
  //
  // rho_star   = sqrt(gamma) *  D
  // tilde(tau) = sqrt(gamma) * tau
  // tilde(S)_i = sqrt(gamma) * S_i
  //
  // Therefore the conversion between the two is straightfoward.
  //
  // One important note, however, is that the Noble2D con2prim
  // routine does not use the variable tau, but instead the
  // energy variable u which is related to IllinoisGRMHD's
  // conservatives via the relation:
  //
  // u = -alpha*tau - (alpha-1)*rho_star + beta^{i}tilde(S)_{i}
  //
  // The magnetic fields in IllinoisGRMHD also need to be
  // rescaled by a factor of sqrt(4pi). In the case of
  // the Noble2D routine, an extra factor of the lapse
  // is also required.
  //
  // Now let us begin the conversion.

  // First compute alpha*sqrt(gamma) = alpha*psi^(6)

  // Finally compute the remaining quantities (which are routine specific)

  // Note that in the c2p routine we used to
  // multiply the cons vector by alpha/detg.
  // But this is equivalent to:
  //
  // alpha/detg = alpha/( alpha*sqrt(gamma) )
  //            = alpha/( alpha*psi^6 )
  //            = psi^{-6}
  const double psim6 = 1.0/metric->psi6;

  cons_undens->rho = cons->rho * psim6;
  cons_undens->S_x = cons->S_x * psim6;
  cons_undens->S_y = cons->S_y * psim6;
  cons_undens->S_z = cons->S_z * psim6;
  cons_undens->tau = cons->tau * psim6;
  cons_undens->Y_e = cons->Y_e * psim6;
  cons_undens->entropy = cons->entropy * psim6;
}
