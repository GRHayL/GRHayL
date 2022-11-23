// Thorn      : IllinoisGRMHD
// File       : con2prim_set_cons_and_prim_from_CONSERVS_and_PRIMS.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This provides functions which 1. convert IllinoisGRMHD's set
//              of conservative variables into the appropriate variables
//              required by the C2P routines and 2. set appropriate primitive
//              guesses.

#include "con2prim_header.h"

void undensitize_conservatives( const eos_parameters *restrict eos,
                                const int c2p_key,
                                const metric_quantities *restrict metric,
                                const primitive_quantities *restrict prims,
                                const conservative_quantities *restrict cons,
                                conservative_quantities *restrict cons_undens ) {

  // (TODO: library name)'s variable is the "densitized" version of
  // the standard conservative variables (D,tau,S_{i}). In
  // other words, we have the relationships:
  //
  // rho_star   = sqrt(gamma) *  D
  // tilde(tau) = sqrt(gamma) * tau
  // tilde(S)_i = sqrt(gamma) * S_i
  //
  // Therefore the conversion between the two is straightfoward.

  const double psim6 = 1.0/metric->psi6;

  cons_undens->rho = cons->rho * psim6;
  cons_undens->S_x = cons->S_x * psim6;
  cons_undens->S_y = cons->S_y * psim6;
  cons_undens->S_z = cons->S_z * psim6;
  cons_undens->tau = cons->tau * psim6;
  cons_undens->Y_e = cons->Y_e * psim6;
  cons_undens->entropy = cons->entropy * psim6;
}
