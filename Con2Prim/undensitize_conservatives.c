#include "con2prim.h"

/* Function    : undensitize_conservatives()
 * Authors     : Leo Werneck and Samuel Cupp
 * Description : this function takes a conservative_quantities struct cons with
 *               densitized variables and computes the undensitized variables,
 *               which are stored in the struct cons_undens
 * Dependencies: None
 *
 * Inputs      : metric         - an initialized metric_quantities struct
 *                                with data for the same gridpoint as cons
 * Inputs      : cons           - an initialized conservative_quantities
 *                                struct with data for the same gridpoint as
 *                                metric
 *
 * Outputs     : cons_undens    - conservative_quantities struct filled with
 *                                the undensitized variables associated with
 *                                cons
 */

void undensitize_conservatives( const metric_quantities *restrict metric,
                                const conservative_quantities *restrict cons,
                                conservative_quantities *restrict cons_undens ) {

  // (TODO: library name)'s variables are the "densitized" versions
  // of the standard conservative variables (D,tau,S_{i}). In
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
