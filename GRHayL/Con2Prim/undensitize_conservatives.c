#include "ghl_con2prim.h"

/* Function     : ghl_undensitize_conservatives()
 * Description  : Computes undensitized conservatives using the metric
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_undensitize_conservatives
*/

void ghl_undensitize_conservatives(
      const double psi6,
      const ghl_conservative_quantities *restrict cons,
      ghl_conservative_quantities *restrict cons_undens) {

  /*
     Many codes evolve the"densitized" conservative
     variables (D,tau,S_{i}). We have the relationships:

     rho_star   = sqrt(gamma) *  D
     tilde(tau) = sqrt(gamma) * tau
     tilde(S)_i = sqrt(gamma) * S_i
   */

  const double psim6 = 1.0/psi6;

  cons_undens->rho     = cons->rho * psim6;
  cons_undens->SD[0]   = cons->SD[0] * psim6;
  cons_undens->SD[1]   = cons->SD[1] * psim6;
  cons_undens->SD[2]   = cons->SD[2] * psim6;
  cons_undens->tau     = cons->tau * psim6;
  cons_undens->Y_e     = cons->Y_e * psim6;
  cons_undens->entropy = cons->entropy * psim6;
}
