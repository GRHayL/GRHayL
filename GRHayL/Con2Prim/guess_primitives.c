#include "ghl_con2prim.h"

/* Function     : ghl_guess_primitives()
 * Description  : Computes initial guesses for the primitives
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_guess_primitives
*/

void ghl_guess_primitives(
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims) {

  // Use atmosphere as initial guess:
  prims->rho         = cons_undens->rho;
  prims->u0          = 1.0;
  prims->vU[0]       = -ADM_metric->betaU[0];
  prims->vU[1]       = -ADM_metric->betaU[1];
  prims->vU[2]       = -ADM_metric->betaU[2];
  prims->Y_e         = cons_undens->Y_e/cons_undens->rho;
  prims->temperature = eos->T_max;

  // If using hybrid or ideal fluid EOS, compute P_cold and eps_cold
  if( eos->eos_type == ghl_eos_hybrid || eos->eos_type == ghl_eos_simple )
    ghl_hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &prims->press, &prims->eps);
}
