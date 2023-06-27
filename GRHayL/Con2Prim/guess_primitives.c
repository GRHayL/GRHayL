#include "con2prim.h"

/* Function    : ghl_guess_primitives()
 * Description : Computes initial guesses for the primitives
 *               by assuming the values are atmosphere
 *
 * Inputs      : eos            - eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : cons           - conservative_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : prims_guess    - returns initial guesses for primitives
 */

void ghl_guess_primitives(
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims) {

  // Use atmosphere as initial guess:
  prims->rho = cons->rho/ADM_metric->sqrt_detgamma;
  prims->vU[0] = -ADM_metric->betaU[0];
  prims->vU[1] = -ADM_metric->betaU[1];
  prims->vU[2] = -ADM_metric->betaU[2];
  prims->Y_e = cons->Y_e/cons->rho;
  prims->temperature = eos->T_max;

  // If using Hybrid EOS, compute P_cold
  if( eos->eos_type == ghl_eos_hybrid )
    ghl_hybrid_compute_P_cold(eos, prims->rho, &prims->press);
}
