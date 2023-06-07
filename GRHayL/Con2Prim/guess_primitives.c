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
      const metric_quantities *restrict metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons,
      primitive_quantities *restrict prims) {

  // Use atmosphere as initial guess:
  prims->rho = cons->rho/metric_aux->psi6;
  prims->vU[0] = -metric->betaU[0];
  prims->vU[1] = -metric->betaU[1];
  prims->vU[2] = -metric->betaU[2];
  prims->Y_e = cons->Y_e/cons->rho;
  prims->temperature = eos->T_max;

  // If using Hybrid EOS, compute P_cold
  if( eos->eos_type == grhayl_eos_hybrid )
    ghl_hybrid_compute_P_cold(eos, prims->rho, &prims->press);
}
