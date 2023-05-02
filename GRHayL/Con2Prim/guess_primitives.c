#include "con2prim.h"

/* Function    : guess_primitives()
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

void guess_primitives( const eos_parameters *restrict eos,
                       const metric_quantities *restrict metric,
                       const conservative_quantities *restrict cons,
                       primitive_quantities *restrict prims ) {

  // Use atmosphere as initial guess:
  prims->rho = cons->rho/metric->psi6;
  prims->vx = -metric->betax;
  prims->vy = -metric->betay;
  prims->vz = -metric->betaz;
  prims->Y_e = cons->Y_e/cons->rho;
  prims->temperature = eos->T_atm;

  // If using Hybrid EOS, compute P_cold
  if( eos->eos_type == grhayl_eos_hybrid )
    eos->hybrid_compute_P_cold(eos, prims->rho, &prims->press);
}
