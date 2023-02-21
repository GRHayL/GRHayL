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
                       const primitive_quantities *restrict prims,
                       const conservative_quantities *restrict cons,
                       primitive_quantities *restrict prims_guess ) {

  //Use atmosphere as initial guess:
  prims_guess->rho = cons->rho/metric->psi6;
  prims_guess->vx = -metric->betax;
  prims_guess->vy = -metric->betay;
  prims_guess->vz = -metric->betaz;
  eos->hybrid_compute_P_cold(eos, prims_guess->rho, &prims_guess->press);
  prims_guess->Y_e = cons->Y_e/cons->rho;
  prims_guess->temperature = eos->T_atm;
}
