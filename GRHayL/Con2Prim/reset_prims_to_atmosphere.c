#include "con2prim.h"

/* Function    : reset_prims_to_atmosphere()
 * Description : Uses the EOS data to reset the primitives to atmospheric
 *               values.
 *
 * Inputs      : eos            - eos_parameters struct with data for the
 *                                EOS of the simulation
 *
 * Outputs     : prims          - returns with all primitives set to atmospheric values
 */

void reset_prims_to_atmosphere(
      const eos_parameters *restrict eos,
      primitive_quantities *restrict prims) {

  // Just a simple reset to atmospheric values.
  // Velocities are set to zero. Keeping it
  // inside a single function ensures that
  // resets are consistent throughout the code.
  prims->rho         = eos->rho_atm;
  prims->press       = eos->press_atm;
  prims->eps         = eos->eps_atm;
  prims->entropy     = eos->entropy_atm;
  prims->Y_e         = eos->Y_e_atm;
  prims->temperature = eos->T_atm;

  prims->vx = 0.0;
  prims->vy = 0.0;
  prims->vz = 0.0;
}
