#include "con2prim_header.h"

/* Function    : reset_prims_to_atmosphere()
 * Authors     : ? and Samuel Cupp
 * Description : Uses the EOS data to reset the primitives to atmospheric
 *               values.
 *
 * Inputs      : eos             - an initialized eos_parameters struct
 *                                 with data for the EOS of the simulation
 *
 * Outputs     : prims           - A primitive_quantities struct which has
 *                                 been floored to atmospheric values
 */

void reset_prims_to_atmosphere( const eos_parameters *restrict eos,
                                primitive_quantities *restrict prims,
                                con2prim_diagnostics *restrict diagnostics ) {

  // Just a simple reset to atmospheric values.
  // Velocities are set to zero. Keeping it
  // inside a single function ensures that
  // resets are consistent throughout the code.
  prims->rho = eos->rho_atm;
  prims->press = eos->press_atm;
  prims->eps = eos->eps_atm;
  prims->entropy = eos->entropy_atm;
  prims->Y_e = eos->Ye_atm;
  prims->temp = eos->temp_atm;
  
  prims->vx = 0.0;
  prims->vy = 0.0;
  prims->vz = 0.0;  
}
