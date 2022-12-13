#include "con2prim_gem.h"

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

void reset_prims_to_atmosphere( const GRHayL_parameters *restrict params,
                                const eos_parameters *restrict eos,
                                const metric_quantities *restrict metric,
                                primitive_quantities *restrict prims,
                                con2prim_diagnostics *restrict diagnostics ) {

  // Just a simple reset to atmospheric values.
  // Velocities are set to zero. Keeping it
  // inside a single function ensures that
  // resets are consistent throughout the code.
  prims->rho         = eos->rho_atm;
  prims->press       = eos->press_atm;
  prims->eps         = eos->eps_atm;
  prims->entropy     = eos->entropy_atm;
  prims->Y_e         = eos->Ye_atm;
  prims->temperature = eos->T_atm;
  
  if(params->Cupp_Fix) { //For atmosphere, ET IGM sets v=-beta, Leo's IGM sets v=0
    prims->vx = 0.0;
    prims->vy = 0.0;
    prims->vz = 0.0;
  } else {
    prims->vx = -metric->betax;
    prims->vy = -metric->betay;
    prims->vz = -metric->betaz;
  } 
}
