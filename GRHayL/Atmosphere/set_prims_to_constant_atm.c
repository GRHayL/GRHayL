#include "ghl_atmosphere.h"

/*
 * Function      : ghl_set_prims_to_constant_atm()
 * Description   : Uses the EOS data to reset the primitives to
 *                 constant density atmosphere
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_set_prims_to_constant_atm
*/

GHL_DEVICE
void ghl_set_prims_to_constant_atm(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims) {

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

  prims->vU[0] = 0.0;
  prims->vU[1] = 0.0;
  prims->vU[2] = 0.0;
}
