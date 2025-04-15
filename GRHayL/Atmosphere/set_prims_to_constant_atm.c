#include "ghl_atmosphere.h"

/**
 * @ingroup Atmosphere
 * @brief Set the primitives to a constant density atmosphere
 * @details
 * This function sets the given primitive struct to a constant
 * density atmosphere using the provided EOS information. Velocities
 * are set to zero.
 *
 * @param[in] eos:     pointer to ghl_eos_parameters struct.
 *
 * @param[out] params: pointer to ghl_primitive_quantities struct.
*/

void ghl_set_prims_to_constant_atm(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims) {

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
