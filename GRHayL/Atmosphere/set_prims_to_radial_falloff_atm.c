#include "con2prim.h"

/* Function    : grhayl_set_prims_to_radial_falloff_atm()
 * Description : Uses the EOS data to reset the primitives to atmospheric
 *               values.
 *
 * Inputs      : eos            - eos_parameters struct with data for the
 *                                EOS of the simulation
 *
 * Outputs     : prims          - returns with all primitives set to atmospheric values
 */

void grhayl_set_prims_to_radial_falloff_atm(
      const eos_parameters *restrict eos,
      const double r, //not sure what is actually needed
      primitive_quantities *restrict prims) {

//const double r_atmo = isotropic_atmo ? 1.0 : std::max(HARM3D_atmo_min_r,r[ijk]);
//
//const double rho_floor   = atmo_rho * std::pow(r_atmo,-1.5);
//const double uu_floor    = atmo_uu  * std::pow(r_atmo,-1.5*gl_gamma);
//const double press_floor = uu_floor * (gl_gamma - 1.0);
//
//rho[ijk]   = rho_floor;
//press[ijk] = press_floor;
//
//  prims->rho         = ;
//  prims->press       = ;
//  prims->eps         = press->press / prims->rho / (gl_gamma - 1.0);
  prims->entropy     = eos->entropy_atm;
  prims->Y_e         = eos->Y_e_atm;
  prims->temperature = eos->T_atm;

  prims->vU[0] = 0.0;
  prims->vU[1] = 0.0;
  prims->vU[2] = 0.0;
}
