#include "con2prim_header.h"

/* Function    : ()
 * Authors     : ? and Samuel Cupp
 * Description : Initialize the conservative struct from user
 *               input
 * Dependencies: None
 *
 * Inputs      : rho             - value of rho_star (densitized density)
 *             : tau             - value of tau tilde (densitized energy
 *                                 variable)
 *             : S_x             - value of the x component of S tilde
 *                                 (densitized momentum variable)
 *             : S_y             - value of the y component of S tilde
 *                                 (densitized momentum variable)
 *             : S_z             - value of the z component of S tilde
 *                                 (densitized momentum variable)
 *             : entropy         - value of densitized entropy
 *
 * Outputs     : cons            - fully initialized conservative_quantities
 *                                 struct containing the input data
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
