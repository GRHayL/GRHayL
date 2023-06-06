#include "grhayl.h"

/* Function    : ghl_initialize_params()
 * Description : Initialize the grhayl_parameters struct from user input
 *
 * Inputs      : main            - selection of Con2Prim routine for the
 *                                 ghl_con2prim_multi_method() function
 *             : backup[3]       - array for up to three backup routines for
 *                                 the ghl_con2prim_multi_method() function
 *             : evolve_entropy  - sets whether entropy is being evolved in
 *                                 the simulation
 *             : evolve_temp     - sets whether temperature is being evolved in
 *                                 the simulation
 *             : calc_prim_guess - sets whether ghl_con2prim_multi_method() should compute
 *                                 initial guesses for the primitives; if set to 0,
 *                                 the initial guess will be whatever has been loaded into
 *                                 the primitive_quantities struct passed to ghl_con2prim_multi_method()
 *             : psi6threshold   - TODO: comment
 *             : Cupp_Fix        - If active, will change behavior of Noble2D which removes
 *                                 cause of most Font fixes and also change the choice of velocity
 *                                 in the case of atmospheric reset
 *
 * Outputs     : params          - fully initialized grhayl_parameters struct
 *                                 containing the input parameters
 */

void ghl_initialize_params(
      const int main_routine,
      const int backup_routine[3],
      const bool evolve_entropy,
      const bool evolve_temp,
      const bool calc_prim_guess,
      const double psi6threshold,
      const bool Cupp_Fix,
      const double Lorenz_damping_factor,
      grhayl_parameters *restrict params ) {

  // This comment is new
  params->main_routine          = main_routine;
  params->backup_routine[0]     = backup_routine[0];
  params->backup_routine[1]     = backup_routine[1];
  params->backup_routine[2]     = backup_routine[2];
  params->evolve_entropy        = evolve_entropy;
  params->evolve_temp           = evolve_temp;
  params->calc_prim_guess       = calc_prim_guess;
  params->psi6threshold         = psi6threshold;
  params->Cupp_Fix              = Cupp_Fix;
  params->Lorenz_damping_factor = Lorenz_damping_factor;
}
