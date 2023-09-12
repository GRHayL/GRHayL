#include "ghl.h"

/*
 * Function     : ghl_initialize_params()
 * Description  : Initialize the ghl_parameters struct from user input
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_initialize_params
*/

void ghl_initialize_params(
      const ghl_con2prim_method_t main_routine,
      const ghl_con2prim_method_t backup_routine[3],
      const bool evolve_entropy,
      const bool evolve_temp,
      const bool calc_prim_guess,
      const double psi6threshold,
      const bool Cupp_Fix,
      const double max_lorenz_factor,
      const double Lorenz_damping_factor,
      ghl_parameters *restrict params) {

  params->main_routine              = main_routine;
  params->backup_routine[0]         = backup_routine[0];
  params->backup_routine[1]         = backup_routine[1];
  params->backup_routine[2]         = backup_routine[2];
  params->evolve_entropy            = evolve_entropy;
  params->evolve_temp               = evolve_temp;
  params->calc_prim_guess           = calc_prim_guess;
  params->psi6threshold             = psi6threshold;
  params->Cupp_Fix                  = Cupp_Fix;
  params->max_lorenz_factor         = max_lorenz_factor;
  params->inv_sq_max_lorenz_factor  = 1.0/SQR(max_lorenz_factor);
  params->Lorenz_damping_factor     = Lorenz_damping_factor;
  params->con2prim_max_iterations   = 30;
  params->con2prim_solver_tolerance = 1e-10;
}
