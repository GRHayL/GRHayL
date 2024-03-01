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
      const double max_Lorentz_factor,
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
  params->max_Lorentz_factor         = max_Lorentz_factor;
  params->inv_sq_max_Lorentz_factor  = 1.0/SQR(max_Lorentz_factor);
  params->Lorenz_damping_factor     = Lorenz_damping_factor;

  // Initialize default Con2Prim values
  params->con2prim_max_iterations   = 30;
  params->con2prim_solver_tolerance = 1e-10;

  // Initialize default PPM values (from Colella & Woodward)
  params->ppm_flattening_epsilon = 0.33;
  params->ppm_flattening_omega1  = 0.75;
  params->ppm_flattening_omega2  = 10.0;
  params->ppm_shock_k0           = 0.1;
  params->ppm_shock_eta1         = 20.0;
  params->ppm_shock_eta2         = 0.05;
  params->ppm_shock_epsilon      = 0.01;
}
