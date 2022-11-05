#include "con2prim_header.h"

/* This function fills the struct GRMHD_parameters with data.
   For more information on the arguments, see the definition of the
   struct in new_header.h. */
void initialize_parameters(GRMHD_parameters *restrict params, const int main, const int backup[3],
                const int evolve_entropy, const int evolve_temp, const int calc_prim_guess,
                const double gamma_speed_limit, const double psi6threshold, const int update_Tmunu) {
  params->main_routine = main;
  params->backup_routine[0] = backup[0];
  params->backup_routine[1] = backup[1];
  params->backup_routine[2] = backup[2];
  params->evolve_entropy = evolve_entropy;
  params->evolve_temp = evolve_temp;
  params->calc_prim_guess = calc_prim_guess;
  params->gamma_speed_limit = gamma_speed_limit;
  params->psi6threshold = psi6threshold;
  params->update_Tmunu = update_Tmunu;
  return;
}
