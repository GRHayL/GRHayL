#include "GRHayL.h"

/* Function    : initialize_GRHayL()
 * Authors     : Samuel Cupp
 * Description : Initialize the GRHayL_parameters struct from user
 *               input
 *
 * Inputs      : main             - selection of Con2Prim routine for the
 *                                  C2P_Select_Hybrid_Method() function
 *             : backup[3]        - array for up to three backup routines for
 *                                  the C2P_Select_Hybrid_Method() function
 *             : evolve_entropy   - 
 *             : evolve_temp      - 
 *             : calc_prim_guess  - Codes which have valid primitive data at
 *                                  the start of Con2Prim routines do not need
 *                                  to generate guesses and can immediately run
 *                                  C2P. Those that don't need to calculate primitive
 *                                  guesses. Setting this to 1 triggers the calculation
 *                                  of guesses. TODO: idk if this is how we want to handle this
 *             : psi6threshold    - 
 *             : update_Tmunu     - This will update the T_{\mu\nu} values of the stress_energy
 *                                  struct after obtaining the new prims and cons
 *
 * Outputs     : params           - fully initialized GRHayL_parameters struct
 *                                  containing the input parameters
 */

void initialize_GRHayL(const int main, const int backup[3],
                const int evolve_entropy, const int evolve_temp, const int calc_prim_guess,
                const double psi6threshold, const int update_Tmunu, const int Cupp_Fix,
                GRHayL_parameters *restrict params) {
  params->main_routine = main;
  params->backup_routine[0] = backup[0];
  params->backup_routine[1] = backup[1];
  params->backup_routine[2] = backup[2];
  params->evolve_entropy = evolve_entropy;
  params->evolve_temp = evolve_temp;
  params->calc_prim_guess = calc_prim_guess;
  params->psi6threshold = psi6threshold;
  params->update_Tmunu = update_Tmunu;
  params->Cupp_Fix = Cupp_Fix;
  return;
}
