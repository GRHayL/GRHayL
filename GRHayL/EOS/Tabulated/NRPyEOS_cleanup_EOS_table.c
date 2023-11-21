#include "nrpyeos_tabulated.h"

/// @brief Adjust and clean up the table.
/// @param eos Struct with equation of state parameters.
/// @details This function performs the following tasks: adjust the reference
/// mass of the table s.t. energy_shift = 0; ensure the sound speed is
/// relativistic; and enforce a floor and ceiling to the sound speed squared to
/// zero and one, respectively.
void NRPyEOS_cleanup_EOS_table(ghl_eos_parameters *restrict eos) {

  ghl_info("Adjusting and cleaning the EOS table...\n");

  const int npoints = eos->N_rho * eos->N_T * eos->N_Y_e;
  const double delta = 1.1 * fabs(eos->table_eps_min);

  for(int i = 0; i < npoints; i++) {

    const int press_idx = NRPyEOS_press_key + NRPyEOS_ntablekeys * i;
    const int cs2_idx = NRPyEOS_cs2_key + NRPyEOS_ntablekeys * i;

    // Read from the table
    double rho = exp(eos->table_logrho[i]);
    double eps = exp(eos->table_logeps[i]) + eos->energy_shift;
    double cs2 = eos->table_all[cs2_idx];
    const double press = eos->table_all[press_idx];

    // Adjust rho and eps, if needed
    if(eos->table_epsmin < 0.0) {
      // Compute the total internal energy: e = (1 + eps)*rho
      const double e = (1.0 + eps) * rho;

      // Adjust eps
      eps += delta;

      // Recompute rho: rho = e / (1 + eps)
      rho = e / (1.0 + eps);
    }

    // Adjust cs2, if needed
    if(*cs2 > 1.0 || isnan(*cs2) || isinf(*cs2)) {
      *cs2 = 1.0;
    }
    else if(*cs2 < 0.0) {
      *cs2 = 0.0;
    }

    // Set cs2 to the relativistic sound speed, if needed
    if(!eos->cs2_is_relativistic) {
      const double h = 1.0 + eps + press / rho;
      cs2 /= h;
    }

    // Write to the table
    eos->table_logrho[i] = log(rho);
    eos->table_logeps[i] = log(eps);
    eos->table_all[cs2_idx] = cs2;
  }

  // Adjust epsmin: epsmin -> epsmin + delta
  if(eos->table_eps_min < 0.0) {
    eos->table_epsmin = eos->table_epsmin + delta;
  }

  // Adjust energy_shift to zero
  if(eos->energy_shift != 0.0) {
    eos->energy_shift = 0.0;
  }

  ghl_info("Finished adjusting and cleaning the EOS table\n");
}
