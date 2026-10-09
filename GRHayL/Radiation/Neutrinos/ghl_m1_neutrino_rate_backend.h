#ifndef GHL_M1_NEUTRINO_RATE_BACKEND_H_
#define GHL_M1_NEUTRINO_RATE_BACKEND_H_

#include "ghl_radiation.h"

/* Private boundary between the provider's shared cache/orchestration and
 * tabulated-EOS plus NRPyLeakage physics. Disabled-HDF5 builds provide the
 * same callbacks from a separate stub branch in the implementation TU. */
ghl_error_codes_t ghl_m1_neutrino_rate_backend_initialize(void);

ghl_error_codes_t ghl_m1_neutrino_rate_backend_temperature_from_eps(
      const ghl_eos_parameters *restrict eos,
      double rho,
      double Ye,
      double eps,
      double *restrict T);

ghl_error_codes_t ghl_m1_neutrino_rate_backend_thermo_from_T(
      const ghl_eos_parameters *restrict eos,
      double rho,
      double Ye,
      double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p);

ghl_error_codes_t ghl_m1_neutrino_rate_backend_assemble(
      int channel_mask,
      double nu_x_multiplicity,
      double rho,
      double T,
      double Ye,
      double muhat,
      double mu_e,
      double mu_p,
      double mu_n,
      double X_n,
      double X_p,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count],
      double beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_species_count],
      bool beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_species_count]);

#endif
