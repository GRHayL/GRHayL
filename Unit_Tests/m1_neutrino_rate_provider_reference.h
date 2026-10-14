#ifndef M1_NEUTRINO_RATE_PROVIDER_REFERENCE_H
#define M1_NEUTRINO_RATE_PROVIDER_REFERENCE_H

/* Test-local synthetic context and cache. The test translation unit includes
 * the shared GRHayL rate and diagnostic types before this header. */
typedef struct {
  int channel_mask;
  ghl_neutrino_rate_failure_policy_t failure_policy;
  ghl_neutrino_rate_table_bounds_policy_t table_bounds_policy;
  double nu_x_multiplicity;
  uint64_t eos_generation;
  double temperature_code_to_mev;
  double baryon_mass_code;
  double charged_current_scale;
  double scattering_scale;
  double pair_scale;
  double bremsstrahlung_scale;
  double plasmon_scale;
  double min_mean_energy;
  double equilibrium_recovery_rate;
} m1_test_reference_provider_context;

typedef struct {
  bool thermo_valid;
  bool rates_valid;
  double thermo_rho;
  double thermo_T;
  double thermo_Ye;
  double rho;
  double T;
  double Ye;
  double muhat;
  double mu_e;
  double mu_p;
  double mu_n;
  double X_n;
  double X_p;
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  m1_test_reference_provider_context provider_snapshot;
  const ghl_eos_parameters *eos_snapshot;
  double beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_species_count];
  bool beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_species_count];
} m1_test_reference_provider_cache;

ghl_error_codes_t m1_test_reference_provider_initialize(
      m1_test_reference_provider_context *restrict provider);
void m1_test_reference_provider_cache_initialize(
      m1_test_reference_provider_cache *restrict cache);
ghl_error_codes_t m1_test_reference_provider_compute_cell(
      const m1_test_reference_provider_context *restrict provider,
      m1_test_reference_provider_cache *restrict cache,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict unused_eos,
      const ghl_primitive_quantities *restrict prims,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]);

#endif
