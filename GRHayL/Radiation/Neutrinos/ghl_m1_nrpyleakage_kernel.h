#ifndef GHL_M1_NRPYLEAKAGE_KERNEL_H_
#define GHL_M1_NRPYLEAKAGE_KERNEL_H_

#include "ghl_radiation.h"
#include "ghl_eos_functions.h"

/*
 * Private M1 adapter types.  These intentionally do not extend the public
 * NRPyLeakage API: legacy leakage keeps its original source and header
 * contract, while the M1 provider consumes channel-resolved rates here.
 */
typedef enum {
  ghl_m1_nrpyleakage_nue = 0,
  ghl_m1_nrpyleakage_anue = 1,
  ghl_m1_nrpyleakage_nux = 2,
  ghl_m1_nrpyleakage_species_count = 3
} ghl_m1_nrpyleakage_species_t;

typedef struct {
  double rho;
  double rho_cgs;
  double T;
  double Ye;
  double muhat;
  double mu_e;
  double mu_p;
  double mu_n;
  double X_n;
  double X_p;
  double Y_np;
  double Y_pn;
} ghl_m1_nrpyleakage_thermo_state;

typedef struct {
  double neutrino_degeneracy;
  double F2;
  double F3;
  double F4;
  double F5;

  double n_eq_cgs;
  double J_eq_mev_cgs;
  double mean_energy_mev;

  double eta_N_beta_cgs;
  double eta_N_pair_cgs;
  double eta_N_plasmon_cgs;
  double eta_N_brems_cgs;

  double eta_E_beta_mev_cgs;
  double eta_E_pair_mev_cgs;
  double eta_E_plasmon_mev_cgs;
  double eta_E_brems_mev_cgs;

  double kappa_a_N_cc_cgs;
  double kappa_a_E_cc_cgs;

  double kappa_s_N_neutron_cgs;
  double kappa_s_N_proton_cgs;
  double kappa_s_E_neutron_cgs;
  double kappa_s_E_proton_cgs;
} ghl_m1_nrpyleakage_species_raw_rates;

typedef struct {
  ghl_m1_nrpyleakage_species_raw_rates
      species[ghl_m1_nrpyleakage_species_count];
  int nux_single_species_multiplicity;
} ghl_m1_nrpyleakage_raw_rates;

typedef struct {
  ghl_m1_nrpyleakage_raw_rates raw;
  double eta_N[ghl_m1_nrpyleakage_species_count];
  double eta_E[ghl_m1_nrpyleakage_species_count];
  double kappa_N[ghl_m1_nrpyleakage_species_count];
  double kappa_E[ghl_m1_nrpyleakage_species_count];
  double n_eq[ghl_m1_nrpyleakage_species_count];
  double J_eq[ghl_m1_nrpyleakage_species_count];
} ghl_m1_nrpyleakage_legacy_kernel_result;

ghl_error_codes_t ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
      double rho,
      double Ye,
      double T,
      double muhat,
      double mu_e,
      double mu_p,
      double mu_n,
      double X_n,
      double X_p,
      bool strict_validation,
      ghl_m1_nrpyleakage_thermo_state *restrict thermo);

ghl_error_codes_t ghl_m1_nrpyleakage_compute_thermo_state(
      const ghl_eos_parameters *restrict eos,
      double rho,
      double Ye,
      double T,
      ghl_m1_nrpyleakage_thermo_state *restrict thermo);

ghl_error_codes_t ghl_m1_nrpyleakage_compute_thermo_state_legacy(
      const ghl_eos_parameters *restrict eos,
      double rho,
      double Ye,
      double T,
      ghl_m1_nrpyleakage_thermo_state *restrict thermo);

ghl_error_codes_t ghl_m1_nrpyleakage_compute_legacy_kernel_from_thermo(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double neutrino_degeneracy[ghl_m1_nrpyleakage_species_count],
      ghl_m1_nrpyleakage_legacy_kernel_result *restrict result);

ghl_error_codes_t ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_legacy(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double neutrino_degeneracy[ghl_m1_nrpyleakage_species_count],
      ghl_m1_nrpyleakage_raw_rates *restrict raw);

ghl_error_codes_t ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double neutrino_degeneracy[ghl_m1_nrpyleakage_species_count],
      ghl_m1_nrpyleakage_raw_rates *restrict raw);

#endif

