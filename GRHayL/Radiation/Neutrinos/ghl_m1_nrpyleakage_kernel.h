#ifndef GHL_M1_NRPYLEAKAGE_KERNEL_H_
#define GHL_M1_NRPYLEAKAGE_KERNEL_H_

#include "ghl.h"
#include "ghl_radiation.h"

/*
 * Private M1 adapter types.  These intentionally do not extend the public
 * NRPyLeakage API: legacy leakage keeps its original source and header
 * contract, while the M1 provider consumes channel-resolved rates here.
 * The primitive thermodynamic state is kept separate from density-derived
 * NRPyLeakage blocking populations; the latter are evaluated at the raw-rate
 * boundary so callers cannot publish stale derived values.
 */
typedef enum {
  ghl_m1_nrpyleakage_nue = 0,
  ghl_m1_nrpyleakage_anue = 1,
  ghl_m1_nrpyleakage_nux = 2,
  ghl_m1_nrpyleakage_species_count = 3
} ghl_m1_nrpyleakage_species_t;

typedef struct {
  double rho;
  double T;
  double Ye;
  double muhat;
  double mu_e;
  double mu_p;
  double mu_n;
  double X_n;
  double X_p;
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
  ghl_m1_nrpyleakage_species_raw_rates species[ghl_m1_nrpyleakage_species_count];
  int nux_single_species_multiplicity;
} ghl_m1_nrpyleakage_raw_rates;

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
      ghl_m1_nrpyleakage_thermo_state *restrict thermo);

ghl_error_codes_t ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double neutrino_degeneracy[ghl_m1_nrpyleakage_species_count],
      ghl_m1_nrpyleakage_raw_rates *restrict raw);

#endif
