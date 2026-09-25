#ifndef GHL_M1_NRPYLEAKAGE_KERNEL_H_
#define GHL_M1_NRPYLEAKAGE_KERNEL_H_

#include <float.h>
#include <math.h>

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

/* Current NRPyLeakage keeps this interpolation-boundary normalization inside
 * its blocking evaluator.  The M1 adapter also needs the normalized fractions
 * before caching them and before using them in its separate raw-rate algebra,
 * so retain the same conservative 27*gamma_64 acceptance envelope here.
 * Callers validate finiteness first and supply local output addresses. */
static inline ghl_error_codes_t ghl_m1_nrpyleakage_normalize_nucleon_fractions(
      const double X_n,
      const double X_p,
      double *restrict normalized_X_n,
      double *restrict normalized_X_p) {
  const double gamma_64 = 64.0 * DBL_EPSILON / (1.0 - 64.0 * DBL_EPSILON);
  const double fraction_roundoff = 27.0 * gamma_64;
  if(X_n < -fraction_roundoff || X_n > 1.0 + fraction_roundoff
     || X_p < -fraction_roundoff
     || X_p > 1.0 + fraction_roundoff) {
    return ghl_error_nrpyleakage_blocking;
  }
  *normalized_X_n = fmin(1.0, fmax(0.0, X_n));
  *normalized_X_p = fmin(1.0, fmax(0.0, X_p));
  return ghl_success;
}

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

/* Strict raw calculation with every channel enabled. */
ghl_error_codes_t ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double neutrino_degeneracy[ghl_m1_nrpyleakage_species_count],
      ghl_m1_nrpyleakage_raw_rates *restrict raw);

/* Provider-only variant. Pair, plasmon, bremsstrahlung, and scattering raw
 * calculations are performed only for enabled channels; equilibrium moments
 * and beta/Kirchhoff diagnostics remain available independently of the mask.
 * Disabled thermal-channel and scattering fields are zero. */
ghl_error_codes_t ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_with_mask(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double neutrino_degeneracy[ghl_m1_nrpyleakage_species_count],
      int channel_mask,
      ghl_m1_nrpyleakage_raw_rates *restrict raw);

#endif
