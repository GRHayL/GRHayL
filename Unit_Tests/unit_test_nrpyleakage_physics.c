#include <float.h>

// clang-format off: the private header requires GRHayL's public type setup.
#include "ghl_unit_tests.h"
#include "../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h"
// clang-format on

/**
 * @file unit_test_nrpyleakage_physics.c
 * @brief Table-free physical checks for NRPyLeakage charged-current helpers.
 *
 * These checks use high-precision blocking references and exact equilibrium
 * identities. They complement EOS-table golden replays, which only detect
 * changes relative to one earlier implementation.
 */

/**
 * Report a failed floating-point comparison.
 *
 * @param[in] quantity Name printed when the comparison fails.
 * @param[in] computed Value produced by the implementation.
 * @param[in] expected Independently calculated reference value.
 * @param[in] tolerance Maximum permitted absolute difference.
 */
static void check_close(
      const char *restrict quantity,
      const double computed,
      const double expected,
      const double tolerance) {
  if(!isfinite(computed) || fabs(computed - expected) > tolerance) {
    ghl_error(
          "%s mismatch: expected %.17e, got %.17e, tolerance %.17e\n", quantity,
          expected, computed, tolerance);
  }
}

/**
 * Check one high-precision nucleon-blocking reference state.
 *
 * @param[in] rho_cgs Rest-mass density in g/cm^3.
 * @param[in] T Temperature in MeV.
 * @param[in] X_n Free-neutron mass fraction.
 * @param[in] X_p Free-proton mass fraction.
 * @param[in] expected_B_n Expected neutron vacancy fraction.
 * @param[in] expected_B_p Expected proton vacancy fraction.
 * @param[in] expected_Y_np Expected neutron-to-proton transition population.
 * @param[in] expected_Y_pn Expected proton-to-neutron transition population.
 */
static void check_blocking_state(
      const double rho_cgs,
      const double T,
      const double X_n,
      const double X_p,
      const double expected_B_n,
      const double expected_B_p,
      const double expected_Y_np,
      const double expected_Y_pn) {
  double B_n, B_p, Y_np, Y_pn, eta_n_minus_eta_p;
  const ghl_error_codes_t error = NRPyLeakage_compute_nucleon_blocking(
        rho_cgs, T, X_n, X_p, &B_n, &B_p, &Y_np, &Y_pn, &eta_n_minus_eta_p);
  ghl_abort_if_error(error);

  /*
   * Fukushima quotes at most seven units in the last place for each inverse
   * fit. The 128-epsilon allowance covers two inversions, the forward fit,
   * libm operations, and ordinary platform rounding. It is a numerical error
   * budget, not a fitted physical tolerance.
   */
  const double reference_factor = 128.0 * DBL_EPSILON;
  check_close(
        "B_n", B_n, expected_B_n, reference_factor * fmax(1.0, fabs(expected_B_n)));
  check_close(
        "B_p", B_p, expected_B_p, reference_factor * fmax(1.0, fabs(expected_B_p)));
  check_close(
        "Y_np", Y_np, expected_Y_np, reference_factor * fmax(1.0, fabs(expected_Y_np)));
  check_close(
        "Y_pn", Y_pn, expected_Y_pn, reference_factor * fmax(1.0, fabs(expected_Y_pn)));

  if(B_n < 0.0 || B_n > X_n || B_p < 0.0 || B_p > X_p || Y_np < 0.0 || Y_np > X_n
     || Y_pn < 0.0 || Y_pn > X_p) {
    ghl_error("Nucleon blocking violates population bounds\n");
  }

  const double identity_tolerance = 16.0 * DBL_EPSILON * fmax(X_n, X_p);
  check_close(
        "transition-population identity", Y_np - Y_pn, X_n - X_p, identity_tolerance);
}

/** Check threshold orientation and its zero-shift limit. */
static void check_shifted_moments(void) {
  const double eta = -2.0;
  nrpyl_beta_moments zero_shift, neutrino_threshold, lepton_shift;
  ghl_abort_if_error(nrpyl_compute_shifted_fermi_moments(0.0, 0.0, eta, &zero_shift));
  ghl_abort_if_error(
        nrpyl_compute_shifted_fermi_moments(2.0, 0.0, eta, &neutrino_threshold));
  ghl_abort_if_error(nrpyl_compute_shifted_fermi_moments(0.0, 2.0, eta, &lepton_shift));

  double F4, F5;
  ghl_abort_if_error(NRPyLeakage_Fermi_Dirac_integrals(4, eta, &F4));
  ghl_abort_if_error(NRPyLeakage_Fermi_Dirac_integrals(5, eta, &F5));
  check_close(
        "zero-shift number moment", zero_shift.number, F4,
        8.0 * DBL_EPSILON * fmax(1.0, fabs(F4)));
  check_close(
        "zero-shift energy moment", zero_shift.energy, F5,
        8.0 * DBL_EPSILON * fmax(1.0, fabs(F5)));

  if(!(neutrino_threshold.number > 0.0) || !(neutrino_threshold.energy > 0.0)
     || !(lepton_shift.number > 0.0) || !(lepton_shift.energy > 0.0)
     || neutrino_threshold.number != lepton_shift.number
     || !(neutrino_threshold.energy > lepton_shift.energy)) {
    ghl_error("Shifted beta moments violate threshold identities\n");
  }
}

/** Check spectral detailed balance for both charged-current orientations. */
static void check_spectral_detailed_balance(void) {
  const double rho_cgs = 1.0e14;
  const double T = 3.0;
  const double X_n = 0.8;
  const double X_p = 0.1;
  const double muhat = 5.0;
  const double mu_e = 2.0;
  double B_n, B_p, Y_np, Y_pn, eta_n_minus_eta_p;
  ghl_abort_if_error(NRPyLeakage_compute_nucleon_blocking(
        rho_cgs, T, X_n, X_p, &B_n, &B_p, &Y_np, &Y_pn, &eta_n_minus_eta_p));

  const double q = nrpyl_compute_reaction_shift(T, muhat, eta_n_minus_eta_p);
  const double E_nu = fabs(q) + 4.0 * T;
  const double nue_lepton_exponent = (mu_e - (E_nu + q)) / T;
  const double anue_lepton_exponent = (-mu_e - (E_nu - q)) / T;
  const double nue_ratio = (Y_pn / Y_np) * exp(nue_lepton_exponent);
  const double anue_ratio = (Y_np / Y_pn) * exp(anue_lepton_exponent);
  const double expected_nue_ratio = exp((mu_e - muhat - E_nu) / T);
  const double expected_anue_ratio = exp((muhat - mu_e - E_nu) / T);

  const double balance_factor = 64.0 * DBL_EPSILON;
  check_close(
        "nue spectral detailed balance", nue_ratio, expected_nue_ratio,
        balance_factor * fmax(fabs(nue_ratio), fabs(expected_nue_ratio)));
  check_close(
        "anue spectral detailed balance", anue_ratio, expected_anue_ratio,
        balance_factor * fmax(fabs(anue_ratio), fabs(expected_anue_ratio)));
}

/**
 * Check the production beta-moment helpers against independent references.
 *
 * The references are 100-digit evaluations of the published Takahashi
 * complete-Fermi fits followed by the documented shifted-moment and vacancy
 * algebra.  The 128-epsilon allowance covers the fit evaluation, the short
 * arithmetic chain, and ordinary libm/platform rounding; it is not a fitted
 * physical tolerance.
 */
static void check_beta_moment_references(void) {
  const double T = 3.0;
  const double mu_e = 2.0;
  const double q = 1.2;
  const double tolerance_factor = 128.0 * DBL_EPSILON;
  nrpyl_beta_moments moments;

  ghl_abort_if_error(nrpyl_compute_beta_emission_moments(T, mu_e, 0.3, 1, q, &moments));
  check_close(
        "nue emission number moment", moments.number, 36.084128434093910748,
        tolerance_factor * 36.084128434093910748);
  check_close(
        "nue emission energy moment", moments.energy, 177.37422178631469128,
        tolerance_factor * 177.37422178631469128);

  ghl_abort_if_error(
        nrpyl_compute_beta_absorption_moments(T, mu_e, 0.3, 1, q, &moments));
  check_close(
        "nue absorption number moment", moments.number, 15.775326802419575840,
        tolerance_factor * 15.775326802419575840);
  check_close(
        "nue absorption energy moment", moments.energy, 24.290732422591624328,
        tolerance_factor * 24.290732422591624328);

  ghl_abort_if_error(
        nrpyl_compute_beta_emission_moments(T, mu_e, -0.2, -1, q, &moments));
  check_close(
        "anue emission number moment", moments.number, 14.618150280375924694,
        tolerance_factor * 14.618150280375924694);
  check_close(
        "anue emission energy moment", moments.energy, 76.860999352712340697,
        tolerance_factor * 76.860999352712340697);

  ghl_abort_if_error(
        nrpyl_compute_beta_absorption_moments(T, mu_e, -0.2, -1, q, &moments));
  check_close(
        "anue absorption number moment", moments.number, 10.405884737783058106,
        tolerance_factor * 10.405884737783058106);
  check_close(
        "anue absorption energy moment", moments.energy, 17.498996065988470906,
        tolerance_factor * 17.498996065988470906);
}

/**
 * Check the two floating-point limits of strongly blocked channels.
 *
 * Cold EOS states can retain a positive trace population while the complete
 * shifted Fermi moments all lie below double precision. Emission then reaches
 * its representable zero-rate limit, while absorption can retain a finite
 * normalized Boltzmann ratio. A second state checks a finite ratio of paired
 * subnormal moments.
 */
static void check_strongly_blocked_channel_limit(void) {
  // Cold trace-proton state exposed by the existing luminosity replay.
  const double rho_cgs = 8.65629955239776193e-9 * NRPyLeakage_units_geom_to_cgs_D;
  const double T = 1.13314052232918786e-3;
  const double muhat = 2.49813562281790986e1;
  const double mu_e = 5.85176444146208308;
  const double X_n = 1.68105605964749355e-1;
  const double X_p = 3.62177668230181617e-99;
  double B_n, B_p, Y_np, Y_pn, eta_n_minus_eta_p;
  ghl_abort_if_error(NRPyLeakage_compute_nucleon_blocking(
        rho_cgs, T, X_n, X_p, &B_n, &B_p, &Y_np, &Y_pn, &eta_n_minus_eta_p));

  const double q = nrpyl_compute_reaction_shift(T, muhat, eta_n_minus_eta_p);
  const double eta_nue = (mu_e - muhat) / T;
  const double eta_anue = -eta_nue;
  nrpyl_beta_moments moments[4];
  ghl_abort_if_error(
        nrpyl_compute_beta_emission_moments(T, mu_e, eta_nue, 1, q, &moments[0]));
  ghl_abort_if_error(
        nrpyl_compute_beta_emission_moments(T, mu_e, eta_anue, -1, q, &moments[1]));
  ghl_abort_if_error(
        nrpyl_compute_beta_absorption_moments(T, mu_e, eta_nue, 1, q, &moments[2]));
  ghl_abort_if_error(
        nrpyl_compute_beta_absorption_moments(T, mu_e, eta_anue, -1, q, &moments[3]));

  for(int i = 0; i < 4; i++) {
    if(!isfinite(moments[i].number) || !isfinite(moments[i].energy)
       || moments[i].number < 0.0 || moments[i].energy < 0.0) {
      ghl_error("Strongly blocked channel produced invalid moments\n");
    }
  }
  if(moments[0].number != 0.0 || moments[0].energy != 0.0 || moments[1].number != 0.0
     || moments[1].energy != 0.0 || !(moments[2].number > 0.0)) {
    ghl_error("Strongly blocked channel has the wrong limiting behavior\n");
  }

  /* Exact Boltzmann-limit oracle after numerator and denominator underflow. */
  nrpyl_beta_moments boltzmann_limit;
  ghl_abort_if_error(nrpyl_compute_beta_absorption_moments(
        1.0, 0.0, -750.0, 1, 2.0, &boltzmann_limit));
  const double limit_tolerance = 128.0 * DBL_EPSILON;
  check_close(
        "underflow-limit absorption number", boltzmann_limit.number, 27.947940644986954,
        limit_tolerance * 27.947940644986954);
  check_close(
        "underflow-limit absorption energy", boltzmann_limit.energy, 39.925629492838503,
        limit_tolerance * 39.925629492838503);

  /*
   * A large threshold can underflow only the shifted numerator. These
   * references retain the independently evaluated, ordinary F2(10) and
   * F3(10) denominators. The absolute allowance is 16 subnormal ulps because
   * a relative tolerance itself rounds to zero at this scale.
   */
  nrpyl_beta_moments mixed_underflow;
  ghl_abort_if_error(nrpyl_compute_beta_absorption_moments(
        1.0, -1000.0, 10.0, 1, -756.0, &mixed_underflow));
  const double subnormal_tolerance = 16.0 * nextafter(0.0, 1.0);
  check_close(
        "mixed-underflow absorption number", mixed_underflow.number,
        3.266449103024899e-321, subnormal_tolerance);
  check_close(
        "mixed-underflow absorption energy", mixed_underflow.energy,
        3.021736809070256e-319, subnormal_tolerance);

  // Cold trace-neutron state that exercises a finite subnormal-moment ratio.
  const double trace_neutron_rho_cgs
        = 2.50071396560124250e-7 * NRPyLeakage_units_geom_to_cgs_D;
  const double trace_neutron_T = 4.83765808590545313e-2;
  const double trace_neutron_muhat = -1.33640031126570094e1;
  const double trace_neutron_mu_e = 2.23964956356413083e1;
  ghl_abort_if_error(NRPyLeakage_compute_nucleon_blocking(
        trace_neutron_rho_cgs, trace_neutron_T, 4.84260997836237311e-97,
        1.29396803282191461e-10, &B_n, &B_p, &Y_np, &Y_pn, &eta_n_minus_eta_p));
  const double trace_neutron_q = nrpyl_compute_reaction_shift(
        trace_neutron_T, trace_neutron_muhat, eta_n_minus_eta_p);
  const double trace_neutron_eta_anue
        = (trace_neutron_muhat - trace_neutron_mu_e) / trace_neutron_T;
  nrpyl_beta_moments subnormal_ratio;
  ghl_abort_if_error(nrpyl_compute_beta_absorption_moments(
        trace_neutron_T, trace_neutron_mu_e, trace_neutron_eta_anue, -1, trace_neutron_q,
        &subnormal_ratio));
  if(!(subnormal_ratio.number > 0.0) || !(subnormal_ratio.energy > 0.0)
     || !isfinite(subnormal_ratio.number) || !isfinite(subnormal_ratio.energy)) {
    ghl_error("Paired subnormal beta moments did not produce a finite ratio\n");
  }
}

/**
 * Run all table-free NRPyLeakage physics checks.
 *
 * @return Zero after every check passes.
 */
int main(void) {
  check_blocking_state(
        1.0e14, 10.0, 0.5, 0.5, 0.260084127753509933, 0.260084127753509933,
        0.299113719948796048, 0.299113719948796048);
  check_blocking_state(
        1.0e14, 10.0, 0.5, 0.50000000000001, 0.260084127753509933, 0.260084127753512119,
        0.299113719948792897, 0.299113719948802897);
  check_blocking_state(
        1.0e15, 0.1, 0.5, 0.5, 0.000839824910039479330, 0.000839824910039479330,
        0.000841236154182995395, 0.000841236154182995395);
  check_blocking_state(
        1.0e13, 1.0, 0.9, 0.01, 0.180423946447079889, 0.01, 0.890297469419761934,
        0.000297469419761934053);
  check_shifted_moments();
  check_spectral_detailed_balance();
  check_beta_moment_references();
  check_strongly_blocked_channel_limit();
  ghl_info("NRPyLeakage physics checks passed\n");
  return 0;
}
