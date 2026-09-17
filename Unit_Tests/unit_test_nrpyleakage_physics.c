#include <float.h>

// clang-format off: the private header requires GRHayL's public type setup.
#include "ghl_unit_tests.h"
#include "../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h"
#include "../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h"
// clang-format on

/**
 * @file unit_test_nrpyleakage_physics.c
 * @brief Table-free physical checks for NRPyLeakage kernels and stencils.
 *
 * These checks use high-precision blocking references, exact equilibrium
 * identities, stable order-zero Fermi values, and a hand-calculated
 * asymmetric optical-depth stencil. They complement EOS-table golden replays,
 * which only detect changes relative to one earlier implementation.
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
  if(!robust_isfinite(computed) || fabs(computed - expected) > tolerance) {
    ghl_error(
          "%s mismatch: expected %.17e, got %.17e, tolerance %.17e\n", quantity,
          expected, computed, tolerance);
  }
}

/**
 * Construct a double from its object representation.
 *
 * @param[in] bits IEEE binary64 bit pattern.
 * @return Double with the requested representation.
 */
static double double_from_bits(const uint64_t bits) {
  double value;
  memcpy(&value, &bits, sizeof(value));
  return value;
}

/** Check finite and NaN classification from exact IEEE binary64 patterns. */
static void check_robust_classifiers(void) {
  static const struct {
    const char *name;
    uint64_t bits;
    int finite;
    int nan;
  } cases[] = {
    { "positive finite", UINT64_C(0x3ff0000000000000), 1, 0 },
    { "maximum finite", UINT64_C(0x7fefffffffffffff), 1, 0 },
    { "positive zero", UINT64_C(0x0000000000000000), 1, 0 },
    { "negative zero", UINT64_C(0x8000000000000000), 1, 0 },
    { "positive subnormal", UINT64_C(0x0000000000000001), 1, 0 },
    { "negative subnormal", UINT64_C(0x8000000000000001), 1, 0 },
    { "positive infinity", UINT64_C(0x7ff0000000000000), 0, 0 },
    { "negative infinity", UINT64_C(0xfff0000000000000), 0, 0 },
    { "quiet NaN", UINT64_C(0x7ff8000000000001), 0, 1 },
    { "signaling NaN", UINT64_C(0x7ff0000000000001), 0, 1 },
  };

  for(size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    const double value = double_from_bits(cases[i].bits);
    if(robust_isfinite(value) != cases[i].finite || robust_isnan(value) != cases[i].nan) {
      ghl_error("Robust classifier failed for %s\n", cases[i].name);
    }
  }
}

/**
 * Check analytic identities behind corrected production emission rates.
 *
 * These manufactured values independently expose the diffusion coefficient,
 * bremsstrahlung density power, and four-species heavy-lepton convention.
 */
static void check_emission_rate_identities(void) {
  check_close(
        "zero free-rate endpoint",
        nrpyl_effective_emission_rate(0.0, 0.0, 0.0, 0.0), 0.0, 0.0);
  check_close(
        "transparent endpoint",
        nrpyl_effective_emission_rate(2.0, 0.0, 0.0, 0.0), 2.0, 0.0);
  check_close(
        "zero-opacity endpoint",
        nrpyl_effective_emission_rate(2.0, 1.0, 0.0, 1.0), 0.0, 0.0);

  const double free_rate = NRPyLeakage_c_light / 6.0;
  const double effective_rate
        = nrpyl_effective_emission_rate(free_rate, 2.0, 2.0, 1.0);
  check_close(
        "diffusion suppression factor", effective_rate, free_rate / 3.0,
        16.0 * DBL_EPSILON * free_rate);

  const double Q_free_nux = 11.0;
  const double Q_eff_nux
        = nrpyl_effective_emission_rate(Q_free_nux, 3.0, 5.0, 7.0);
  const double expected_Q_eff_nux
        = Q_free_nux
          / (1.0 + 9.0 * (6.0 / NRPyLeakage_c_light) * Q_free_nux / 35.0);
  check_close(
        "heavy-lepton self-rate suppression", Q_eff_nux, expected_Q_eff_nux,
        16.0 * DBL_EPSILON * Q_free_nux);

  const double brems_rate
        = nrpyl_bremsstrahlung_number_rate(1.0, 2.0, 1.0, 0.0);
  const double doubled_density_rate
        = nrpyl_bremsstrahlung_number_rate(1.0, 4.0, 1.0, 0.0);
  const double expected_brems_rate
        = 4.0 * NRPyLeakage_Brems_C1 * NRPyLeakage_Brems_zeta;
  check_close(
        "bremsstrahlung analytic rate", brems_rate, expected_brems_rate,
        16.0 * DBL_EPSILON * expected_brems_rate);
  check_close(
        "bremsstrahlung density scaling", doubled_density_rate, 4.0 * brems_rate,
        16.0 * DBL_EPSILON * doubled_density_rate);
  const double brems_energy = nrpyl_bremsstrahlung_energy_rate(3.0, 7.0);
  const double expected_brems_energy
        = NRPyLeakage_Brems_C2 * 3.0 * 7.0 / NRPyLeakage_Brems_C1;
  check_close(
        "bremsstrahlung energy conversion", brems_energy,
        expected_brems_energy, 16.0 * DBL_EPSILON * expected_brems_energy);

  const double source = nrpyl_matter_energy_source(2.0, 3.0, 5.0);
  const double expected_source = -25.0 * NRPyLeakage_units_cgs_to_geom_Q;
  check_close(
        "four-species heavy-lepton source", source, expected_source,
        16.0 * DBL_EPSILON * fabs(expected_source));
}

/** Check the stable order-zero Fermi expression at numerical edge cases. */
static void check_fermi_dirac_zero_order(void) {
  static const double z[] = { 710.0, -40.0, 0.0, 1.0e-3 };
  static const double expected[] = {
    7.10000000000000000e2,
    4.24835425529158887e-18,
    6.93147180559945286e-1,
    6.93647305559940142e-1,
  };

  for(int i = 0; i < 4; i++) {
    double computed;
    ghl_abort_if_error(NRPyLeakage_Fermi_Dirac_integrals(0, z[i], &computed));
    /* Relative scaling remains meaningful for the tiny z=-40 result. */
    check_close(
          "order-zero Fermi integral", computed, expected[i],
          8.0 * DBL_EPSILON * fabs(expected[i]));
  }
}

/** Check every supported Fermi key on both sides of the fit branch. */
static void check_all_fermi_dirac_keys(void) {
  static const double z[] = { 1.0e-2, 1.0e-4 };
  static const double expected[2][6] = {
    { 6.98159680507862257e-1, 8.29406243971260060e-1,
      1.81959808582090532e0, 5.73653915904569534e0,
      2.35590999834500465e1, 1.19438368623386765e2 },
    { 6.93197181809945384e-1, 8.22505367332569515e-1,
      1.80326583833903631e0, 5.68289726190304290e0,
      2.33326899268987837e1, 1.18273220285188671e2 },
  };

  for(int branch = 0; branch < 2; branch++) {
    for(int key = 0; key < 6; key++) {
      double computed;
      ghl_abort_if_error(
            NRPyLeakage_Fermi_Dirac_integrals(key, z[branch], &computed));
      check_close(
            "Fermi integral branch/key", computed, expected[branch][key],
            16.0*DBL_EPSILON*fabs(expected[branch][key]));
    }
  }
}

/** Check every public leakage result's deterministic nonfinite fallback. */
static void check_output_fallbacks(void) {
  ghl_neutrino_opacities kappa = {
    .nue = { NAN, INFINITY },
    .anue = { -INFINITY, NAN },
    .nux = { INFINITY, -INFINITY },
  };
  ghl_neutrino_luminosities lum = {
    .nue = NAN,
    .anue = INFINITY,
    .nux = -INFINITY,
  };
  double R_source = NAN;
  double Q_source = -INFINITY;

  if(!nrpyl_sanitize_opacities(&kappa)
     || !nrpyl_sanitize_luminosities(&lum)
     || !nrpyl_sanitize_sources(&R_source, &Q_source)) {
    ghl_error("Nonfinite sanitizer did not report a replaced output\n");
  }

  const double opacity_outputs[] = {
    kappa.nue[0], kappa.nue[1], kappa.anue[0], kappa.anue[1],
    kappa.nux[0], kappa.nux[1],
  };
  for(size_t i = 0; i < sizeof(opacity_outputs)/sizeof(opacity_outputs[0]); i++) {
    if(opacity_outputs[i]
       != NRPyLeakage_units_geom_to_cgs_L * 1.0e-15)
      ghl_error("Opacity did not use the nonfinite floor\n");
  }
  const double rate_outputs[] = {
    lum.nue, lum.anue, lum.nux, R_source, Q_source,
  };
  for(size_t i = 0; i < sizeof(rate_outputs)/sizeof(rate_outputs[0]); i++) {
    if(rate_outputs[i] != 0.0)
      ghl_error("Emission or source output did not use the neutral fallback\n");
  }
}

#ifndef GHL_DISABLE_HDF5
static double mock_X_n = 0.7;
static double mock_X_p = 0.2;

static ghl_error_codes_t mock_eos_success(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p) {
  (void)eos;
  (void)rho;
  (void)Y_e;
  (void)T;
  *muhat = 5.0;
  *mu_e = 2.0;
  *mu_p = 1.0;
  *mu_n = 6.0;
  *X_n = mock_X_n;
  *X_p = mock_X_p;
  return ghl_success;
}

static ghl_error_codes_t mock_eos_shifted_reference(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p) {
  const ghl_error_codes_t error = mock_eos_success(
        eos, rho, Y_e, T, muhat, mu_e, mu_p, mu_n, X_n, X_p);
  *mu_p += 100.0;
  *mu_n += 100.0;
  return error;
}

static ghl_error_codes_t mock_eos_failure(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p) {
  (void)eos;
  (void)rho;
  (void)Y_e;
  (void)T;
  (void)muhat;
  (void)mu_e;
  (void)mu_p;
  (void)mu_n;
  (void)X_n;
  (void)X_p;
  return ghl_error_table_max_rho;
}
#endif

/** Check public leakage error, consistency, and finite-output contracts. */
static void check_public_api_contracts(void) {
  const ghl_eos_parameters eos = { 0 };
  const ghl_neutrino_optical_depths tau = { 0 };
  ghl_neutrino_opacities standalone = { 0 }, combined = { 0 };
  ghl_neutrino_luminosities lum = { 0 };
  double R_source = 0.0, Q_source = 0.0;

#ifdef GHL_DISABLE_HDF5
  if(NRPyLeakage_compute_neutrino_opacities(
           &eos, 1.0e-5, 0.1, 3.0, &tau, &standalone)
           != ghl_error_used_disabled_hdf5
     || NRPyLeakage_compute_neutrino_luminosities(
              &eos, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0e-5, 0.1,
              3.0, 1.0, &tau, &lum)
              != ghl_error_used_disabled_hdf5
     || NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
              &eos, 1.0e-5, 0.1, 3.0, &tau, &combined, &R_source, &Q_source)
              != ghl_error_used_disabled_hdf5) {
    ghl_error("Leakage API did not reject a disabled-HDF5 build\n");
  }
#else
  ghl_error_codes_t (*saved_eos_callback)(
        const ghl_eos_parameters *restrict, double, double, double,
        double *restrict, double *restrict, double *restrict, double *restrict,
        double *restrict, double *restrict)
        = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T;

  ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = mock_eos_failure;
  if(NRPyLeakage_compute_neutrino_opacities(
           &eos, 1.0e-5, 0.1, 3.0, &tau, &standalone)
           != ghl_error_table_max_rho
     || NRPyLeakage_compute_neutrino_luminosities(
              &eos, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0e-5, 0.1,
              3.0, 1.0, &tau, &lum)
              != ghl_error_table_max_rho
     || NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
              &eos, 1.0e-5, 0.1, 3.0, &tau, &combined, &R_source, &Q_source)
              != ghl_error_table_max_rho) {
    ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = saved_eos_callback;
    ghl_error("Leakage API did not propagate its EOS error\n");
  }

  ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = mock_eos_success;
  ghl_error_codes_t error = NRPyLeakage_compute_neutrino_opacities(
        &eos, 1.0e-5, 0.1, 3.0, &tau, &standalone);
  if(error != ghl_success)
    ghl_error("Standalone opacity mock returned error %d\n", error);
  error = NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
        &eos, 1.0e-5, 0.1, 3.0, &tau, &combined, &R_source, &Q_source);
  if(error != ghl_success)
    ghl_error("Combined leakage mock returned error %d\n", error);
  const double standalone_values[6] = {
    standalone.nue[0], standalone.nue[1], standalone.anue[0],
    standalone.anue[1], standalone.nux[0], standalone.nux[1],
  };
  const double combined_values[6] = {
    combined.nue[0], combined.nue[1], combined.anue[0],
    combined.anue[1], combined.nux[0], combined.nux[1],
  };
  for(int i = 0; i < 6; i++) {
    check_close(
          "combined opacity", combined_values[i], standalone_values[i],
          128.0*DBL_EPSILON*fmax(DBL_MIN, fabs(standalone_values[i])));
  }

  const ghl_neutrino_optical_depths thick_tau = {
    .nue = { 2.0, 3.0 },
    .anue = { 5.0, 7.0 },
    .nux = { 11.0, 13.0 },
  };
  error = NRPyLeakage_compute_neutrino_luminosities(
        &eos, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0e-5, 0.1,
        3.0, 1.0, &thick_tau, &lum);
  if(error != ghl_success)
    ghl_error("Finite-depth luminosity mock returned error %d\n", error);
  check_close(
        "public heavy-lepton luminosity suppression", lum.nux,
        7.42012378328398818e-15, 1.0e-12*7.42012378328398818e-15);
  error = NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
        &eos, 1.0e-5, 0.1, 3.0, &thick_tau, &combined,
        &R_source, &Q_source);
  if(error != ghl_success)
    ghl_error("Finite-depth combined leakage mock returned error %d\n", error);
  check_close(
        "public heavy-lepton source suppression", Q_source,
        -8.44271649383971011e-13, 1.0e-12*8.44271649383971011e-13);

  ghl_neutrino_luminosities free_lum = { 0 }, oracle_lum = { 0 };
  ghl_neutrino_opacities free_kappa = { 0 }, oracle_kappa = { 0 };
  double free_R = 0.0, free_Q = 0.0, oracle_R = 0.0, oracle_Q = 0.0;
  const ghl_neutrino_optical_depths oracle_tau = { .nux = { 0.0, 1000.0 } };
  error = NRPyLeakage_compute_neutrino_luminosities(
        &eos, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0e-5, 0.1,
        3.0, 1.0, &tau, &free_lum);
  error |= NRPyLeakage_compute_neutrino_luminosities(
        &eos, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0e-5, 0.1,
        3.0, 1.0, &oracle_tau, &oracle_lum);
  error |= NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
        &eos, 1.0e-5, 0.1, 3.0, &tau, &free_kappa, &free_R, &free_Q);
  error |= NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
        &eos, 1.0e-5, 0.1, 3.0, &oracle_tau, &oracle_kappa,
        &oracle_R, &oracle_Q);
  if(error != ghl_success)
    ghl_error("Public diffusion-oracle setup returned error %d\n", error);
  double F3_zero;
  ghl_abort_if_error(NRPyLeakage_Fermi_Dirac_integrals(3, 0.0, &F3_zero));
  const double Q_free_nux = free_lum.nux / NRPyLeakage_units_cgs_to_geom_Q;
  const double equilibrium_energy_density
        = 4.0*M_PI*pow(3.0, 4)*F3_zero/NRPyLeakage_hc3;
  const double kappa_cgs
        = free_kappa.nux[1]/NRPyLeakage_units_geom_to_cgs_L;
  const double Q_diff_nux = kappa_cgs*equilibrium_energy_density
                            * NRPyLeakage_c_light/(6.0*1000.0*1000.0);
  const double Q_eff_nux = Q_free_nux/(1.0 + Q_free_nux/Q_diff_nux);
  const double expected_lum_nux
        = NRPyLeakage_units_cgs_to_geom_Q*Q_eff_nux;
  const double expected_source_change
        = 4.0*NRPyLeakage_units_cgs_to_geom_Q*(Q_free_nux - Q_eff_nux);
  check_close(
        "public heavy-lepton diffusion oracle", oracle_lum.nux,
        expected_lum_nux, 1.0e-12*fabs(expected_lum_nux));
  check_close(
        "public heavy-lepton source diffusion oracle", oracle_Q - free_Q,
        expected_source_change, 1.0e-12*fabs(expected_source_change));

  const ghl_neutrino_optical_depths depths[2] = { tau, thick_tau };
  ghl_neutrino_opacities reference_standalone[2] = { 0 };
  ghl_neutrino_opacities reference_combined[2] = { 0 };
  ghl_neutrino_opacities shifted_standalone[2] = { 0 };
  ghl_neutrino_opacities shifted_combined[2] = { 0 };
  ghl_neutrino_luminosities reference_lum[2] = { 0 };
  ghl_neutrino_luminosities shifted_lum[2] = { 0 };
  double reference_R_source[2] = { 0.0 }, reference_Q_source[2] = { 0.0 };
  double shifted_R_source[2] = { 0.0 }, shifted_Q_source[2] = { 0.0 };
  bool reference_error = false, shifted_error = false;
  for(int d = 0; d < 2; d++) {
    reference_error |= NRPyLeakage_compute_neutrino_opacities(
          &eos, 1.0e-5, 0.1, 3.0, &depths[d], &reference_standalone[d])
          != ghl_success;
    reference_error |= NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
          &eos, 1.0e-5, 0.1, 3.0, &depths[d], &reference_combined[d],
          &reference_R_source[d], &reference_Q_source[d]) != ghl_success;
    reference_error |= NRPyLeakage_compute_neutrino_luminosities(
          &eos, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0e-5, 0.1,
          3.0, 1.0, &depths[d], &reference_lum[d]) != ghl_success;
  }
  ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T
        = mock_eos_shifted_reference;
  for(int d = 0; d < 2; d++) {
    shifted_error |= NRPyLeakage_compute_neutrino_opacities(
          &eos, 1.0e-5, 0.1, 3.0, &depths[d], &shifted_standalone[d])
          != ghl_success;
    shifted_error |= NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
          &eos, 1.0e-5, 0.1, 3.0, &depths[d], &shifted_combined[d],
          &shifted_R_source[d], &shifted_Q_source[d]) != ghl_success;
    shifted_error |= NRPyLeakage_compute_neutrino_luminosities(
          &eos, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0e-5, 0.1,
          3.0, 1.0, &depths[d], &shifted_lum[d]) != ghl_success;
  }
  ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = mock_eos_success;
  if(reference_error || shifted_error) {
    ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = saved_eos_callback;
    ghl_error("Common chemical-potential shift returned an error\n");
  }
  for(int d = 0; d < 2; d++) {
    const double reference_standalone_values[6] = {
      reference_standalone[d].nue[0], reference_standalone[d].nue[1],
      reference_standalone[d].anue[0], reference_standalone[d].anue[1],
      reference_standalone[d].nux[0], reference_standalone[d].nux[1],
    };
    const double shifted_standalone_values[6] = {
      shifted_standalone[d].nue[0], shifted_standalone[d].nue[1],
      shifted_standalone[d].anue[0], shifted_standalone[d].anue[1],
      shifted_standalone[d].nux[0], shifted_standalone[d].nux[1],
    };
    const double reference_combined_values[6] = {
      reference_combined[d].nue[0], reference_combined[d].nue[1],
      reference_combined[d].anue[0], reference_combined[d].anue[1],
      reference_combined[d].nux[0], reference_combined[d].nux[1],
    };
    const double shifted_combined_values[6] = {
      shifted_combined[d].nue[0], shifted_combined[d].nue[1],
      shifted_combined[d].anue[0], shifted_combined[d].anue[1],
      shifted_combined[d].nux[0], shifted_combined[d].nux[1],
    };
    for(int i = 0; i < 6; i++) {
      check_close(
            "standalone common-reference invariance",
            shifted_standalone_values[i], reference_standalone_values[i], 0.0);
      check_close(
            "combined common-reference invariance",
            shifted_combined_values[i], reference_combined_values[i], 0.0);
    }
    check_close(
          "number-source common-reference invariance",
          shifted_R_source[d], reference_R_source[d], 0.0);
    check_close(
          "energy-source common-reference invariance",
          shifted_Q_source[d], reference_Q_source[d], 0.0);
    check_close(
          "nue luminosity common-reference invariance",
          shifted_lum[d].nue, reference_lum[d].nue, 0.0);
    check_close(
          "anue luminosity common-reference invariance",
          shifted_lum[d].anue, reference_lum[d].anue, 0.0);
    check_close(
          "nux luminosity common-reference invariance",
          shifted_lum[d].nux, reference_lum[d].nux, 0.0);
  }

  error = NRPyLeakage_compute_neutrino_luminosities(
        &eos, NAN, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0e-5, 0.1,
        3.0, 1.0, &tau, &lum);
  if(error != ghl_error_nrpyleakage_nonfinite_output)
    ghl_error("Nonfinite luminosity mock returned error %d\n", error);
  ghl_neutrino_optical_depths nonfinite_tau = tau;
  nonfinite_tau.nue[1] = NAN;
  nonfinite_tau.anue[1] = NAN;
  nonfinite_tau.nux[1] = NAN;
  error = NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
        &eos, 1.0e-5, 0.1, 3.0, &nonfinite_tau, &combined, &R_source, &Q_source);
  if(error != ghl_error_nrpyleakage_nonfinite_output)
    ghl_error("Nonfinite combined leakage mock returned error %d\n", error);
  ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = saved_eos_callback;

  const double finite_outputs[] = {
    standalone.nue[0], standalone.nue[1], standalone.anue[0],
    standalone.anue[1], standalone.nux[0], standalone.nux[1],
    lum.nue, lum.anue, lum.nux, R_source, Q_source,
    combined.nue[0], combined.nue[1], combined.anue[0],
    combined.anue[1], combined.nux[0], combined.nux[1],
  };
  for(size_t i = 0; i < sizeof(finite_outputs)/sizeof(finite_outputs[0]); i++) {
    if(!robust_isfinite(finite_outputs[i]))
      ghl_error("Leakage API exposed a nonfinite final output\n");
  }
  if(lum.nue != 0.0 || lum.anue != 0.0 || lum.nux != 0.0
     || Q_source != 0.0) {
    ghl_error("Leakage API did not apply its nonfinite-output fallback\n");
  }
#endif
}

#ifndef GHL_DISABLE_HDF5
/** Evaluate all three public leakage paths with the current mock fractions. */
static void evaluate_mock_public_outputs(double outputs[17]) {
  const ghl_eos_parameters eos = { 0 };
  const ghl_neutrino_optical_depths tau = { 0 };
  ghl_neutrino_opacities standalone, combined;
  ghl_neutrino_luminosities lum;
  double R_source, Q_source;
  if(NRPyLeakage_compute_neutrino_opacities(
           &eos, 1.0e-5, 0.1, 3.0, &tau, &standalone) != ghl_success
     || NRPyLeakage_compute_neutrino_luminosities(
              &eos, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
              1.0e-5, 0.1, 3.0, 1.0, &tau, &lum) != ghl_success
     || NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(
              &eos, 1.0e-5, 0.1, 3.0, &tau, &combined, &R_source, &Q_source)
              != ghl_success) {
    ghl_error("Public leakage API rejected an endpoint-limit state\n");
  }

  const double values[17] = {
    standalone.nue[0], standalone.nue[1], standalone.anue[0],
    standalone.anue[1], standalone.nux[0], standalone.nux[1],
    lum.nue, lum.anue, lum.nux,
    combined.nue[0], combined.nue[1], combined.anue[0],
    combined.anue[1], combined.nux[0], combined.nux[1], R_source, Q_source,
  };
  for(size_t i = 0; i < 17; i++) {
    if(!robust_isfinite(values[i]))
      ghl_error("Public endpoint-limit output %zu is nonfinite\n", i);
    outputs[i] = values[i];
  }
}

/** Check the analytic public limit at single-species free-nucleon endpoints. */
static void check_public_single_species_endpoints(void) {
  ghl_error_codes_t (*saved_eos_callback)(
        const ghl_eos_parameters *restrict, double, double, double,
        double *restrict, double *restrict, double *restrict, double *restrict,
        double *restrict, double *restrict)
        = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T;
  ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = mock_eos_success;

  const double endpoint_pairs[][2] = {
    { 0.8, 0.0 },
    { 0.0, 0.2 },
  };
  const double minority_trace[] = { 1.0e-12, 1.0e-40, 1.0e-100, 1.0e-200, 1.0e-300 };
  for(size_t i = 0; i < sizeof(endpoint_pairs)/sizeof(endpoint_pairs[0]); i++) {
    double endpoint[17];
    mock_X_n = endpoint_pairs[i][0];
    mock_X_p = endpoint_pairs[i][1];
    evaluate_mock_public_outputs(endpoint);

    const double opacity_floor
          = NRPyLeakage_units_geom_to_cgs_L * 1.0e-15;
    for(size_t j = 0; j < 6; j++) {
      if(endpoint[j] <= opacity_floor)
        ghl_error("Standalone endpoint opacity %zu lost occupied-species scattering\n", j);
    }
    for(size_t j = 9; j < 15; j++) {
      if(endpoint[j] <= opacity_floor)
        ghl_error("Combined endpoint opacity %zu lost occupied-species scattering\n", j - 9);
    }
    if(!(endpoint[6] > 0.0) || !(endpoint[7] > 0.0) || !(endpoint[8] > 0.0)
       || endpoint[15] != 0.0 || endpoint[16] == 0.0)
      ghl_error("Endpoint lost non-beta emission or gained a beta source\n");

    double previous_error = INFINITY;
    double near_endpoint[17];
    for(size_t k = 0; k < sizeof(minority_trace)/sizeof(minority_trace[0]); k++) {
      mock_X_n = endpoint_pairs[i][0] == 0.0
                       ? minority_trace[k] : endpoint_pairs[i][0];
      mock_X_p = endpoint_pairs[i][1] == 0.0
                       ? minority_trace[k] : endpoint_pairs[i][1];
      evaluate_mock_public_outputs(near_endpoint);
      double trace_error = 0.0;
      for(size_t j = 0; j < 17; j++) {
        const double scaled_error
              = fabs(near_endpoint[j] - endpoint[j]) / fmax(1.0, fabs(endpoint[j]));
        trace_error = fmax(trace_error, scaled_error);
      }
      if(trace_error > previous_error + 512.0 * DBL_EPSILON)
        ghl_error("Single-species limit trace did not converge monotonically\n");
      previous_error = trace_error;
    }
    for(size_t j = 0; j < 17; j++) {
      check_close(
            "single-species analytic limit", endpoint[j], near_endpoint[j],
            512.0 * DBL_EPSILON * fmax(1.0, fabs(endpoint[j])));
    }
  }

  double normalized[17], exact[17];
  mock_X_n = -8.237492054256009e-17;
  mock_X_p = 0.2;
  evaluate_mock_public_outputs(normalized);
  mock_X_n = 0.0;
  evaluate_mock_public_outputs(exact);
  for(size_t i = 0; i < 17; i++) {
    check_close("normalized endpoint", normalized[i], exact[i], 0.0);
  }

  mock_X_n = 0.0;
  mock_X_p = 0.0;
  evaluate_mock_public_outputs(exact);

  mock_X_n = 0.7;
  mock_X_p = 0.2;
  ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T = saved_eos_callback;
}
#endif

/** Check roundoff normalization and exact endpoint overlap limits. */
static void check_mass_fraction_validation(void) {
  double B_n, B_p, Y_np, Y_pn, eta_n_minus_eta_p;
  ghl_error_codes_t error = NRPyLeakage_compute_nucleon_blocking(
        1.0e14, 1.0, -8.237492054256009e-17, 0.1, &B_n, &B_p, &Y_np, &Y_pn,
        &eta_n_minus_eta_p);
  if(error != ghl_success || B_n != 0.0 || !(B_p > 0.0) || Y_np != 0.0
     || Y_pn != 0.1 || eta_n_minus_eta_p != 0.0) {
    ghl_error("Roundoff-sized endpoint did not use the proton-only limit\n");
  }

  error = NRPyLeakage_compute_nucleon_blocking(
        1.0e14, 1.0, 0.8, 0.0, &B_n, &B_p, &Y_np, &Y_pn,
        &eta_n_minus_eta_p);
  if(error != ghl_success || !(B_n > 0.0) || B_p != 0.0 || Y_np != 0.8
     || Y_pn != 0.0 || eta_n_minus_eta_p != 0.0) {
    ghl_error("Neutron-only endpoint did not use its analytic overlap limit\n");
  }

  error = NRPyLeakage_compute_nucleon_blocking(
        1.0e14, 1.0, -1.0e-6, 0.1, &B_n, &B_p, &Y_np, &Y_pn, &eta_n_minus_eta_p);
  if(error != ghl_error_nrpyleakage_blocking) {
    ghl_error("Materially negative nucleon fraction returned error code %d\n", error);
  }

  error = NRPyLeakage_compute_nucleon_blocking(
        1.0e14, 1.0, 0.0, 0.0, &B_n, &B_p, &Y_np, &Y_pn,
        &eta_n_minus_eta_p);
  if(error != ghl_success || B_n != 0.0 || B_p != 0.0 || Y_np != 0.0
     || Y_pn != 0.0 || eta_n_minus_eta_p != 0.0) {
    ghl_error("Both-zero nucleon endpoint did not remain the neutral state\n");
  }
}

/** Check neighbor ordering with distinct data and unequal face metrics. */
static void check_asymmetric_optical_depth_stencil(void) {
  const double dxx[3] = { 1.0, 1.0, 1.0 };
  const double gxx[3] = { 1.0, 4.0, 16.0 };
  const double gyy[3] = { 4.0, 9.0, 25.0 };
  const double gzz[3] = { 9.0, 16.0, 36.0 };
  const ghl_neutrino_opacities kappa_im1
        = { .nue = { 4.0, 4.0 }, .anue = { 4.0, 4.0 }, .nux = { 4.0, 4.0 } };
  const ghl_neutrino_opacities kappa_ip1
        = { .nue = { 6.0, 6.0 }, .anue = { 6.0, 6.0 }, .nux = { 6.0, 6.0 } };
  const ghl_neutrino_opacities kappa_jm1
        = { .nue = { 8.0, 8.0 }, .anue = { 8.0, 8.0 }, .nux = { 8.0, 8.0 } };
  const ghl_neutrino_opacities kappa_jp1
        = { .nue = { 10.0, 10.0 }, .anue = { 10.0, 10.0 }, .nux = { 10.0, 10.0 } };
  const ghl_neutrino_opacities kappa_km1
        = { .nue = { 12.0, 12.0 }, .anue = { 12.0, 12.0 }, .nux = { 12.0, 12.0 } };
  const ghl_neutrino_opacities kappa_kp1
        = { .nue = { 14.0, 14.0 }, .anue = { 14.0, 14.0 }, .nux = { 14.0, 14.0 } };
  const ghl_neutrino_opacities kappa_center
        = { .nue = { 2.0, 2.0 }, .anue = { 2.0, 2.0 }, .nux = { 2.0, 2.0 } };
  /* Each field selects a different direction; 1000 excludes every other path. */
  const ghl_neutrino_optical_depths tau_im1 = { .nue = { 0.0, 1000.0 },
                                                .anue = { 1000.0, 1000.0 },
                                                .nux = { 1000.0, 1000.0 } };
  const ghl_neutrino_optical_depths tau_ip1 = { .nue = { 1000.0, 0.0 },
                                                .anue = { 1000.0, 1000.0 },
                                                .nux = { 1000.0, 1000.0 } };
  const ghl_neutrino_optical_depths tau_jm1 = { .nue = { 1000.0, 1000.0 },
                                                .anue = { 0.0, 1000.0 },
                                                .nux = { 1000.0, 1000.0 } };
  const ghl_neutrino_optical_depths tau_jp1 = { .nue = { 1000.0, 1000.0 },
                                                .anue = { 1000.0, 0.0 },
                                                .nux = { 1000.0, 1000.0 } };
  const ghl_neutrino_optical_depths tau_km1 = { .nue = { 1000.0, 1000.0 },
                                                .anue = { 1000.0, 1000.0 },
                                                .nux = { 0.0, 1000.0 } };
  const ghl_neutrino_optical_depths tau_kp1 = { .nue = { 1000.0, 1000.0 },
                                                .anue = { 1000.0, 1000.0 },
                                                .nux = { 1000.0, 0.0 } };
  ghl_neutrino_optical_depths computed;

  NRPyLeakage_optical_depths_PathOfLeastResistance(
        dxx, gxx, gyy, gzz, &kappa_im1, &kappa_ip1, &kappa_jm1, &kappa_jp1, &kappa_km1,
        &kappa_kp1, &tau_im1, &tau_ip1, &tau_jm1, &tau_jp1, &tau_km1, &tau_kp1,
        &kappa_center, &computed);

  const double expected[6] = {
    3.0 * sqrt(2.5),  4.0 * sqrt(10.0), 5.0 * sqrt(6.5),
    6.0 * sqrt(17.0), 7.0 * sqrt(12.5), 8.0 * sqrt(26.0),
  };
  const double computed_values[6]
        = { computed.nue[0],  computed.nue[1], computed.anue[0],
            computed.anue[1], computed.nux[0], computed.nux[1] };
  for(int i = 0; i < 6; i++) {
    check_close(
          "asymmetric optical depth", computed_values[i], expected[i],
          16.0 * DBL_EPSILON * expected[i]);
  }
}

/**
 * Check one high-precision nucleon-blocking reference state.
 *
 * @param[in] rho_cgs Rest-mass density in g/cm^3.
 * @param[in] T Temperature in MeV.
 * @param[in] X_n Free-neutron mass fraction.
 * @param[in] X_p Free-proton mass fraction.
 * @param[in] expected_B_n Expected effective neutron scattering population.
 * @param[in] expected_B_p Expected effective proton scattering population.
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
 * Check the production beta-moment helpers against high-precision arithmetic
 * references.
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
    if(!robust_isfinite(moments[i].number) || !robust_isfinite(moments[i].energy)
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
     || !robust_isfinite(subnormal_ratio.number)
     || !robust_isfinite(subnormal_ratio.energy)) {
    ghl_error("Paired subnormal beta moments did not produce a finite ratio\n");
  }
}

/**
 * Run all table-free NRPyLeakage physics checks.
 *
 * @return Zero after every check passes.
 */
int main(void) {
  check_robust_classifiers();
  check_emission_rate_identities();
  check_fermi_dirac_zero_order();
  check_all_fermi_dirac_keys();
  check_output_fallbacks();
  check_public_api_contracts();
#ifndef GHL_DISABLE_HDF5
  check_public_single_species_endpoints();
#endif
  check_mass_fraction_validation();
  check_asymmetric_optical_depth_stencil();
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
