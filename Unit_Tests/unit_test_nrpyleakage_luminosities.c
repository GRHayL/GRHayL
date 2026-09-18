#include <float.h>

#include "ghl_unit_tests.h"

static inline bool luminosity_pert_test_fail(
      const double trusted,
      const double computed,
      const double perturbed) {
  if(trusted == 0.0 && perturbed == 0.0) {
    return computed != 0.0;
  }
  const double rel_tol = 1024.0 * DBL_EPSILON;
  const double abs_tol = 0.0;
  return ghl_pert_test_fail_with_tolerance(
        trusted, computed, perturbed, rel_tol, abs_tol);
}

static inline void check_fermi_dirac_integral(
    const int k,
    const double z,
    const double expected) {

  double integral = 0.0;
  const ghl_error_codes_t error =
    NRPyLeakage_Fermi_Dirac_integrals(k, z, &integral);
  if(error != ghl_success) {
    ghl_error("NRPyLeakage_Fermi_Dirac_integrals(k=%d, z=%.17e) returned error code %d\n",
              k, z, error);
  }

  const double atol = 1e-15;
  const double rtol = 1e-14;
  if(!isfinite(expected)
     || !(fabs(integral - expected) <= atol + rtol * fabs(expected))) {
    ghl_error("Incorrect Fermi-Dirac integral for k=%d, z=%.17e: expected %.17e, got %.17e\n",
              k, z, expected, integral);
  }
}

static inline void test_fermi_dirac_integrals(void) {
  const double z_large = 1.0e-2;
  const double z_small = 1.0e-4;

  check_fermi_dirac_integral(0, z_large, log(exp(z_large) + 1.0));
  check_fermi_dirac_integral(1, z_large,
                             (0.5*z_large*z_large + 1.6449)/(1.0 + exp(-1.6855*z_large)));
  check_fermi_dirac_integral(2, z_large,
                             (((1.0/3.0)*z_large*z_large*z_large) + 3.2898999999999998*z_large)
                             /(1.0 - exp(-1.8246*z_large)));

  check_fermi_dirac_integral(0, z_small, log(exp(z_small) + 1.0));
  check_fermi_dirac_integral(1, z_small,
                             exp(z_small)/(0.21590000000000001*exp(0.88570000000000004*z_small) + 1.0));
}

static double bounded_perturbation(
      const double base,
      const double delta,
      const double minimum,
      const double maximum) {
  double perturbed = base * (1.0 + delta);
  if(perturbed < minimum || perturbed > maximum) {
    perturbed = base * (1.0 - delta);
  }
  if(!robust_isfinite(perturbed) || perturbed < minimum || perturbed > maximum) {
    ghl_error("Could not perturb bounded luminosity input\n");
  }
  return perturbed;
}

static void check_input_pair(
      const char *restrict name,
      const double base,
      const double perturbed,
      const double minimum,
      const double maximum) {
  const double relative_tolerance = 1.0e-14 + 8.0 * DBL_EPSILON;
  const int invalid_zero_pair = base == 0.0 && perturbed != 0.0;
  const int invalid_relative_pair
        = base != 0.0 && fabs(perturbed / base - 1.0) > relative_tolerance;
  if(!robust_isfinite(base) || !robust_isfinite(perturbed) || perturbed < minimum
     || perturbed > maximum || invalid_zero_pair || invalid_relative_pair) {
    ghl_error(
          "Invalid luminosity input pair %s: base %.17e, perturbed %.17e\n", name, base,
          perturbed);
  }
}

void
generate_test_data(const ghl_eos_parameters *restrict eos) {

  srand(100);

  const int npoints = 1024;

  FILE *fp[2];
  fp[0] = fopen_with_check("nrpyleakage_luminosities_unperturbed.bin", "wb");
  fp[1] = fopen_with_check("nrpyleakage_luminosities_perturbed.bin", "wb");
  fwrite(&npoints, sizeof(int), 1, fp[0]);
  fwrite(&npoints, sizeof(int), 1, fp[1]);

  for(int n = 0; n < npoints; n++) {

    // Get random metric values
    double base_alpha;
    __attribute__((unused)) double betax, betay, betaz;
    double base_gammaxx, base_gammaxy, base_gammaxz;
    double base_gammayy, base_gammayz, base_gammazz;
    ghl_randomize_metric(
          &base_alpha, &betax, &betay, &betaz, &base_gammaxx, &base_gammaxy,
          &base_gammaxz, &base_gammayy, &base_gammayz, &base_gammazz);

    // Get random primitive values
    const double base_rho = fmin(
          eos->rho_max,
          fmax(eos->rho_min, pow(10, randf(log10(eos->rho_min), log10(eos->rho_max)))));
    const double base_Y_e = randf(eos->Y_e_min, eos->Y_e_max);
    const double base_T = fmin(
          eos->T_max,
          fmax(eos->T_min, pow(10, randf(log10(eos->T_min), log10(eos->T_max)))));
    const double base_W = randf(1, 10);

    // Get random optical depths (not sure these are reasonable values)
    ghl_neutrino_optical_depths base_tau;
    base_tau.nue[0] = randf(1, 1000);
    base_tau.nue[1] = randf(1, 1000);
    base_tau.anue[0] = randf(1, 1000);
    base_tau.anue[1] = randf(1, 1000);
    base_tau.nux[0] = randf(1, 1000);
    base_tau.nux[1] = randf(1, 1000);

    // Evaluate this base state, then its small perturbation, so that the two
    // output rows always form a matched pair.
    for(int perturb = 0; perturb <= 1; perturb++) {
      double alpha = base_alpha;
      double gammaxx = base_gammaxx, gammaxy = base_gammaxy;
      double gammaxz = base_gammaxz, gammayy = base_gammayy;
      double gammayz = base_gammayz, gammazz = base_gammazz;
      double rho = base_rho, Y_e = base_Y_e, T = base_T, W = base_W;
      ghl_neutrino_optical_depths tau = base_tau;

      if( perturb ) {
        alpha       *= (1+randf(-1,1)*1e-14);
        gammaxx     *= (1+randf(-1,1)*1e-14);
        gammaxy     *= (1+randf(-1,1)*1e-14);
        gammaxz     *= (1+randf(-1,1)*1e-14);
        gammayy     *= (1+randf(-1,1)*1e-14);
        gammayz     *= (1+randf(-1,1)*1e-14);
        gammazz     *= (1+randf(-1,1)*1e-14);
        rho = bounded_perturbation(
              base_rho, randf(-1, 1) * 1e-14, eos->rho_min, eos->rho_max);
        Y_e = bounded_perturbation(
              base_Y_e, randf(-1, 1) * 1e-14, eos->Y_e_min, eos->Y_e_max);
        T = bounded_perturbation(base_T, randf(-1, 1) * 1e-14, eos->T_min, eos->T_max);
        W = bounded_perturbation(base_W, randf(-1, 1) * 1e-14, 1.0, 10.0);
        tau.nue[0]
              = bounded_perturbation(base_tau.nue[0], randf(-1, 1) * 1e-14, 1.0, 1000.0);
        tau.nue[1]
              = bounded_perturbation(base_tau.nue[1], randf(-1, 1) * 1e-14, 1.0, 1000.0);
        tau.anue[0] = bounded_perturbation(
              base_tau.anue[0], randf(-1, 1) * 1e-14, 1.0, 1000.0);
        tau.anue[1] = bounded_perturbation(
              base_tau.anue[1], randf(-1, 1) * 1e-14, 1.0, 1000.0);
        tau.nux[0]
              = bounded_perturbation(base_tau.nux[0], randf(-1, 1) * 1e-14, 1.0, 1000.0);
        tau.nux[1]
              = bounded_perturbation(base_tau.nux[1], randf(-1, 1) * 1e-14, 1.0, 1000.0);
      }

      check_input_pair("alpha", base_alpha, alpha, -DBL_MAX, DBL_MAX);
      check_input_pair("gammaxx", base_gammaxx, gammaxx, -DBL_MAX, DBL_MAX);
      check_input_pair("gammaxy", base_gammaxy, gammaxy, -DBL_MAX, DBL_MAX);
      check_input_pair("gammaxz", base_gammaxz, gammaxz, -DBL_MAX, DBL_MAX);
      check_input_pair("gammayy", base_gammayy, gammayy, -DBL_MAX, DBL_MAX);
      check_input_pair("gammayz", base_gammayz, gammayz, -DBL_MAX, DBL_MAX);
      check_input_pair("gammazz", base_gammazz, gammazz, -DBL_MAX, DBL_MAX);
      check_input_pair("rho", base_rho, rho, eos->rho_min, eos->rho_max);
      check_input_pair("Y_e", base_Y_e, Y_e, eos->Y_e_min, eos->Y_e_max);
      check_input_pair("T", base_T, T, eos->T_min, eos->T_max);
      check_input_pair("W", base_W, W, 1.0, 10.0);
      check_input_pair("tau.nue[0]", base_tau.nue[0], tau.nue[0], 1.0, 1000.0);
      check_input_pair("tau.nue[1]", base_tau.nue[1], tau.nue[1], 1.0, 1000.0);
      check_input_pair("tau.anue[0]", base_tau.anue[0], tau.anue[0], 1.0, 1000.0);
      check_input_pair("tau.anue[1]", base_tau.anue[1], tau.anue[1], 1.0, 1000.0);
      check_input_pair("tau.nux[0]", base_tau.nux[0], tau.nux[0], 1.0, 1000.0);
      check_input_pair("tau.nux[1]", base_tau.nux[1], tau.nux[1], 1.0, 1000.0);

      // Compute luminosities
      ghl_neutrino_luminosities lum;
      ghl_error_codes_t error = NRPyLeakage_compute_neutrino_luminosities(eos, alpha,
                                                                          gammaxx, gammaxy, gammaxz,
                                                                          gammayy, gammayz, gammazz,
                                                                          rho, Y_e, T, W,
                                                                          &tau, &lum);
      ghl_abort_if_error(error);

      // Output to file
      if( !perturb ) {
        fwrite(&alpha, sizeof(double), 1, fp[perturb]);
        fwrite(&gammaxx, sizeof(double), 1, fp[perturb]);
        fwrite(&gammaxy, sizeof(double), 1, fp[perturb]);
        fwrite(&gammaxz, sizeof(double), 1, fp[perturb]);
        fwrite(&gammayy, sizeof(double), 1, fp[perturb]);
        fwrite(&gammayz, sizeof(double), 1, fp[perturb]);
        fwrite(&gammazz, sizeof(double), 1, fp[perturb]);
        fwrite(&rho, sizeof(double), 1, fp[perturb]);
        fwrite(&Y_e, sizeof(double), 1, fp[perturb]);
        fwrite(&T, sizeof(double), 1, fp[perturb]);
        fwrite(&W, sizeof(double), 1, fp[perturb]);
        fwrite(&tau, sizeof(ghl_neutrino_optical_depths), 1, fp[perturb]);
      }
      fwrite(&lum, sizeof(ghl_neutrino_luminosities), 1, fp[perturb]);
    }
  }
  fclose(fp[0]);
  fclose(fp[1]);
}

void
run_unit_test(const ghl_eos_parameters *restrict eos) {
  if(!luminosity_pert_test_fail(0.0, DBL_MIN, 0.0)
     || luminosity_pert_test_fail(0.0, 0.0, 0.0)) {
    ghl_error("Luminosity comparison did not enforce exact-zero references\n");
  }
  test_fermi_dirac_integrals();

  const double boundary_cases[][2] = {
    { 1.0, -1.0e-14 },
    { 1.0, 1.0e-14 },
    { 10.0, -1.0e-14 },
    { 10.0, 1.0e-14 },
  };
  for(size_t i = 0; i < sizeof(boundary_cases) / sizeof(boundary_cases[0]); i++) {
    const double perturbed
          = bounded_perturbation(boundary_cases[i][0], boundary_cases[i][1], 1.0, 10.0);
    check_input_pair("boundary endpoint", boundary_cases[i][0], perturbed, 1.0, 10.0);
  }

  int n1, n2;

  FILE *fp_unpert = fopen_with_check("nrpyleakage_luminosities_unperturbed.bin", "rb");
  FILE *fp_pert   = fopen_with_check("nrpyleakage_luminosities_perturbed.bin", "rb");

  int err = 0;
  err += fread(&n1, sizeof(int), 1, fp_unpert);
  err += fread(&n2, sizeof(int), 1, fp_pert  );
  if( err != 2 || n1 != n2 ) {
    fclose(fp_unpert);
    fclose(fp_pert);
    ghl_error("Problem reading number of points from file (err: %d, n1: %d, n2: %d)\n",
                 err, n1, n2);
  }

  const int npoints=n1;
  for(int n=0;n<npoints;n++) {

    // Read metric and primitive quantities from the unperturbed data file
    double alpha;
    double gammaxx, gammaxy, gammaxz, gammayy, gammayz, gammazz;
    double rho, Y_e, T, W;
    ghl_neutrino_optical_depths tau;

    err  = 0;
    err += fread(&alpha  , sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammaxx, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammaxy, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammaxz, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammayy, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammayz, sizeof(double)                 , 1, fp_unpert);
    err += fread(&gammazz, sizeof(double)                 , 1, fp_unpert);
    err += fread(&rho    , sizeof(double)                 , 1, fp_unpert);
    err += fread(&Y_e    , sizeof(double)                 , 1, fp_unpert);
    err += fread(&T      , sizeof(double)                 , 1, fp_unpert);
    err += fread(&W      , sizeof(double)                 , 1, fp_unpert);
    err += fread(&tau    , sizeof(ghl_neutrino_optical_depths), 1, fp_unpert);

    if( err != 12 ) {
      fclose(fp_unpert); fclose(fp_pert);
      ghl_error("Failed to read inputs from unperturbed data file\n");
    }

    // Compute luminosities
    ghl_neutrino_luminosities lum;
    ghl_error_codes_t error = NRPyLeakage_compute_neutrino_luminosities(eos, alpha,
                                                                        gammaxx, gammaxy, gammaxz,
                                                                        gammayy, gammayz, gammazz,
                                                                        rho, Y_e, T, W,
                                                                        &tau, &lum);
    ghl_abort_if_error(error);

    // Now read luminosities from unperturbed and perturbed data files
    ghl_neutrino_luminosities lum_trusted, lum_pert;
    if( 1 != fread(&lum_trusted, sizeof(ghl_neutrino_luminosities), 1, fp_unpert) ) {
      fclose(fp_unpert); fclose(fp_pert);
      ghl_error("Failed to read luminosities from unperturbed data file\n");
    }

    if( 1 != fread(&lum_pert, sizeof(ghl_neutrino_luminosities), 1, fp_pert) ) {
      fclose(fp_unpert); fclose(fp_pert);
      ghl_error("Failed to read luminosities from perturbed data file\n");
    }

    if(luminosity_pert_test_fail(lum_trusted.nue, lum.nue, lum_pert.nue)) {
      fclose(fp_unpert);
      fclose(fp_pert);
      ghl_error(
            "Validation failed for lum.nue at row %d: trusted %.17e, computed %.17e, "
            "perturbed %.17e\n",
            n, lum_trusted.nue, lum.nue, lum_pert.nue);
    }
    if(luminosity_pert_test_fail(lum_trusted.anue, lum.anue, lum_pert.anue)) {
      fclose(fp_unpert);
      fclose(fp_pert);
      ghl_error(
            "Validation failed for lum.anue at row %d: trusted %.17e, computed %.17e, "
            "perturbed %.17e\n",
            n, lum_trusted.anue, lum.anue, lum_pert.anue);
    }
    if(luminosity_pert_test_fail(lum_trusted.nux, lum.nux, lum_pert.nux)) {
      fclose(fp_unpert);
      fclose(fp_pert);
      ghl_error(
            "Validation failed for lum.nux at row %d: trusted %.17e, computed %.17e, "
            "perturbed %.17e\n",
            n, lum_trusted.nux, lum.nux, lum_pert.nux);
    }
  }
  fclose(fp_unpert);
  fclose(fp_pert);
}

#include "nrpyleakage_main.h"
