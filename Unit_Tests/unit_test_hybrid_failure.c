#include "ghl_unit_tests.h"
#include "../GRHayL/Con2Prim/utils_Noble.h"

static void check_close(
      const char *restrict name,
      double expected,
      double actual,
      double tolerance);

static void test_Noble_pressure_validation(void) {
  volatile double press[] = {1.0, 0.0, -1.0, INFINITY, -INFINITY, NAN};
  const bool expected[]   = {true, false, false, false, false, false};

  for(size_t i=0; i<sizeof(press)/sizeof(press[0]); i++) {
    if(ghl_Noble_pressure_is_valid(press[i]) != expected[i]) {
      ghl_error("Noble finalized-pressure validation failed case %zu\n", i);
    }
  }
}

static void test_Noble_finalizer_speed_limit(void) {
  const ghl_con2prim_id_t backups[3] = {
    ghl_con2prim_id_None, ghl_con2prim_id_None, ghl_con2prim_id_None
  };
  ghl_parameters params;
  ghl_initialize_params(
        ghl_con2prim_id_Noble1D, backups, true, false, false, 1e100, 1.1, 0.0,
        &params);

  const double rho_ppoly[1] = { 0.0 };
  const double Gamma_ppoly[1] = { 2.0 };
  ghl_eos_parameters eos = { 0 };
  ghl_initialize_hybrid_eos_functions_and_params(
        1e-8, 1e-8, 1e6, 1, rho_ppoly, Gamma_ppoly, 0.7, 2.0, &eos);

  ghl_metric_quantities metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&metric, &metric_aux);

  ghl_conservative_quantities cons;
  ghl_initialize_conservatives(1.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, &cons);
  harm_aux_vars_struct harm_aux = { 0 };
  harm_aux.D = cons.rho;
  harm_aux.QU[1] = 10.0;

  ghl_primitive_quantities prims = { 0 };
  const double Z = 2.0;
  const double input_W = 2.0;
  const double input_vsq = 1.0 - 1.0/(input_W*input_W);
  if(!ghl_finalize_Noble(
        &params, &eos, &metric, &metric_aux, &cons, &harm_aux, Z, input_vsq,
        &prims)) {
    ghl_error("ordinary Noble finalizer did not trigger its speed limiter\n");
  }
  const double limited_W = metric.lapse*prims.u0;
  const double expected_rho = cons.rho/limited_W;
  const double expected_w = Z/(limited_W*limited_W);
  check_close("limited Noble rho", expected_rho, prims.rho, 1e-14);
  check_close(
        "limited Noble pressure", ghl_pressure_rho0_w(&eos, expected_rho, expected_w),
        prims.press, 1e-14);
  if(prims.press <= 0.0) {
    ghl_error("limited Noble closure test did not recover positive pressure\n");
  }
  const double recovered_Z = (prims.rho*(1.0 + prims.eps) + prims.press)
                           * limited_W*limited_W;
  check_close("limited Noble Z closure", Z, recovered_Z, 1e-14);

  prims = (ghl_primitive_quantities){ 0 };
  prims.rho = cons.rho/input_W;
  if(!ghl_finalize_Noble_entropy(
        &params, &eos, &metric, &metric_aux, &cons, &harm_aux, Z, input_W,
        &prims)) {
    ghl_error("entropy Noble finalizer did not trigger its speed limiter\n");
  }
  const double limited_entropy_W = metric.lapse*prims.u0;
  const double expected_entropy_rho = cons.rho/limited_entropy_W;
  const double expected_entropy_press
        = cons.entropy*pow(expected_entropy_rho, Gamma_ppoly[0] - 1.0)
        / limited_entropy_W;
  check_close("limited entropy Noble rho", expected_entropy_rho, prims.rho, 1e-14);
  check_close(
        "limited entropy Noble pressure", expected_entropy_press, prims.press, 1e-14);
}

static void check_close(
      const char *restrict name,
      const double expected,
      const double actual,
      const double tolerance) {

  const double scale = fmax(1.0, fabs(expected));
  if(!isfinite(actual) || fabs(actual - expected) > tolerance * scale) {
    ghl_error(
          "%s mismatch: expected %.17e, got %.17e\n", name, expected, actual);
  }
}

static void test_utilde_speed_limit_monotonicity(void) {
  const ghl_con2prim_id_t backups[3] = {
    ghl_con2prim_id_None, ghl_con2prim_id_None, ghl_con2prim_id_None
  };
  ghl_parameters params;
  ghl_initialize_params(
        ghl_con2prim_id_None, backups, false, false, false, 1e100, 2.0, 0.0,
        &params);
  ghl_metric_quantities metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);

  const double W_cases[3] = {1.9999998, 1.99999995, 2.0000001};
  for(int i = 0; i < 3; ++i) {
    double utU[3] = {sqrt(W_cases[i]*W_cases[i] - 1.0), 0.0, 0.0};
    const double input_ut = utU[0];
    ghl_primitive_quantities prims = { 0 };
    const bool limited = ghl_limit_utilde_and_compute_v(
          &params, &metric, utU, &prims);
    if((i < 2 && (limited || utU[0] != input_ut))
          || (i == 2 && (!limited || utU[0] >= input_ut))
          || utU[0] > input_ut) {
      ghl_error("utilde limiter was not monotone in boundary case %d\n", i);
    }
    const double output_W = metric.lapse*prims.u0;
    if(i == 2 && fabs(output_W - params.max_Lorentz_factor) > 1e-14) {
      ghl_error("utilde limiter did not reach W_max\n");
    }
  }
}

static void test_Font1D_roundtrips(void) {
  const ghl_con2prim_id_t None = ghl_con2prim_id_None;
  const ghl_con2prim_id_t backups[3] = { None, None, None };
  ghl_parameters params;
  ghl_initialize_params(
        ghl_con2prim_id_Font1D, backups, false, false, false, 1e100, 20.0, 0.0,
        &params);

  const double rho_ppoly[1] = { 0.0 };
  const double Gamma_ppoly[1] = { 2.0 };
  ghl_eos_parameters eos = { 0 };
  ghl_initialize_hybrid_eos_functions_and_params(
        1e-8, 1e-8, 1e6, 1, rho_ppoly, Gamma_ppoly, 0.7, 2.0, &eos);

  const double piecewise_rho_ppoly[2] = { 1.0, 0.0 };
  const double piecewise_Gamma_ppoly[2] = { 2.0, 2.5 };
  ghl_eos_parameters piecewise_eos = { 0 };
  ghl_initialize_hybrid_eos_functions_and_params(
        1e-8, 1e-8, 1e6, 2, piecewise_rho_ppoly,
        piecewise_Gamma_ppoly, 0.7, 2.0, &piecewise_eos);

  for(int test = 0; test < 5; ++test) {
    ghl_metric_quantities metric;
    if(test == 3) {
      ghl_initialize_metric(
            1.0, 0.0, 0.0, 0.0,
            1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
    }
    else {
      ghl_initialize_metric(
            0.91, 0.0, 0.0, 0.0,
            1.25, 0.04, 0.02, 1.12, -0.03, 0.97, &metric);
    }
    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&metric, &metric_aux);

    const ghl_eos_parameters *test_eos = test == 2 ? &piecewise_eos : &eos;
    const double rho = test == 0 ? 0.4 : (test == 4 ? 0.5*eos.rho_min : 1.3);
    double press, eps;
    ghl_hybrid_compute_P_cold_and_eps_cold(test_eos, rho, &press, &eps);
    ghl_primitive_quantities source;
    ghl_initialize_primitives(
          rho, press, eps,
          test == 0 ? -metric.betaU[0] : 0.16,
          test == 0 ? -metric.betaU[1] : -0.09,
          test == 0 ? -metric.betaU[2] : 0.06,
          test == 0 ? 0.0 : 0.25,
          test == 0 ? 0.0 : -0.17,
          test == 0 ? 0.0 : 0.11,
          0.0, 0.1, 0.0, &source);
    bool limited = false;
    ghl_abort_if_error(ghl_limit_v_and_compute_u0(&params, &metric, &source, &limited));
    if(limited) {
      ghl_error("Font1D source state was unexpectedly speed limited\n");
    }

    ghl_conservative_quantities cons, cons_undens;
    ghl_compute_conservs(&metric, &metric_aux, &source, &cons);
    ghl_undensitize_conservatives(metric.sqrt_detgamma, &cons, &cons_undens);

    ghl_primitive_quantities recovered = source;
    recovered.rho *= 1.2;
    recovered.press *= 0.8;
    recovered.eps *= 0.8;
    recovered.vU[0] *= 0.5;
    recovered.vU[1] *= 0.5;
    recovered.vU[2] *= 0.5;
    limited = false;
    ghl_abort_if_error(ghl_limit_v_and_compute_u0(&params, &metric, &recovered, &limited));

    ghl_con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);
    const ghl_error_codes_t error = ghl_hybrid_Font1D(
          &params, test_eos, &metric, &metric_aux, &cons_undens, &recovered, &diagnostics);
    if(error != ghl_success || diagnostics.which_routine != ghl_con2prim_id_Font1D
          || diagnostics.speed_limited || (test == 0 && diagnostics.n_iter != 0)
          || (test != 0 && diagnostics.n_iter < 1)) {
      ghl_error("Font1D independent roundtrip failed for case %d\n", test);
    }

    const double tolerance = test == 4
                           ? params.con2prim_solver_tolerance
                           : 1e-12;
    check_close("Font rho", source.rho, recovered.rho, tolerance);
    check_close("Font press", source.press, recovered.press, tolerance);
    check_close("Font eps", source.eps, recovered.eps, tolerance);
    for(int i = 0; i < 3; ++i) {
      check_close("Font velocity", source.vU[i], recovered.vU[i], tolerance);
      check_close("Font magnetic field", source.BU[i], recovered.BU[i], 0.0);
    }

    ghl_conservative_quantities recovered_cons, recovered_undens;
    ghl_compute_conservs(&metric, &metric_aux, &recovered, &recovered_cons);
    ghl_undensitize_conservatives(
          metric.sqrt_detgamma, &recovered_cons, &recovered_undens);
    check_close("Font conservative rho", cons_undens.rho, recovered_undens.rho, tolerance);
    check_close("Font conservative tau", cons_undens.tau, recovered_undens.tau, tolerance);
    for(int i = 0; i < 3; ++i) {
      check_close("Font conservative momentum", cons_undens.SD[i],
                  recovered_undens.SD[i], tolerance);
    }
  }

  // This finite, strongly magnetized state exhausts the first 300-iteration
  // attempt and converges after 319 iterations in the second attempt. The
  // diagnostic is the total over both attempts, not the final-attempt count.
  ghl_metric_quantities metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&metric, &metric_aux);
  ghl_conservative_quantities cons = {
    .rho = 5.9555223760409987e-4,
    .SD = {2.4062456840432316e5, 1.1216822303263594e1,
           -6.3928140990643856e-1},
  };
  ghl_primitive_quantities prims = {
    .BU = {2.0276492249373113e-5, 1.5166494025966969e-4,
           4.5941056264157982e2},
  };
  ghl_con2prim_diagnostics diagnostics;
  ghl_initialize_diagnostics(&diagnostics);
  const ghl_error_codes_t error = ghl_hybrid_Font1D(
        &params, &eos, &metric, &metric_aux, &cons, &prims, &diagnostics);
  if(error != ghl_success || diagnostics.which_routine != ghl_con2prim_id_Font1D
        || diagnostics.n_iter != 619) {
    ghl_error("Font1D multi-attempt iteration total mismatch: error=%d, routine=%d, n_iter=%d\n",
              (int)error, (int)diagnostics.which_routine, diagnostics.n_iter);
  }
}

static void check_entropy2_roundtrip(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict source,
      const ghl_conservative_quantities *restrict cons_undens) {

  ghl_primitive_quantities recovered = *source;
  recovered.rho *= 1.1;
  recovered.press *= 0.8;
  recovered.eps *= 0.9;
  for(int i = 0; i < 3; i++) {
    recovered.vU[i] *= 0.7;
  }
  bool speed_limited = false;
  ghl_error_codes_t error
        = ghl_limit_v_and_compute_u0(params, metric_adm, &recovered, &speed_limited);
  ghl_abort_if_error(error);
  const ghl_primitive_quantities initial_recovered = recovered;

  ghl_con2prim_diagnostics diagnostics;
  ghl_initialize_diagnostics(&diagnostics);
  error = ghl_con2prim_hybrid_select_method(
        ghl_con2prim_id_Noble1D_entropy2, params, eos, metric_adm, metric_aux,
        cons_undens, &recovered, &diagnostics);
  if(error != ghl_success) {
    ghl_error("Noble1D_entropy2 round trip failed with error %d\n", error);
  }
  if(diagnostics.which_routine != ghl_con2prim_id_Noble1D_entropy2
     || diagnostics.n_iter < 1 || diagnostics.speed_limited) {
    ghl_error("Noble1D_entropy2 returned inconsistent success diagnostics\n");
  }

  recovered = initial_recovered;
  ghl_initialize_diagnostics(&diagnostics);
  diagnostics.speed_limited = true;
  error = ghl_con2prim_hybrid_select_method(
        ghl_con2prim_id_Noble1D_entropy2, params, eos, metric_adm, metric_aux,
        cons_undens, &recovered, &diagnostics);
  if(error != ghl_success || !diagnostics.speed_limited) {
    ghl_error("Noble1D_entropy2 did not preserve an incoming speed-limit diagnostic\n");
  }

  check_close("rho", source->rho, recovered.rho, 1e-12);
  check_close("press", source->press, recovered.press, 1e-12);
  check_close("eps", source->eps, recovered.eps, 1e-12);
  check_close("entropy", source->entropy, recovered.entropy, 1e-12);
  for(int i = 0; i < 3; i++) {
    check_close("velocity", source->vU[i], recovered.vU[i], 1e-12);
    check_close("magnetic field", source->BU[i], recovered.BU[i], 0.0);
  }

  ghl_conservative_quantities recovered_cons;
  ghl_compute_conservs(metric_adm, metric_aux, &recovered, &recovered_cons);
  ghl_conservative_quantities recovered_cons_undens;
  ghl_undensitize_conservatives(
        metric_adm->sqrt_detgamma, &recovered_cons, &recovered_cons_undens);
  check_close("conservative rho", cons_undens->rho, recovered_cons_undens.rho, 1e-12);
  check_close("conservative tau", cons_undens->tau, recovered_cons_undens.tau, 1e-12);
  check_close(
        "conservative entropy", cons_undens->entropy, recovered_cons_undens.entropy,
        1e-12);
  for(int i = 0; i < 3; i++) {
    check_close(
          "conservative momentum", cons_undens->SD[i], recovered_cons_undens.SD[i],
          1e-12);
  }
}

static void test_Noble1D_entropy2(void) {
  const ghl_con2prim_id_t None = ghl_con2prim_id_None;
  const ghl_con2prim_id_t backups[3] = { None, None, None };
  ghl_parameters params;
  ghl_initialize_params(
        ghl_con2prim_id_Noble1D_entropy2, backups, true, false, false, 1e100, 10.0, 0.0,
        &params);
  params.con2prim_solver_tolerance = 1e-12;

  // Distinct cold and thermal Gammas make this a genuine hybrid-EOS test.
  const double rho_ppoly[2] = { 1.0, 0.0 };
  const double Gamma_ppoly[2] = { 1.8, 2.2 };
  ghl_eos_parameters eos = { 0 };
  ghl_initialize_hybrid_eos_functions_and_params(
        1e-6, 1e-6, 1e6, 2, rho_ppoly, Gamma_ppoly, 0.4, 1.6, &eos);

  for(int test = 0; test < 3; test++) {
    ghl_metric_quantities metric_adm;
    if(test == 0) {
      ghl_initialize_metric(
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric_adm);
    }
    else {
      ghl_initialize_metric(
            0.93, 0.02, -0.01, 0.015, 1.2, 0.04, 0.02, 1.1, 0.03, 0.95, &metric_adm);
    }
    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&metric_adm, &metric_aux);

    eos.rho_max = test == 2 ? 1.0 : 1e6;
    const double rho = test == 0 ? 0.3 : 3.0;
    double P_cold, eps_cold;
    ghl_hybrid_compute_P_cold_and_eps_cold(&eos, rho, &P_cold, &eps_cold);
    const double press = (test == 0 ? 1.35 : 1.8) * P_cold;
    const double eps = eps_cold + (press - P_cold) / ((eos.Gamma_th - 1.0) * rho);
    const double entropy = ghl_hybrid_compute_entropy_function(&eos, rho, press);
    const double vx = test == 0 ? 0.0 : 0.18;
    const double vy = test == 0 ? 0.0 : -0.11;
    const double vz = test == 0 ? 0.0 : 0.07;
    const double Bx = test == 0 ? 0.0 : 0.3;
    const double By = test == 0 ? 0.0 : -0.2;
    const double Bz = test == 0 ? 0.0 : 0.15;

    ghl_primitive_quantities source;
    ghl_initialize_primitives(
          rho, press, eps, vx, vy, vz, Bx, By, Bz, entropy, 0.1, 0.0, &source);
    bool speed_limited = false;
    ghl_error_codes_t error
          = ghl_limit_v_and_compute_u0(&params, &metric_adm, &source, &speed_limited);
    ghl_abort_if_error(error);
    if(speed_limited) {
      ghl_error("Noble1D_entropy2 source state was unexpectedly speed limited\n");
    }

    ghl_conservative_quantities cons;
    ghl_compute_conservs(&metric_adm, &metric_aux, &source, &cons);
    ghl_conservative_quantities cons_undens;
    ghl_undensitize_conservatives(metric_adm.sqrt_detgamma, &cons, &cons_undens);

    check_entropy2_roundtrip(
          &params, &eos, &metric_adm, &metric_aux, &source, &cons_undens);

    if(test == 0) {
      ghl_eos_parameters bounded_eos = eos;
      bounded_eos.rho_max = 0.1;
      ghl_primitive_quantities recovered = source;
      ghl_con2prim_diagnostics diagnostics;
      ghl_initialize_diagnostics(&diagnostics);
      error = ghl_hybrid_Noble1D_entropy2(
            &params, &bounded_eos, &metric_adm, &metric_aux, &cons_undens,
            &recovered, &diagnostics);
      if(error != ghl_success || recovered.rho <= bounded_eos.rho_max) {
        ghl_error("Noble1D_entropy2 did not return its unclamped density root\n");
      }
      bool limited = false;
      error = ghl_enforce_primitive_limits_and_compute_u0(
            &params, &bounded_eos, &metric_adm, &recovered, &limited);
      if(error != ghl_success || recovered.rho != bounded_eos.rho_max) {
        ghl_error("post-recovery limiter did not enforce the entropy2 density ceiling\n");
      }
    }

    if(test == 1) {
      ghl_parameters backup_params = params;
      backup_params.calc_prim_guess = false;
      backup_params.main_routine = ghl_con2prim_id_Noble1D_entropy2;
      backup_params.backup_routine[0] = ghl_con2prim_id_Noble1D;
      ghl_conservative_quantities negative_entropy = cons_undens;
      negative_entropy.entropy = -fabs(negative_entropy.entropy);

      ghl_primitive_quantities pristine_guess = source;
      pristine_guess.rho *= 1.05;
      pristine_guess.press *= 0.9;
      pristine_guess.eps *= 0.9;
      ghl_primitive_quantities failed_attempt = pristine_guess;
      ghl_con2prim_diagnostics diagnostics;
      ghl_initialize_diagnostics(&diagnostics);
      error = ghl_hybrid_Noble1D_entropy2(
            &backup_params, &eos, &metric_adm, &metric_aux, &negative_entropy,
            &failed_attempt, &diagnostics);
      if(error != ghl_error_neg_pressure
            || (failed_attempt.rho == pristine_guess.rho
                && failed_attempt.press == pristine_guess.press
                && failed_attempt.vU[0] == pristine_guess.vU[0])) {
        ghl_error("entropy2 failure did not provide a mutating retry test\n");
      }

      ghl_primitive_quantities expected = pristine_guess;
      ghl_initialize_diagnostics(&diagnostics);
      error = ghl_hybrid_Noble1D(
            &backup_params, &eos, &metric_adm, &metric_aux, &cons_undens,
            &expected, &diagnostics);
      if(error != ghl_success) {
        ghl_error("Noble1D reference backup failed with error %d\n", error);
      }
      const int expected_n_iter = diagnostics.n_iter;

      ghl_primitive_quantities recovered = pristine_guess;
      ghl_initialize_diagnostics(&diagnostics);
      error = ghl_con2prim_hybrid_multi_method(
            &backup_params, &eos, &metric_adm, &metric_aux, &negative_entropy,
            &recovered, &diagnostics);
      if(error != ghl_success || !diagnostics.backup[0] || diagnostics.backup[1]
            || diagnostics.backup[2]
            || diagnostics.which_routine != ghl_con2prim_id_Noble1D
            || diagnostics.n_iter != expected_n_iter) {
        ghl_error("hybrid backup restoration routing failed\n");
      }
      check_close("backup rho", expected.rho, recovered.rho, 2e-12);
      check_close("backup press", expected.press, recovered.press, 2e-12);
      check_close("backup eps", expected.eps, recovered.eps, 2e-12);
      for(int i = 0; i < 3; ++i) {
        check_close("backup velocity", expected.vU[i], recovered.vU[i], 2e-12);
        check_close("backup magnetic field", expected.BU[i], recovered.BU[i], 0.0);
      }

      ghl_primitive_quantities bad_guess = source;
      bad_guess.u0 = NAN;
      ghl_initialize_diagnostics(&diagnostics);
      error = ghl_con2prim_hybrid_select_method(
            ghl_con2prim_id_Noble1D_entropy2, &params, &eos, &metric_adm, &metric_aux,
            &cons_undens, &bad_guess, &diagnostics);
      if(error != ghl_error_invalid_utsq) {
        ghl_error("Noble1D_entropy2 did not propagate its initialization error\n");
      }

      ghl_conservative_quantities zero_density = cons_undens;
      zero_density.rho = 0.0;
      bad_guess = source;
      ghl_initialize_diagnostics(&diagnostics);
      error = ghl_con2prim_hybrid_select_method(
            ghl_con2prim_id_Noble1D_entropy2, &params, &eos, &metric_adm, &metric_aux,
            &zero_density, &bad_guess, &diagnostics);
      if(error != ghl_error_neg_rho) {
        ghl_error("Noble1D_entropy2 did not reject nonpositive density\n");
      }

      ghl_parameters iteration_params = params;
      iteration_params.con2prim_max_iterations = 1;
      iteration_params.con2prim_solver_tolerance = 1e-30;
      bad_guess = source;
      bad_guess.vU[0] *= 0.5;
      speed_limited = false;
      error = ghl_limit_v_and_compute_u0(
            &iteration_params, &metric_adm, &bad_guess, &speed_limited);
      ghl_abort_if_error(error);
      ghl_initialize_diagnostics(&diagnostics);
      error = ghl_con2prim_hybrid_select_method(
            ghl_con2prim_id_Noble1D_entropy2, &iteration_params, &eos, &metric_adm,
            &metric_aux, &cons_undens, &bad_guess, &diagnostics);
      if(error != ghl_error_c2p_max_iter) {
        ghl_error("Noble1D_entropy2 did not propagate its iteration failure\n");
      }

      ghl_conservative_quantities singular_entropy = cons_undens;
      singular_entropy.entropy = NAN;
      bad_guess = source;
      ghl_initialize_diagnostics(&diagnostics);
      error = ghl_con2prim_hybrid_select_method(
            ghl_con2prim_id_Noble1D_entropy2, &params, &eos, &metric_adm, &metric_aux,
            &singular_entropy, &bad_guess, &diagnostics);
      if(error != ghl_error_c2p_singular) {
        ghl_error("Noble1D_entropy2 did not propagate its singular solve failure\n");
      }

      bad_guess = source;
      ghl_initialize_diagnostics(&diagnostics);
      error = ghl_con2prim_hybrid_select_method(
            ghl_con2prim_id_Noble1D_entropy2, &params, &eos, &metric_adm, &metric_aux,
            &negative_entropy, &bad_guess, &diagnostics);
      if(error != ghl_error_neg_pressure) {
        ghl_error("Noble1D_entropy2 did not reject negative recovered pressure\n");
      }
    }
  }

  ghl_eos_parameters simple_eos = { 0 };
  ghl_initialize_simple_eos_functions_and_params(
        1e-6, 1e-6, 1e6, 1e-8, 1e-10, 1e6, 1.7, &simple_eos);
  ghl_metric_quantities simple_metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &simple_metric);
  ghl_ADM_aux_quantities simple_aux;
  ghl_compute_ADM_auxiliaries(&simple_metric, &simple_aux);
  ghl_primitive_quantities simple_source;
  const double simple_rho = 0.7;
  const double simple_press = 0.2;
  ghl_initialize_primitives(
        simple_rho, simple_press,
        simple_press / (simple_rho * (simple_eos.Gamma_th - 1.0)), 0.12, -0.04, 0.08,
        0.15, 0.07, -0.11,
        ghl_hybrid_compute_entropy_function(&simple_eos, simple_rho, simple_press), 0.1,
        0.0, &simple_source);
  bool simple_speed_limited = false;
  ghl_error_codes_t simple_error = ghl_limit_v_and_compute_u0(
        &params, &simple_metric, &simple_source, &simple_speed_limited);
  ghl_abort_if_error(simple_error);
  if(simple_speed_limited) {
    ghl_error("Noble1D_entropy2 simple-EOS source was speed limited\n");
  }
  ghl_conservative_quantities simple_cons, simple_cons_undens;
  ghl_compute_conservs(&simple_metric, &simple_aux, &simple_source, &simple_cons);
  ghl_undensitize_conservatives(
        simple_metric.sqrt_detgamma, &simple_cons, &simple_cons_undens);
  check_entropy2_roundtrip(
        &params, &simple_eos, &simple_metric, &simple_aux, &simple_source,
        &simple_cons_undens);

  const double high_W = 4.0;
  ghl_primitive_quantities high_W_source;
  ghl_initialize_primitives(
        simple_rho, simple_press,
        simple_press / (simple_rho * (simple_eos.Gamma_th - 1.0)),
        sqrt(1.0 - 1.0 / SQR(high_W)), 0.0, 0.0, 0.0, 0.0, 0.0,
        ghl_hybrid_compute_entropy_function(&simple_eos, simple_rho, simple_press), 0.1,
        0.0, &high_W_source);
  simple_speed_limited = false;
  simple_error = ghl_limit_v_and_compute_u0(
        &params, &simple_metric, &high_W_source, &simple_speed_limited);
  ghl_abort_if_error(simple_error);
  ghl_conservative_quantities high_W_cons, high_W_cons_undens;
  ghl_compute_conservs(&simple_metric, &simple_aux, &high_W_source, &high_W_cons);
  ghl_undensitize_conservatives(
        simple_metric.sqrt_detgamma, &high_W_cons, &high_W_cons_undens);
  ghl_con2prim_diagnostics high_W_diagnostics;
  ghl_initialize_diagnostics(&high_W_diagnostics);
  simple_error = ghl_hybrid_Noble1D_entropy2(
        &params, &simple_eos, &simple_metric, &simple_aux, &high_W_cons_undens,
        &high_W_source, &high_W_diagnostics);
  if(simple_error != ghl_success) {
    ghl_error("Noble1D_entropy2 rejected an exact W=4 guess\n");
  }
  check_close("exact W=4 rho", simple_rho, high_W_source.rho, 1e-12);
  check_close("exact W=4 press", simple_press, high_W_source.press, 1e-12);
  check_close("exact W=4 Lorentz factor", high_W, high_W_source.u0, 1e-12);

  // This malformed state has a positive-entropy momentum-equation root at
  // rho=2D, where the inferred velocity norm is negative.
  const double invalid_rho_ppoly[1] = { 0.0 };
  const double invalid_Gamma_ppoly[1] = { 2.0 };
  ghl_eos_parameters invalid_eos = { 0 };
  ghl_initialize_hybrid_eos_functions_and_params(
        1e-6, 1e-6, 1e6, 1, invalid_rho_ppoly, invalid_Gamma_ppoly, 1.0, 1.5,
        &invalid_eos);
  ghl_primitive_quantities invalid_guess;
  ghl_initialize_primitives(
        1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0 / 12.0, 0.1, 0.0,
        &invalid_guess);
  bool invalid_speed_limited = false;
  ghl_error_codes_t invalid_error = ghl_limit_v_and_compute_u0(
        &params, &simple_metric, &invalid_guess, &invalid_speed_limited);
  ghl_abort_if_error(invalid_error);
  ghl_conservative_quantities invalid_cons;
  ghl_initialize_conservatives(1.0, 0.0, 1.0, 0.0, 0.0, 1.0 / 12.0, 0.0, &invalid_cons);
  ghl_parameters invalid_params = params;
  invalid_params.con2prim_max_iterations = 100;
  ghl_con2prim_diagnostics invalid_diagnostics;
  ghl_initialize_diagnostics(&invalid_diagnostics);
  invalid_error = ghl_con2prim_hybrid_select_method(
        ghl_con2prim_id_Noble1D_entropy2, &invalid_params, &invalid_eos, &simple_metric,
        &simple_aux, &invalid_cons, &invalid_guess, &invalid_diagnostics);
  if(invalid_error != ghl_error_neg_vsq) {
    ghl_error("Noble1D_entropy2 did not reject a negative velocity norm\n");
  }
}

int main(int argc, char **argv) {

  const int arraylength = 2;

  ghl_error_codes_t *expected_errors = (ghl_error_codes_t*) malloc(sizeof(ghl_error_codes_t)*arraylength);
  expected_errors[0] = ghl_error_c2p_max_iter;
  expected_errors[1] = ghl_error_c2p_singular;

  const double poison = 0.0/0.0;

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const ghl_con2prim_id_t Noble2D = ghl_con2prim_id_Noble2D;
  const ghl_con2prim_id_t backup_routine[3] = {Noble2D, Noble2D, Noble2D};
  const bool evolve_entropy = false;
  const bool evolve_temperature = false;
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;

  const int neos = 1;
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(
        ghl_con2prim_id_Palenzuela1D, backup_routine,
        evolve_entropy, evolve_temperature, calc_prims_guess,
        Psi6threshold, W_max, 0.0, &params);

  ghl_eos_parameters eos = { 0 };
  ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_min, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &eos);

  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  double *gxx = (double*) malloc(sizeof(double)*arraylength);
  double *gxy = (double*) malloc(sizeof(double)*arraylength);
  double *gxz = (double*) malloc(sizeof(double)*arraylength);
  double *gyy = (double*) malloc(sizeof(double)*arraylength);
  double *gyz = (double*) malloc(sizeof(double)*arraylength);
  double *gzz = (double*) malloc(sizeof(double)*arraylength);

  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  double *rho_star = (double*) malloc(sizeof(double)*arraylength);
  double *tau = (double*) malloc(sizeof(double)*arraylength);
  double *S_x = (double*) malloc(sizeof(double)*arraylength);
  double *S_y = (double*) malloc(sizeof(double)*arraylength);
  double *S_z = (double*) malloc(sizeof(double)*arraylength);

  /*
     Hybrid_Noble2D failures
     1) rho=1e15: triggers failure to converge
     2) rho=0: nans as inputs to Newton-Raphson
  */

  for(int i=0; i<arraylength; i++) {
    lapse[i] = 1.0;
    betax[i] = 0.0;
    betay[i] = 0.0;
    betaz[i] = 0.0;
    gxx[i] = 1.0;
    gxy[i] = 0.0;
    gxz[i] = 0.0;
    gyy[i] = 1.0;
    gyz[i] = 0.0;
    gzz[i] = 1.0;
    Bx[i] = By[i] = Bz[i] = 0.0;
    rho_star[i] = 1e-2;
    tau[i]   = 1e-2;
    S_x[i] = S_y[i] = S_z[i] = 1000.0*tau[i]*(tau[i] + 2.0e-2);
  }
  rho_star[0] = 1e15;
  rho_star[1] = 0.0;

  for(int i=0; i<arraylength; i++) {
    ghl_con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);
    ghl_metric_quantities metric_adm;
    ghl_primitive_quantities prims;
    ghl_conservative_quantities cons, cons_undens;

    ghl_initialize_metric(lapse[i],
                      betax[i], betay[i], betaz[i],
                      gxx[i], gxy[i], gxz[i],
                      gyy[i], gyz[i], gzz[i],
                      &metric_adm);

    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&metric_adm, &metric_aux);

    ghl_initialize_primitives(
                        poison, poison, poison,
                        poison, poison, poison,
                        Bx[i], By[i], Bz[i],
                        poison, poison, poison, &prims);

    ghl_initialize_conservatives(rho_star[i], tau[i],
             S_x[i], S_y[i], S_z[i],
             poison, poison, &cons);

    ghl_undensitize_conservatives(metric_adm.sqrt_detgamma, &cons, &cons_undens);
    int check = ghl_con2prim_hybrid_multi_method(&params, &eos, &metric_adm, &metric_aux, &cons_undens, &prims, &diagnostics);
    if(check != expected_errors[i]
          || diagnostics.which_routine != ghl_con2prim_id_None)
      ghl_error(
            "hybrid all-method failure contract failed: expected error %d/routine %d, "
            "got error %d/routine %d",
            expected_errors[i], ghl_con2prim_id_None, check,
            diagnostics.which_routine);

    if(i == 1) {
      ghl_initialize_diagnostics(&diagnostics);
      check = ghl_hybrid_Palenzuela1D_entropy(
            &params, &eos, &metric_adm, &metric_aux, &cons_undens, &prims,
            &diagnostics);
      if(check == ghl_success
            || diagnostics.which_routine != ghl_con2prim_id_None) {
        ghl_error(
              "hybrid Palenzuela1D entropy failure set routine %d (error=%d)",
              diagnostics.which_routine, check);
      }
    }

  }
  test_Font1D_roundtrips();
  test_Noble_pressure_validation();
  test_Noble_finalizer_speed_limit();
  test_utilde_speed_limit_monotonicity();
  test_Noble1D_entropy2();
  return 0;
}
