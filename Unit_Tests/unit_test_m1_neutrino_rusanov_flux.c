#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit.h"
#include "m1_neutrino_seeded_test_utils.h"
#include "m1_thcm1_rusanov_fixture.h"

/*
 * Boundary and reference coverage for the grey neutrino Rusanov path. The
 * test checks number-current and physical-flux construction, validation and
 * transactional failures, plus seeded and retained trusted/perturbed pairs.
 */

static void fail_case(const char *message, const int case_index) {
  ghl_error("unit_test_m1_neutrino_rusanov_flux case %d: %s\n", case_index, message);
}

static void check_close(
      const double actual,
      const double expected,
      const char *label,
      const int case_index) {
  if(!m1_nearly_equal(actual, expected, 3.0e-12, 2.0e-13)) {
    fail_case(label, case_index);
  }
}

static void check_number_product_overflow(const ghl_m1_parameters *params) {
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  metric.gammaDD[0][0] = 0.01;
  metric.gammaDD[1][1] = metric.gammaDD[2][2] = 10.0;
  metric.gammaUU[0][0] = 100.0;
  metric.gammaUU[1][1] = metric.gammaUU[2][2] = 0.1;
  const ghl_m1_closure closure = {
    .chi = 1.0 / 3.0,
    .P = { { 100.0 / 3.0, 0.0, 0.0 }, { 0.0, 0.1 / 3.0, 0.0 }, { 0.0, 0.0, 0.1 / 3.0 } }
  };
  const ghl_m1_neutrino_parameters nu = { .N_floor = 0.0 };
  /* A coordinate speed of two is subluminal in this metric. The number
   * product still overflows for finite N=DBL_MAX on either face side. */
  for(int side = 0; side < 2; ++side) {
    ghl_m1_neutrino_state state[2] = { { .N = 1.0, .E = 1.0 }, { .N = 1.0, .E = 1.0 } };
    double velocity[2][3] = { { 0.0 } };
    double number[2][3] = { { 0.0 } };
    state[side].N = DBL_MAX;
    velocity[side][0] = 2.0;
    double N = -1.0, E = -2.0, F[3] = { -3.0, -4.0, -5.0 };
    if(ghl_m1_compute_neutrino_rusanov_flux(
             params, &nu, &metric, ghl_m1_dirn0, &state[0], &state[1], &closure,
             &closure, number[0], number[1], velocity[0], velocity[1], 1.0, &N, &E, F)
             != ghl_error_m1_invalid_state
       || N != -1.0 || E != -2.0 || F[0] != -3.0 || F[1] != -4.0 || F[2] != -5.0) {
      fail_case("finite number-product overflow did not reject transactionally", side);
    }
  }
}

static void check_neutrino_number_transport_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims) {
  const ghl_m1_neutrino_state state = { .N = 1.0, .E = 1.0, .F = { 0.0, 0.0, 0.0 } };
  const ghl_m1_comoving comoving
        = { .J = 1.0, .HU = { 0.0, 0.0, 0.0 }, .HD = { 0.0, 0.0, 0.0 }, .Hn = 0.0 };
  const double V_con[3] = { 0.0, 0.0, 0.0 };
  ghl_m1_neutrino_current current
        = { .J = 91.0,
            .h_n = 92.0,
            .HU = { 93.0, 94.0, 95.0 },
            .Gamma_N = 96.0,
            .n_com = 97.0,
            .number_flux = { 98.0, 99.0, 100.0 },
            .number_transport_velocity = { 101.0, 102.0, 103.0 } };

  if(ghl_m1_neutrino_build_current_from_moments(
           NULL, nu_params, &state, &comoving, V_con, 1.0, &current)
           != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_build_current_from_moments(
              metric, NULL, &state, &comoving, V_con, 1.0, &current)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_build_current_from_moments(
              metric, nu_params, NULL, &comoving, V_con, 1.0, &current)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_build_current_from_moments(
              metric, nu_params, &state, NULL, V_con, 1.0, &current)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_build_current_from_moments(
              metric, nu_params, &state, &comoving, NULL, 1.0, &current)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_build_current_from_moments(
              metric, nu_params, &state, &comoving, V_con, 1.0, NULL)
              != ghl_error_m1_null_pointer) {
    fail_case("number-current moment NULL boundary failed", 300);
  }

  ghl_m1_neutrino_parameters bad_nu = *nu_params;
  bad_nu.N_floor = NAN;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, &bad_nu, &state, &comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite N floor was accepted", 301);
  }
  bad_nu.N_floor = -1.0;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, &bad_nu, &state, &comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("negative N floor was accepted", 302);
  }
  bad_nu = *nu_params;
  bad_nu.J_floor = NAN;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, &bad_nu, &state, &comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite J floor was accepted", 303);
  }
  bad_nu.J_floor = -1.0;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, &bad_nu, &state, &comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("negative J floor was accepted", 304);
  }
  bad_nu = *nu_params;
  bad_nu.Gamma_N_floor = NAN;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, &bad_nu, &state, &comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite Gamma floor was accepted", 305);
  }
  bad_nu.Gamma_N_floor = -1.0;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, &bad_nu, &state, &comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("negative Gamma floor was accepted", 306);
  }
  ghl_m1_neutrino_state bad_state = state;
  bad_state.N = NAN;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &bad_state, &comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite neutrino number was accepted", 307);
  }
  bad_nu = *nu_params;
  bad_nu.N_floor = 2.0;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, &bad_nu, &state, &comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("number below the configured floor was accepted", 308);
  }

  ghl_m1_comoving bad_comoving = comoving;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, V_con, NAN, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite Lorentz factor was accepted", 309);
  }
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, V_con, 0.5, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("subunit Lorentz factor was accepted", 310);
  }
  bad_comoving = comoving;
  bad_comoving.J = NAN;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite comoving energy was accepted", 311);
  }
  bad_comoving.J = 0.0;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonpositive comoving energy was accepted", 312);
  }
  bad_comoving = comoving;
  bad_comoving.Hn = NAN;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite comoving number flux was accepted", 313);
  }
  bad_comoving = comoving;
  bad_comoving.HU[0] = NAN;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite comoving spatial flux was accepted", 314);
  }

  bad_comoving = comoving;
  bad_comoving.J = DBL_MIN;
  bad_comoving.Hn = DBL_MAX;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite Gamma_N was accepted", 315);
  }
  bad_comoving = comoving;
  bad_comoving.Hn = 1.0;
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, V_con, 1.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonpositive Gamma_N with N>0 was accepted", 316);
  }
  bad_state = state;
  bad_state.N = 0.0;
  ghl_m1_neutrino_parameters zero_nu = *nu_params;
  zero_nu.N_floor = 0.0;
  const ghl_error_codes_t zero_gamma_error = ghl_m1_neutrino_build_current_from_moments(
        metric, &zero_nu, &bad_state, &bad_comoving, V_con, 1.0, &current);
  if(zero_gamma_error != ghl_success) {
    fail_case("zero-N singular Gamma_N returned an error", 317);
  }
  if(current.n_com != 0.0 || current.number_flux[0] != 0.0
     || current.number_flux[1] != 0.0 || current.number_flux[2] != 0.0
     || current.number_transport_velocity[0] != 0.0
     || current.number_transport_velocity[1] != 0.0
     || current.number_transport_velocity[2] != 0.0) {
    fail_case("zero-N singular Gamma_N policy failed", 317);
  }

  bad_comoving = comoving;
  double huge_V[3] = { DBL_MAX, 0.0, 0.0 };
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, huge_V, DBL_MAX, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite number transport velocity was accepted", 318);
  }
  bad_state = state;
  bad_comoving = comoving;
  bad_comoving.Hn = 1.0;
  bad_state.N = DBL_MAX;
  double large_V[3] = { 1.0, 0.0, 0.0 };
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &bad_state, &bad_comoving, large_V, 3.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite number flux was accepted", 319);
  }
  if(ghl_m1_neutrino_build_current_from_moments(
           metric, nu_params, &state, &bad_comoving, large_V, 3.0, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("superluminal number transport velocity was accepted", 320);
  }

  ghl_m1_closure closure = { 0 };
  if(ghl_m1_compute_neutrino_closure(m1_params, metric, prims, &state, &closure)
     != ghl_success) {
    fail_case("number-flux boundary closure construction failed", 321);
  }

  double number_flux[3] = { 111.0, 112.0, 113.0 };
  double number_velocity[3] = { 114.0, 115.0, 116.0 };
  const double number_flux_before[3] = { 111.0, 112.0, 113.0 };
  const double number_velocity_before[3] = { 114.0, 115.0, 116.0 };
  if(ghl_m1_compute_neutrino_number_flux(
           m1_params, nu_params, metric, prims, &state, NULL, number_velocity)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_number_flux(
              m1_params, nu_params, metric, prims, &state, number_flux, NULL)
              != ghl_error_m1_null_pointer) {
    fail_case("number-flux output NULL boundary failed", 322);
  }
  bad_state = state;
  bad_state.F[0] = 2.0;
  if(ghl_m1_compute_neutrino_number_flux(
           m1_params, nu_params, metric, prims, &bad_state, number_flux, number_velocity)
           != ghl_error_m1_invalid_state
     || memcmp(number_flux, number_flux_before, sizeof(number_flux)) != 0
     || memcmp(number_velocity, number_velocity_before, sizeof(number_velocity)) != 0) {
    fail_case("number-flux invalid state was not transactional", 323);
  }

  if(ghl_m1_neutrino_derive_current(NULL, nu_params, metric, prims, &state, &current)
           != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_derive_current(m1_params, NULL, metric, prims, &state, &current)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_derive_current(
              m1_params, nu_params, NULL, prims, &state, &current)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_derive_current(
              m1_params, nu_params, metric, NULL, &state, &current)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_derive_current(m1_params, nu_params, metric, prims, &state, NULL)
              != ghl_error_m1_null_pointer) {
    fail_case("derived-current NULL boundary failed", 324);
  }
  ghl_metric_quantities bad_metric = *metric;
  bad_metric.gammaDD[0][0] = -1.0;
  if(ghl_m1_neutrino_derive_current(
           m1_params, nu_params, &bad_metric, prims, &state, &current)
     != ghl_error_m1_invalid_metric) {
    fail_case("derived-current invalid metric was accepted", 325);
  }
  bad_nu = *nu_params;
  bad_nu.N_floor = NAN;
  if(ghl_m1_neutrino_derive_current(m1_params, &bad_nu, metric, prims, &state, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("derived-current nonfinite N floor was accepted", 326);
  }
  bad_nu.N_floor = -1.0;
  if(ghl_m1_neutrino_derive_current(m1_params, &bad_nu, metric, prims, &state, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("derived-current negative N floor was accepted", 327);
  }
  bad_nu = *nu_params;
  bad_nu.J_floor = NAN;
  if(ghl_m1_neutrino_derive_current(m1_params, &bad_nu, metric, prims, &state, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("derived-current nonfinite J floor was accepted", 328);
  }
  bad_nu.J_floor = -1.0;
  if(ghl_m1_neutrino_derive_current(m1_params, &bad_nu, metric, prims, &state, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("derived-current negative J floor was accepted", 329);
  }
  bad_nu = *nu_params;
  bad_nu.Gamma_N_floor = NAN;
  if(ghl_m1_neutrino_derive_current(m1_params, &bad_nu, metric, prims, &state, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("derived-current nonfinite Gamma floor was accepted", 330);
  }
  bad_nu.Gamma_N_floor = -1.0;
  if(ghl_m1_neutrino_derive_current(m1_params, &bad_nu, metric, prims, &state, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("derived-current negative Gamma floor was accepted", 331);
  }
  bad_state = state;
  bad_state.N = -1.0;
  if(ghl_m1_neutrino_derive_current(
           m1_params, nu_params, metric, prims, &bad_state, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("derived-current negative N was accepted", 332);
  }
  ghl_primitive_quantities bad_prims = *prims;
  bad_prims.vU[0] = NAN;
  if(ghl_m1_neutrino_derive_current(
           m1_params, nu_params, metric, &bad_prims, &state, &current)
     != ghl_error_m1_invalid_state) {
    fail_case("derived-current closure failure was not propagated", 333);
  }

  ghl_m1_closure bad_closure = closure;
  bad_closure.P[0][0] = NAN;
  if(ghl_m1_neutrino_derive_current_from_closure(
           m1_params, nu_params, metric, prims, &state, NULL, &current)
           != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_derive_current_from_closure(
              m1_params, nu_params, metric, prims, &state, &closure, NULL)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_derive_current_from_closure(
              m1_params, nu_params, metric, prims, &state, &bad_closure, &current)
              != ghl_error_m1_invalid_state) {
    fail_case("closure-supplied current boundary failed", 334);
  }

  current = (ghl_m1_neutrino_current){ .number_flux = { 0.25, -0.1, 0.05 },
                                       .number_transport_velocity = { 0.0, 0.0, 0.0 } };
  double physical_number = 117.0;
  if(ghl_m1_neutrino_physical_number_flux_from_current(
           NULL, &state, &current, ghl_m1_dirn0, &physical_number)
           != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_physical_number_flux_from_current(
              metric, NULL, &current, ghl_m1_dirn0, &physical_number)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_physical_number_flux_from_current(
              metric, &state, NULL, ghl_m1_dirn0, &physical_number)
              != ghl_error_m1_null_pointer
     || ghl_m1_neutrino_physical_number_flux_from_current(
              metric, &state, &current, ghl_m1_dirn0, NULL)
              != ghl_error_m1_null_pointer) {
    fail_case("physical number-current NULL boundary failed", 335);
  }
  if(ghl_m1_neutrino_physical_number_flux_from_current(
           metric, &state, &current, (ghl_m1_direction_t)3, &physical_number)
     != ghl_error_m1_invalid_state) {
    fail_case("invalid physical number direction was accepted", 336);
  }
  if(ghl_m1_neutrino_physical_number_flux_from_current(
           &bad_metric, &state, &current, ghl_m1_dirn0, &physical_number)
     != ghl_error_m1_invalid_metric) {
    fail_case("invalid physical number metric was accepted", 337);
  }
  bad_state = state;
  bad_state.N = NAN;
  if(ghl_m1_neutrino_physical_number_flux_from_current(
           metric, &bad_state, &current, ghl_m1_dirn0, &physical_number)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite physical number state was accepted", 338);
  }
  ghl_m1_neutrino_current bad_current = current;
  bad_current.number_flux[0] = NAN;
  if(ghl_m1_neutrino_physical_number_flux_from_current(
           metric, &state, &bad_current, ghl_m1_dirn0, &physical_number)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite physical number flux was accepted", 339);
  }
  bad_current = current;
  bad_current.number_transport_velocity[0] = NAN;
  if(ghl_m1_neutrino_physical_number_flux_from_current(
           metric, &state, &bad_current, ghl_m1_dirn0, &physical_number)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite physical number velocity was accepted", 340);
  }
  bad_current = current;
  bad_current.number_transport_velocity[0] = 2.0;
  if(ghl_m1_neutrino_physical_number_flux_from_current(
           metric, &state, &bad_current, ghl_m1_dirn0, &physical_number)
     != ghl_error_m1_invalid_state) {
    fail_case("superluminal physical number velocity was accepted", 341);
  }
  ghl_metric_quantities huge_lapse_metric = *metric;
  huge_lapse_metric.lapse = DBL_MAX;
  current.number_flux[0] = 2.0;
  if(ghl_m1_neutrino_physical_number_flux_from_current(
           &huge_lapse_metric, &state, &current, ghl_m1_dirn0, &physical_number)
     != ghl_error_m1_invalid_state) {
    fail_case("overflowing physical number flux was accepted", 342);
  }

  ghl_metric_quantities shifted_metric = *metric;
  shifted_metric.lapse = 4.0;
  shifted_metric.betaU[0] = 2.0;
  const ghl_m1_neutrino_state large_number_state = { .N = 1.0e308, .E = 1.0 };
  const ghl_m1_neutrino_current large_number_current
        = { .number_flux = { 5.0e307, 0.0, 0.0 },
            .number_transport_velocity = { 0.5, 0.0, 0.0 } };
  physical_number = 119.0;
  if(ghl_m1_neutrino_physical_number_flux_from_current(
           &shifted_metric, &large_number_state, &large_number_current, ghl_m1_dirn0,
           &physical_number)
           != ghl_success
     || physical_number != 0.0) {
    fail_case("finite cancelled physical number flux was rejected", 342);
  }

  double public_number_flux[3] = { 121.0, 122.0, 123.0 };
  double public_number_flux_before[3] = { 121.0, 122.0, 123.0 };
  if(ghl_m1_compute_neutrino_number_flux_from_closure(
           m1_params, nu_params, metric, prims, &state, NULL, public_number_flux,
           number_velocity)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_number_flux_from_closure(
              m1_params, nu_params, metric, prims, &state, &closure, NULL,
              number_velocity)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_number_flux_from_closure(
              m1_params, nu_params, metric, prims, &state, &closure, public_number_flux,
              NULL)
              != ghl_error_m1_null_pointer) {
    fail_case("closure number-flux output NULL boundary failed", 343);
  }
  if(ghl_m1_compute_neutrino_number_flux_from_closure(
           m1_params, nu_params, metric, prims, &state, &bad_closure, public_number_flux,
           number_velocity)
           != ghl_error_m1_invalid_state
     || memcmp(public_number_flux, public_number_flux_before, sizeof(public_number_flux))
              != 0) {
    fail_case("closure number-flux invalid input was not transactional", 344);
  }

  if(ghl_m1_compute_neutrino_physical_number_flux(
           NULL, nu_params, metric, prims, &state, ghl_m1_dirn0, &physical_number)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_physical_number_flux(
              m1_params, nu_params, NULL, prims, &state, ghl_m1_dirn0, &physical_number)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_physical_number_flux(
              m1_params, nu_params, metric, prims, NULL, ghl_m1_dirn0, &physical_number)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_physical_number_flux(
              m1_params, nu_params, metric, prims, &state, ghl_m1_dirn0, NULL)
              != ghl_error_m1_null_pointer) {
    fail_case("public physical number NULL boundary failed", 345);
  }
  bad_state = state;
  bad_state.F[0] = 2.0;
  if(ghl_m1_compute_neutrino_physical_number_flux(
           m1_params, nu_params, metric, prims, &bad_state, ghl_m1_dirn0,
           &physical_number)
     != ghl_error_m1_invalid_state) {
    fail_case("public physical number invalid state was accepted", 346);
  }
  if(ghl_m1_compute_neutrino_physical_number_flux_from_closure(
           m1_params, nu_params, metric, prims, &state, NULL, ghl_m1_dirn0,
           &physical_number)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_physical_number_flux_from_closure(
              m1_params, nu_params, metric, prims, &state, &closure, ghl_m1_dirn0, NULL)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_physical_number_flux_from_closure(
              m1_params, nu_params, metric, prims, &state, &bad_closure, ghl_m1_dirn0,
              &physical_number)
              != ghl_error_m1_invalid_state) {
    fail_case("closure physical number boundary failed", 347);
  }
}

static void check_neutrino_rusanov_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims) {
  const ghl_m1_neutrino_state state_L = { .N = 1.0, .E = 1.0, .F = { 0.0, 0.0, 0.0 } };
  const ghl_m1_neutrino_state state_R = { .N = 1.5, .E = 1.2, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_closure closure_L = { 0 };
  ghl_m1_closure closure_R = { 0 };
  if(ghl_m1_compute_neutrino_closure(m1_params, metric, prims, &state_L, &closure_L)
           != ghl_success
     || ghl_m1_compute_neutrino_closure(m1_params, metric, prims, &state_R, &closure_R)
              != ghl_success) {
    fail_case("neutrino Rusanov boundary closure construction failed", 350);
  }

  const double number_flux_L[3] = { 0.0, 0.0, 0.0 };
  const double number_flux_R[3] = { 0.0, 0.0, 0.0 };
  const double number_velocity_L[3] = { 0.0, 0.0, 0.0 };
  const double number_velocity_R[3] = { 0.0, 0.0, 0.0 };
  double flux_N = 131.0, flux_E = 132.0;
  double flux_F[3] = { 133.0, 134.0, 135.0 };
  if(ghl_m1_compute_neutrino_rusanov_flux(
           NULL, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, NULL, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, NULL, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, NULL, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, NULL, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, NULL,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              NULL, number_flux_L, number_flux_R, number_velocity_L, number_velocity_R,
              0.5, &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, NULL, number_flux_R, number_velocity_L, number_velocity_R, 0.5,
              &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, NULL, number_velocity_L, number_velocity_R, 0.5,
              &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, NULL, number_velocity_R, 0.5,
              &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L, NULL, 0.5,
              &flux_N, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, 0.5, NULL, &flux_E, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, 0.5, &flux_N, NULL, flux_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, 0.5, &flux_N, &flux_E, NULL)
              != ghl_error_m1_null_pointer) {
    fail_case("neutrino Rusanov NULL boundary failed", 351);
  }

  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, (ghl_m1_direction_t)3, &state_L, &state_R,
           &closure_L, &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
           != ghl_error_m1_invalid_state
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, NAN, &flux_N, &flux_E, flux_F)
              != ghl_error_m1_invalid_state
     || ghl_m1_compute_neutrino_rusanov_flux(
              m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
              &closure_R, number_flux_L, number_flux_R, number_velocity_L,
              number_velocity_R, -1.0, &flux_N, &flux_E, flux_F)
              != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov direction/speed validation failed", 352);
  }
  ghl_metric_quantities bad_metric = *metric;
  bad_metric.gammaDD[0][0] = -1.0;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, &bad_metric, ghl_m1_dirn0, &state_L, &state_R,
           &closure_L, &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_metric) {
    fail_case("neutrino Rusanov invalid metric was accepted", 353);
  }
  ghl_m1_neutrino_parameters bad_nu = *nu_params;
  bad_nu.N_floor = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, &bad_nu, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov nonfinite N floor was accepted", 354);
  }
  bad_nu.N_floor = -1.0;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, &bad_nu, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov negative N floor was accepted", 355);
  }

  ghl_m1_neutrino_state bad_state = state_L;
  bad_state.F[0] = 2.0;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &bad_state, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov left realizability failure was lost", 356);
  }
  bad_state = state_R;
  bad_state.F[0] = 2.0;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &bad_state, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov right realizability failure was lost", 357);
  }
  ghl_m1_closure bad_closure = closure_L;
  bad_closure.P[0][0] = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &bad_closure,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov left closure failure was lost", 358);
  }
  bad_closure = closure_R;
  bad_closure.P[0][0] = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &bad_closure, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov right closure failure was lost", 359);
  }

  bad_state = state_L;
  bad_state.N = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &bad_state, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov left nonfinite N was accepted", 360);
  }
  bad_state = state_R;
  bad_state.N = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &bad_state, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov right nonfinite N was accepted", 361);
  }
  bad_state = state_L;
  bad_state.N = -1.0;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &bad_state, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov left floor violation was accepted", 362);
  }
  bad_state = state_R;
  bad_state.N = -1.0;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &bad_state, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("neutrino Rusanov right floor violation was accepted", 363);
  }

  double bad_number_flux[3] = { NAN, 0.0, 0.0 };
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, bad_number_flux, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite left number flux was accepted", 364);
  }
  bad_number_flux[0] = 0.0;
  bad_number_flux[1] = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, bad_number_flux, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite second left number flux was accepted", 365);
  }
  bad_number_flux[1] = 0.0;
  bad_number_flux[2] = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, bad_number_flux, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite third left number flux was accepted", 366);
  }
  bad_number_flux[2] = 0.0;
  bad_number_flux[0] = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, bad_number_flux, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite right number flux was accepted", 367);
  }
  bad_number_flux[0] = 0.0;
  bad_number_flux[1] = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, bad_number_flux, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite second right number flux was accepted", 368);
  }
  bad_number_flux[1] = 0.0;
  bad_number_flux[2] = NAN;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, bad_number_flux, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite third right number flux was accepted", 369);
  }
  double bad_velocity[3] = { NAN, 0.0, 0.0 };
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, bad_velocity, number_velocity_R,
           0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite left number velocity was accepted", 371);
  }
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L, bad_velocity,
           0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("nonfinite right number velocity was accepted", 372);
  }
  bad_velocity[0] = DBL_MAX;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, bad_velocity, number_velocity_R,
           0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("superluminal left number velocity was accepted", 373);
  }
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, number_flux_L, number_flux_R, number_velocity_L, bad_velocity,
           0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("superluminal right number velocity was accepted", 374);
  }

  double inconsistent_flux[3] = { 1.0, 0.0, 0.0 };
  const double output_before[5] = { 141.0, 142.0, 143.0, 144.0, 145.0 };
  flux_N = output_before[0];
  flux_E = output_before[1];
  memcpy(flux_F, output_before + 2, sizeof(flux_F));
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &state_R, &closure_L,
           &closure_R, inconsistent_flux, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
           != ghl_error_m1_invalid_state
     || flux_N != output_before[0] || flux_E != output_before[1]
     || memcmp(flux_F, output_before + 2, sizeof(flux_F)) != 0) {
    fail_case("positive-N number-current mismatch was not transactional", 375);
  }
  ghl_m1_neutrino_state zero_state = state_L;
  zero_state.N = 0.0;
  ghl_m1_neutrino_parameters zero_nu = *nu_params;
  zero_nu.N_floor = 0.0;
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, &zero_nu, metric, ghl_m1_dirn0, &zero_state, &zero_state,
           &closure_L, &closure_L, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_success) {
    fail_case("zero-N zero-current Rusanov case failed", 376);
  }
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, &zero_nu, metric, ghl_m1_dirn0, &zero_state, &state_R, &closure_L,
           &closure_R, inconsistent_flux, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("zero-N nonzero current mismatch was accepted", 377);
  }

  ghl_metric_quantities huge_lapse_metric = *metric;
  huge_lapse_metric.lapse = DBL_MAX;
  ghl_m1_neutrino_state huge_number_state = state_L;
  huge_number_state.N = DBL_MAX;
  const double huge_number_flux[3] = { DBL_MAX, 0.0, 0.0 };
  const double unit_number_velocity[3] = { 1.0, 0.0, 0.0 };
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, &huge_lapse_metric, ghl_m1_dirn0, &huge_number_state,
           &state_R, &closure_L, &closure_R, huge_number_flux, number_flux_R,
           unit_number_velocity, number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("left physical number overflow was accepted", 378);
  }
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, &huge_lapse_metric, ghl_m1_dirn0, &state_R,
           &huge_number_state, &closure_R, &closure_L, number_flux_R, huge_number_flux,
           number_velocity_R, unit_number_velocity, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("right physical number overflow was accepted", 379);
  }

  const ghl_m1_neutrino_state high_state
        = { .N = 0.0, .E = DBL_MAX, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_closure high_closure = { 0 };
  if(ghl_m1_compute_neutrino_closure(
           m1_params, metric, prims, &high_state, &high_closure)
     != ghl_success) {
    fail_case("high-energy neutrino closure construction failed", 380);
  }
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, &huge_lapse_metric, ghl_m1_dirn0, &high_state, &state_R,
           &high_closure, &closure_R, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("left physical E/F overflow was accepted", 381);
  }
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, &huge_lapse_metric, ghl_m1_dirn0, &state_L, &high_state,
           &closure_L, &high_closure, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, 0.5, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("right physical E/F overflow was accepted", 382);
  }
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, metric, ghl_m1_dirn0, &state_L, &high_state, &closure_L,
           &high_closure, number_flux_L, number_flux_R, number_velocity_L,
           number_velocity_R, DBL_MAX, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("Rusanov candidate overflow was accepted", 383);
  }

  ghl_metric_quantities large_metric = { 0 };
  const double diagonal = 1.0e100;
  large_metric.lapse = 1.0;
  large_metric.lapseinv = 1.0;
  large_metric.lapseinv2 = 1.0;
  large_metric.detgamma = diagonal * diagonal * diagonal;
  large_metric.sqrt_detgamma = sqrt(large_metric.detgamma);
  for(int i = 0; i < 3; ++i) {
    large_metric.gammaDD[i][i] = diagonal;
    large_metric.gammaUU[i][i] = 1.0 / diagonal;
  }
  const ghl_m1_neutrino_state large_state
        = { .N = 0.0, .E = 1.0e200, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_closure large_closure = { 0 };
  if(ghl_m1_compute_neutrino_closure(
           m1_params, &large_metric, prims, &large_state, &large_closure)
     != ghl_success) {
    fail_case("large-metric closure construction failed", 384);
  }
  if(ghl_m1_compute_neutrino_rusanov_flux(
           m1_params, nu_params, &large_metric, ghl_m1_dirn0, &large_state, &large_state,
           &large_closure, &large_closure, number_flux_L, number_flux_R,
           number_velocity_L, number_velocity_R, 0.0, &flux_N, &flux_E, flux_F)
     != ghl_error_m1_invalid_state) {
    fail_case("densitized Rusanov overflow was accepted", 385);
  }
}


static void check_rusanov_arithmetic_boundaries(const ghl_m1_parameters *params) {
  for(int scenario=0; scenario<7; ++scenario) {
    ghl_metric_quantities metric;
    m1_setup_flat_metric(&metric);
    ghl_m1_neutrino_parameters nu = {.N_floor=0};
    ghl_m1_neutrino_state state[2] = {{.N=1, .E=1}, {.N=1, .E=1}};
    ghl_m1_closure closure[2] = {
      {.chi=1.0/3, .P={{1.0/3,0,0},{0,1.0/3,0},{0,0,1.0/3}}},
      {.chi=1.0/3, .P={{1.0/3,0,0},{0,1.0/3,0},{0,0,1.0/3}}}
    };
    double number[2][3]={{0}}, velocity[2][3]={{0}};
    double speed=1;
    if(scenario==0) {state[1].N=0; number[1][0]=1;}
    if(scenario==1) {velocity[1][0]=0.5; number[1][0]=0.25;}
    if(scenario==2 || scenario==3) {
      /* E/F physical flux overflows on exactly one side, while N flux is zero. */
      const int side=scenario-2;
      metric.lapse=DBL_MAX;
      state[side].E=8;
      for(int i=0;i<3;++i) closure[side].P[i][i]=8.0/3;
    }
    if(scenario==4) {state[1].N=DBL_MAX; speed=4;}
    if(scenario==5) {
      /* Each physical flux is finite; densitization alone overflows. */
      metric.gammaDD[0][0]=metric.gammaDD[1][1]=metric.gammaDD[2][2]=4;
      metric.gammaUU[0][0]=metric.gammaUU[1][1]=metric.gammaUU[2][2]=0.25;
      metric.detgamma=64; metric.sqrt_detgamma=8;
      for(int side=0;side<2;++side)
        for(int i=0;i<3;++i) closure[side].P[i][i]=1.0/12;
      state[1].N=DBL_MAX;
    }
    if(scenario==6) {
      /* A rounded supplied current above the product is still within tolerance. */
      velocity[1][0]=0.5; number[1][0]=nextafter(0.5,1.0);
    }
    double N=-1, E=-2, F[3]={-3,-4,-5};
    ghl_error_codes_t error=ghl_m1_compute_neutrino_rusanov_flux(
          params,&nu,&metric,ghl_m1_dirn0,&state[0],&state[1],&closure[0],&closure[1],
          number[0],number[1],velocity[0],velocity[1],speed,&N,&E,F);
    if(scenario==6) {
      if(error!=ghl_success) fail_case("roundoff-sized current discrepancy rejected",scenario);
      check_close(N,0.25,"rounded current flux mismatch",scenario);
    } else if(error!=ghl_error_m1_invalid_state || N!=-1 || E!=-2
              || F[0]!=-3 || F[1]!=-4 || F[2]!=-5) {
      fail_case("Rusanov arithmetic failure was not transactional",scenario);
    }
  }
}

int main(int argc, char **argv) {
  const char *fixture_dir = "Unit_Tests/data/m1_thcm1";
  if(argc == 3 && strcmp(argv[1], "--fixture-dir") == 0) {
    fixture_dir = argv[2];
  }
  else if(argc != 1) {
    fail_case("usage: [--fixture-dir PATH]", -1);
    return 1;
  }
  ghl_m1_parameters m1_params = { 0 };
  if(ghl_m1_initialize(
           1.0e-10, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-12, 20, 1.0e-10, &m1_params)
     != ghl_success) {
    fail_case("M1 initialization failed", -1);
  }
  ghl_m1_neutrino_parameters nu_params;
  m1_neutrino_seeded_default_parameters(&nu_params);
  check_number_product_overflow(&m1_params);
  check_rusanov_arithmetic_boundaries(&m1_params);
  const ghl_m1_neutrino_state *volatile no_state = NULL;
  ghl_m1_closure rejected_closure = { .xi = -1.0 };
  if(ghl_m1_compute_neutrino_closure(&m1_params, NULL, NULL, no_state, &rejected_closure)
           != ghl_error_m1_null_pointer
     || rejected_closure.xi != -1.0) {
    fail_case("NULL neutrino closure state changed output", -1);
  }
  char fixture_error[256] = { 0 };
  if(!m1_thcm1_rusanov_check_neutrino_fixture(
           fixture_dir, &m1_params, &nu_params, fixture_error, sizeof(fixture_error))) {
    fail_case(
          fixture_error[0] != '\0' ? fixture_error
                                   : "neutrino Rusanov fixture evaluation failed",
          -2);
  }
  memset(fixture_error, 0, sizeof(fixture_error));
  if(!m1_thcm1_rusanov_check_current_fixture(
           fixture_dir, &m1_params, &nu_params, fixture_error, sizeof(fixture_error))) {
    fail_case(
          fixture_error[0] != '\0'
                ? fixture_error
                : "current neutrino Rusanov fixture evaluation failed",
          -3);
  }
  else {
    ghl_info(
          "unit_test_m1_neutrino_rusanov_flux: "
          "1024 strict neutrino Rusanov pairs passed\n");
  }

  m1_neutrino_seeded_rng rng
        = { .state = M1_NEUTRINO_SEEDED_PRNG_SEED ^ UINT64_C(0x9e3779b97f4a7c15) };
  for(int case_index = 0; case_index < M1_NEUTRINO_SEEDED_CASE_COUNT; ++case_index) {
    m1_neutrino_seeded_case test_case;
    m1_neutrino_seeded_make_case(&rng, case_index, &test_case);

    ghl_m1_closure closure_L = { 0 };
    ghl_m1_closure closure_R = { 0 };
    if(ghl_m1_compute_neutrino_closure(
             &m1_params, &test_case.metric, &test_case.prims, &test_case.state,
             &closure_L)
             != ghl_success
       || ghl_m1_compute_neutrino_closure(
                &m1_params, &test_case.metric, &test_case.prims,
                &test_case.perturbed_state, &closure_R)
                != ghl_success) {
      fail_case("closure construction failed", case_index);
    }

    double number_flux_L[3] = { 0.0, 0.0, 0.0 };
    double number_flux_R[3] = { 0.0, 0.0, 0.0 };
    double number_velocity_L[3] = { 0.0, 0.0, 0.0 };
    double number_velocity_R[3] = { 0.0, 0.0, 0.0 };
    if(ghl_m1_compute_neutrino_number_flux_from_closure(
             &m1_params, &nu_params, &test_case.metric, &test_case.prims,
             &test_case.state, &closure_L, number_flux_L, number_velocity_L)
             != ghl_success
       || ghl_m1_compute_neutrino_number_flux_from_closure(
                &m1_params, &nu_params, &test_case.metric, &test_case.prims,
                &test_case.perturbed_state, &closure_R, number_flux_R, number_velocity_R)
                != ghl_success) {
      fail_case("number-current construction failed", case_index);
    }

    const ghl_m1_direction_t direction = (ghl_m1_direction_t)(case_index % 3);
    const int d = (int)direction;
    double physical_E_L = NAN, physical_E_R = NAN;
    double physical_F_L[3] = { NAN, NAN, NAN };
    double physical_F_R[3] = { NAN, NAN, NAN };
    const ghl_m1_rad_state rad_L = ghl_m1_neutrino_project_rad_state(&test_case.state);
    const ghl_m1_rad_state rad_R
          = ghl_m1_neutrino_project_rad_state(&test_case.perturbed_state);
    if(ghl_m1_compute_physical_flux(
             &test_case.metric, direction, &rad_L, &closure_L, &physical_E_L,
             physical_F_L)
             != ghl_success
       || ghl_m1_compute_physical_flux(
                &test_case.metric, direction, &rad_R, &closure_R, &physical_E_R,
                physical_F_R)
                != ghl_success) {
      fail_case("physical flux construction failed", case_index);
    }

    const double physical_N_L = test_case.metric.lapse * number_flux_L[d]
                                - test_case.metric.betaU[d] * test_case.state.N;
    const double physical_N_R
          = test_case.metric.lapse * number_flux_R[d]
            - test_case.metric.betaU[d] * test_case.perturbed_state.N;
    const double speed = 0.2 + 0.03 * (double)(case_index % 11);
    double actual_N = NAN, actual_E = NAN;
    double actual_F[3] = { NAN, NAN, NAN };
    if(ghl_m1_compute_neutrino_rusanov_flux(
             &m1_params, &nu_params, &test_case.metric, direction, &test_case.state,
             &test_case.perturbed_state, &closure_L, &closure_R, number_flux_L,
             number_flux_R, number_velocity_L, number_velocity_R, speed, &actual_N,
             &actual_E, actual_F)
       != ghl_success) {
      fail_case("combined neutrino Rusanov rejected generated input", case_index);
    }

    const double state_L[5]
          = { test_case.state.N, test_case.state.E, test_case.state.F[0],
              test_case.state.F[1], test_case.state.F[2] };
    const double state_R[5]
          = { test_case.perturbed_state.N, test_case.perturbed_state.E,
              test_case.perturbed_state.F[0], test_case.perturbed_state.F[1],
              test_case.perturbed_state.F[2] };
    const double physical_L[5] = { physical_N_L, physical_E_L, physical_F_L[0],
                                   physical_F_L[1], physical_F_L[2] };
    const double physical_R[5] = { physical_N_R, physical_E_R, physical_F_R[0],
                                   physical_F_R[1], physical_F_R[2] };
    const double actual[5]
          = { actual_N, actual_E, actual_F[0], actual_F[1], actual_F[2] };
    for(int component = 0; component < 5; ++component) {
      const double expected
            = test_case.metric.sqrt_detgamma
              * (0.5 * (physical_L[component] + physical_R[component])
                 - 0.5 * speed * (state_R[component] - state_L[component]));
      check_close(
            actual[component], expected, "combined Rusanov component mismatch",
            case_index);
    }

    /* The paired state is intentionally changed before this call. */
    if(actual_E == 0.0 || test_case.state.E == test_case.perturbed_state.E) {
      fail_case("baseline/perturbed energy response was not exercised", case_index);
    }

    double rejected_N = 71.0, rejected_E = 72.0;
    double rejected_F[3] = { 73.0, 74.0, 75.0 };
    if(ghl_m1_compute_neutrino_rusanov_flux(
             &m1_params, &nu_params, &test_case.metric, direction, &test_case.state,
             &test_case.perturbed_state, &closure_L, &closure_R, number_flux_L,
             number_flux_R, number_velocity_L, number_velocity_R, -1.0, &rejected_N,
             &rejected_E, rejected_F)
             != ghl_error_m1_invalid_state
       || rejected_N != 71.0 || rejected_E != 72.0 || rejected_F[0] != 73.0
       || rejected_F[1] != 74.0 || rejected_F[2] != 75.0) {
      fail_case("invalid combined Rusanov input was not transactional", case_index);
    }
  }

  ghl_metric_quantities boundary_metric;
  m1_setup_flat_metric(&boundary_metric);
  ghl_primitive_quantities boundary_prims = { 0 };
  boundary_prims.u0 = 1.0;
  check_neutrino_number_transport_boundaries(
        &m1_params, &nu_params, &boundary_metric, &boundary_prims);
  check_neutrino_rusanov_boundaries(
        &m1_params, &nu_params, &boundary_metric, &boundary_prims);

  ghl_info(
        "unit_test_m1_neutrino_rusanov_flux: %d seeded paired cases passed\n",
        M1_NEUTRINO_SEEDED_CASE_COUNT);
  return 0;
}
