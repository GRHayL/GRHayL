#include <math.h>
#include <stdio.h>
#include <string.h>

#include "ghl_flux_source.h"
#include "m1_test_utils.h"
#include "m1_thcm1_rusanov_fixture.h"

/*
 * Boundary and reference coverage for shared and typed M1 Rusanov fluxes.
 * The test checks component-wise arithmetic, physical/scalar wrappers, and
 * transactional validation against the retained generic corpus.
 */

static void fail_case(const char *message, const int case_index) {
  ghl_error("unit_test_rusanov_flux case %d: %s\n", case_index, message);
}

static void setup_primitives(
      ghl_primitive_quantities *restrict prims,
      const double vx,
      const double vy,
      const double vz) {
  *prims = (ghl_primitive_quantities){ 0 };
  prims->rho = 1.0;
  prims->press = 0.4;
  prims->eps = 0.2;
  prims->u0 = 1.0;
  prims->vU[0] = vx;
  prims->vU[1] = vy;
  prims->vU[2] = vz;
  prims->Y_e = 0.2;
  prims->temperature = 1.0;
  prims->entropy = 0.3;
}

static void check_close(
      const double actual,
      const double expected,
      const char *label,
      const int case_index) {
  if(!m1_nearly_equal(actual, expected, 2.0e-12, 2.0e-13)) {
    fail_case(label, case_index);
  }
}

static void check_shared_rusanov_boundaries(void) {
  double state_L[3] = { 1.0, 2.0, 3.0 };
  double state_R[3] = { 1.5, 2.5, 3.5 };
  double flux_L[3] = { 0.1, 0.2, 0.3 };
  double flux_R[3] = { -0.1, -0.2, -0.3 };
  double output[3] = { 151.0, 152.0, 153.0 };
  const double output_before[3] = { 151.0, 152.0, 153.0 };
  if(ghl_calculate_Rusanov_flux(NULL, state_R, flux_L, flux_R, 3, 0.5, output)
           != ghl_error_flux_source_invalid_input
     || ghl_calculate_Rusanov_flux(state_L, NULL, flux_L, flux_R, 3, 0.5, output)
              != ghl_error_flux_source_invalid_input
     || ghl_calculate_Rusanov_flux(state_L, state_R, NULL, flux_R, 3, 0.5, output)
              != ghl_error_flux_source_invalid_input
     || ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, NULL, 3, 0.5, output)
              != ghl_error_flux_source_invalid_input
     || ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, flux_R, 3, 0.5, NULL)
              != ghl_error_flux_source_invalid_input
     || ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, flux_R, 0, 0.5, output)
              != ghl_error_flux_source_invalid_input
     || ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, flux_R, 3, NAN, output)
              != ghl_error_flux_source_invalid_input
     || ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, flux_R, 3, -0.5, output)
              != ghl_error_flux_source_invalid_input) {
    fail_case("shared Rusanov pointer/count/speed validation failed", 400);
  }

  state_L[0] = NAN;
  if(ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, flux_R, 3, 0.5, output)
     != ghl_error_flux_source_invalid_input) {
    fail_case("shared Rusanov nonfinite left state was accepted", 401);
  }
  state_L[0] = 1.0;
  state_R[0] = NAN;
  if(ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, flux_R, 3, 0.5, output)
     != ghl_error_flux_source_invalid_input) {
    fail_case("shared Rusanov nonfinite right state was accepted", 402);
  }
  state_R[0] = 1.5;
  flux_L[0] = NAN;
  if(ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, flux_R, 3, 0.5, output)
     != ghl_error_flux_source_invalid_input) {
    fail_case("shared Rusanov nonfinite left flux was accepted", 403);
  }
  flux_L[0] = 0.1;
  flux_R[0] = NAN;
  if(ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, flux_R, 3, 0.5, output)
     != ghl_error_flux_source_invalid_input) {
    fail_case("shared Rusanov nonfinite right flux was accepted", 404);
  }
  flux_R[0] = -0.1;

  double huge_state_L[1] = { 0.0 };
  double huge_state_R[1] = { DBL_MAX };
  double zero_flux_L[1] = { 0.0 };
  double zero_flux_R[1] = { 0.0 };
  double rejected[1] = { 155.0 };
  if(ghl_calculate_Rusanov_flux(
           huge_state_L, huge_state_R, zero_flux_L, zero_flux_R, 1, DBL_MAX, rejected)
           != ghl_error_flux_source_invalid_input
     || rejected[0] != 155.0) {
    fail_case("shared Rusanov candidate overflow was published", 405);
  }

  if(ghl_calculate_Rusanov_flux(state_L, state_R, flux_L, flux_R, 3, 0.5, output)
     != ghl_success) {
    fail_case("shared Rusanov finite boundary case failed", 406);
  }
  for(int i = 0; i < 3; ++i) {
    const double expected
          = 0.5 * (flux_L[i] + flux_R[i]) - 0.25 * (state_R[i] - state_L[i]);
    check_close(output[i], expected, "shared Rusanov output mismatch", 407);
  }
  const double equal_state_L[1] = { 1.0 };
  const double equal_state_R[1] = { 1.0 };
  const double large_flux_L[1] = { 1.0e308 };
  const double large_flux_R[1] = { 1.0e308 };
  double large_flux_output[1] = { 156.0 };
  if(ghl_calculate_Rusanov_flux(
           equal_state_L, equal_state_R, large_flux_L, large_flux_R, 1, 0.5,
           large_flux_output)
           != ghl_success
     || !isfinite(large_flux_output[0]) || large_flux_output[0] != large_flux_L[0]) {
    fail_case("shared Rusanov rejected representable large flux", 408);
  }
  (void)output_before;
}

static void check_m1_rusanov_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims) {
  const ghl_m1_rad_state state_L = { .E = 1.0, .F = { 0.0, 0.0, 0.0 } };
  const ghl_m1_rad_state state_R = { .E = 1.2, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_closure closure_L = { 0 };
  ghl_m1_closure closure_R = { 0 };
  if(ghl_m1_compute_closure_with_primitives(
           m1_params, metric, prims, &state_L, &closure_L)
           != ghl_success
     || ghl_m1_compute_closure_with_primitives(
              m1_params, metric, prims, &state_R, &closure_R)
              != ghl_success) {
    fail_case("M1 Rusanov boundary closure construction failed", 410);
  }

  double physical_E = 161.0;
  double physical_F[3] = { 162.0, 163.0, 164.0 };
  if(ghl_m1_compute_physical_flux(
           NULL, ghl_m1_dirn0, &state_L, &closure_L, &physical_E, physical_F)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_physical_flux(
              metric, ghl_m1_dirn0, NULL, &closure_L, &physical_E, physical_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_physical_flux(
              metric, ghl_m1_dirn0, &state_L, NULL, &physical_E, physical_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_physical_flux(
              metric, ghl_m1_dirn0, &state_L, &closure_L, NULL, physical_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_physical_flux(
              metric, ghl_m1_dirn0, &state_L, &closure_L, &physical_E, NULL)
              != ghl_error_m1_null_pointer) {
    fail_case("M1 physical flux NULL boundary failed", 411);
  }
  if(ghl_m1_compute_physical_flux(
           metric, (ghl_m1_direction_t)3, &state_L, &closure_L, &physical_E, physical_F)
     != ghl_error_m1_invalid_state) {
    fail_case("M1 physical flux invalid direction was accepted", 412);
  }
  ghl_metric_quantities bad_metric = *metric;
  bad_metric.gammaDD[0][0] = -1.0;
  if(ghl_m1_compute_physical_flux(
           &bad_metric, ghl_m1_dirn0, &state_L, &closure_L, &physical_E, physical_F)
     != ghl_error_m1_invalid_metric) {
    fail_case("M1 physical flux invalid metric was accepted", 413);
  }
  ghl_m1_closure bad_closure = closure_L;
  bad_closure.P[0][0] = NAN;
  if(ghl_m1_compute_physical_flux(
           metric, ghl_m1_dirn0, &state_L, &bad_closure, &physical_E, physical_F)
     != ghl_error_m1_invalid_state) {
    fail_case("M1 physical flux invalid closure was accepted", 414);
  }

  const ghl_m1_rad_state high_state = { .E = DBL_MAX, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_closure high_closure = { 0 };
  if(ghl_m1_compute_closure_with_primitives(
           m1_params, metric, prims, &high_state, &high_closure)
     != ghl_success) {
    fail_case("M1 high-energy closure construction failed", 415);
  }
  ghl_metric_quantities huge_lapse_metric = *metric;
  huge_lapse_metric.lapse = DBL_MAX;
  huge_lapse_metric.betaU[0] = DBL_MAX;
  if(ghl_m1_compute_physical_flux(
           &huge_lapse_metric, ghl_m1_dirn0, &high_state, &high_closure, &physical_E,
           physical_F)
     != ghl_error_m1_invalid_state) {
    fail_case("M1 physical energy overflow was accepted", 416);
  }
  huge_lapse_metric.betaU[0] = 0.0;
  if(ghl_m1_compute_physical_flux(
           &huge_lapse_metric, ghl_m1_dirn0, &high_state, &high_closure, &physical_E,
           physical_F)
     != ghl_error_m1_invalid_state) {
    fail_case("M1 physical momentum overflow was accepted", 417);
  }

  const double flux_F_L[3] = { 0.1, 0.2, 0.3 };
  const double flux_F_R[3] = { -0.1, -0.2, -0.3 };
  double rusanov_E = 171.0;
  double rusanov_F[3] = { 172.0, 173.0, 174.0 };
  if(ghl_m1_compute_rusanov_flux(
           NULL, &state_R, 0.1, flux_F_L, 0.2, flux_F_R, 0.5, &rusanov_E, rusanov_F)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_rusanov_flux(
              &state_L, NULL, 0.1, flux_F_L, 0.2, flux_F_R, 0.5, &rusanov_E, rusanov_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_rusanov_flux(
              &state_L, &state_R, 0.1, NULL, 0.2, flux_F_R, 0.5, &rusanov_E, rusanov_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_rusanov_flux(
              &state_L, &state_R, 0.1, flux_F_L, 0.2, NULL, 0.5, &rusanov_E, rusanov_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_rusanov_flux(
              &state_L, &state_R, 0.1, flux_F_L, 0.2, flux_F_R, 0.5, NULL, rusanov_F)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_rusanov_flux(
              &state_L, &state_R, 0.1, flux_F_L, 0.2, flux_F_R, 0.5, &rusanov_E, NULL)
              != ghl_error_m1_null_pointer) {
    fail_case("typed M1 Rusanov NULL boundary failed", 418);
  }
  if(ghl_m1_compute_rusanov_flux(
           &state_L, &state_R, 0.1, flux_F_L, 0.2, flux_F_R, NAN, &rusanov_E, rusanov_F)
           != ghl_error_m1_invalid_state
     || ghl_m1_compute_rusanov_flux(
              &state_L, &state_R, 0.1, flux_F_L, 0.2, flux_F_R, -1.0, &rusanov_E,
              rusanov_F)
              != ghl_error_m1_invalid_state
     || ghl_m1_compute_rusanov_flux(
              &state_L, &state_R, NAN, flux_F_L, 0.2, flux_F_R, 0.5, &rusanov_E,
              rusanov_F)
              != ghl_error_m1_invalid_state
     || ghl_m1_compute_rusanov_flux(
              &state_L, &state_R, 0.1, flux_F_L, NAN, flux_F_R, 0.5, &rusanov_E,
              rusanov_F)
              != ghl_error_m1_invalid_state) {
    fail_case("typed M1 Rusanov scalar validation failed", 419);
  }
  ghl_m1_rad_state bad_state = state_L;
  bad_state.E = NAN;
  if(ghl_m1_compute_rusanov_flux(
           &bad_state, &state_R, 0.1, flux_F_L, 0.2, flux_F_R, 0.5, &rusanov_E,
           rusanov_F)
     != ghl_error_m1_invalid_state) {
    fail_case("typed M1 Rusanov left state validation failed", 420);
  }
  bad_state = state_R;
  bad_state.E = NAN;
  if(ghl_m1_compute_rusanov_flux(
           &state_L, &bad_state, 0.1, flux_F_L, 0.2, flux_F_R, 0.5, &rusanov_E,
           rusanov_F)
     != ghl_error_m1_invalid_state) {
    fail_case("typed M1 Rusanov right state validation failed", 421);
  }
  double bad_flux_F[3] = { NAN, 0.0, 0.0 };
  const double rusanov_before[4] = { 171.0, 172.0, 173.0, 174.0 };
  rusanov_E = rusanov_before[0];
  memcpy(rusanov_F, rusanov_before + 1, sizeof(rusanov_F));
  if(ghl_m1_compute_rusanov_flux(
           &state_L, &state_R, 0.1, bad_flux_F, 0.2, flux_F_R, 0.5, &rusanov_E,
           rusanov_F)
           != ghl_error_m1_invalid_state
     || rusanov_E != rusanov_before[0]
     || memcmp(rusanov_F, rusanov_before + 1, sizeof(rusanov_F)) != 0) {
    fail_case("typed M1 Rusanov helper failure was not transactional", 422);
  }
  ghl_m1_rad_state huge_state_L = { .E = 0.0, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_rad_state huge_state_R = { .E = DBL_MAX, .F = { 0.0, 0.0, 0.0 } };
  rusanov_E = rusanov_before[0];
  memcpy(rusanov_F, rusanov_before + 1, sizeof(rusanov_F));
  if(ghl_m1_compute_rusanov_flux(
           &huge_state_L, &huge_state_R, 0.0, flux_F_L, 0.0, flux_F_R, DBL_MAX,
           &rusanov_E, rusanov_F)
           != ghl_error_m1_invalid_state
     || rusanov_E != rusanov_before[0]
     || memcmp(rusanov_F, rusanov_before + 1, sizeof(rusanov_F)) != 0) {
    fail_case("typed M1 Rusanov candidate overflow was published", 423);
  }

  double number_flux = 181.0;
  if(ghl_m1_compute_number_rusanov_flux(1.0, 1.2, 0.1, 0.2, 0.5, NULL)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_number_rusanov_flux(NAN, 1.2, 0.1, 0.2, 0.5, &number_flux)
              != ghl_error_m1_invalid_state
     || ghl_m1_compute_number_rusanov_flux(1.0, 1.2, NAN, 0.2, 0.5, &number_flux)
              != ghl_error_m1_invalid_state
     || ghl_m1_compute_number_rusanov_flux(1.0, 1.2, 0.1, NAN, 0.5, &number_flux)
              != ghl_error_m1_invalid_state
     || ghl_m1_compute_number_rusanov_flux(1.0, 1.2, 0.1, 0.2, NAN, &number_flux)
              != ghl_error_m1_invalid_state
     || ghl_m1_compute_number_rusanov_flux(1.0, 1.2, 0.1, 0.2, -1.0, &number_flux)
              != ghl_error_m1_invalid_state) {
    fail_case("scalar M1 Rusanov validation failed", 424);
  }
  number_flux = 181.0;
  if(ghl_m1_compute_number_rusanov_flux(0.0, DBL_MAX, 0.0, 0.0, DBL_MAX, &number_flux)
           != ghl_error_m1_invalid_state
     || number_flux != 181.0) {
    fail_case("scalar M1 Rusanov candidate overflow was published", 425);
  }
}

static void check_physical_flux(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_direction_t direction,
      const ghl_m1_rad_state *restrict state,
      const ghl_m1_closure *restrict closure,
      const int case_index) {

  double actual_E = NAN;
  double actual_F[3] = { NAN, NAN, NAN };
  if(ghl_m1_compute_physical_flux(metric, direction, state, closure, &actual_E, actual_F)
     != ghl_success) {
    fail_case("physical flux rejected an admissible state", case_index);
  }

  double raised_F[3] = { 0.0, 0.0, 0.0 };
  ghl_raise_lower_vector_3D(metric->gammaUU, state->F, raised_F);
  const int d = (int)direction;
  const double expected_E = metric->lapse * raised_F[d] - metric->betaU[d] * state->E;
  check_close(actual_E, expected_E, "physical energy flux mismatch", case_index);

  for(int i = 0; i < 3; ++i) {
    double P_mixed = 0.0;
    for(int k = 0; k < 3; ++k) {
      P_mixed += closure->P[d][k] * metric->gammaDD[k][i];
    }
    const double expected_F = metric->lapse * P_mixed - metric->betaU[d] * state->F[i];
    check_close(actual_F[i], expected_F, "physical momentum flux mismatch", case_index);
  }
}

static void check_case(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict state_L,
      const ghl_m1_rad_state *restrict state_R,
      const int case_index) {

  ghl_m1_closure closure_L = { 0 };
  ghl_m1_closure closure_R = { 0 };
  if(ghl_m1_compute_closure_with_primitives(
           m1_params, metric, prims, state_L, &closure_L)
           != ghl_success
     || ghl_m1_compute_closure_with_primitives(
              m1_params, metric, prims, state_R, &closure_R)
              != ghl_success) {
    fail_case("closure construction failed", case_index);
  }

  for(int direction = 0; direction < 3; ++direction) {
    check_physical_flux(
          metric, (ghl_m1_direction_t)direction, state_L, &closure_L, case_index);
    check_physical_flux(
          metric, (ghl_m1_direction_t)direction, state_R, &closure_R, case_index);
  }

  const double speed = 0.37 + 0.11 * (double)case_index;
  const double flux_E_L = 0.31 - 0.07 * (double)case_index;
  const double flux_E_R = -0.23 + 0.05 * (double)case_index;
  const double flux_F_L[3] = { 0.17, -0.08, 0.12 };
  const double flux_F_R[3] = { -0.04, 0.21, -0.09 };
  double actual_E = NAN;
  double actual_F[3] = { NAN, NAN, NAN };
  if(ghl_m1_compute_rusanov_flux(
           state_L, state_R, flux_E_L, flux_F_L, flux_E_R, flux_F_R, speed, &actual_E,
           actual_F)
     != ghl_success) {
    fail_case("typed Rusanov flux rejected finite input", case_index);
  }

  const double expected_E
        = 0.5 * (flux_E_L + flux_E_R) - 0.5 * speed * (state_R->E - state_L->E);
  check_close(actual_E, expected_E, "typed Rusanov energy mismatch", case_index);
  for(int i = 0; i < 3; ++i) {
    const double expected_F = 0.5 * (flux_F_L[i] + flux_F_R[i])
                              - 0.5 * speed * (state_R->F[i] - state_L->F[i]);
    check_close(actual_F[i], expected_F, "typed Rusanov momentum mismatch", case_index);
  }

  const double N_L = 0.42 + 0.1 * (double)case_index;
  const double N_R = 0.77 - 0.04 * (double)case_index;
  const double nflux_L = -0.19 + 0.02 * (double)case_index;
  const double nflux_R = 0.28 - 0.03 * (double)case_index;
  double actual_N = NAN;
  if(ghl_m1_compute_number_rusanov_flux(N_L, N_R, nflux_L, nflux_R, speed, &actual_N)
     != ghl_success) {
    fail_case("scalar Rusanov flux rejected finite input", case_index);
  }
  const double expected_N = 0.5 * (nflux_L + nflux_R) - 0.5 * speed * (N_R - N_L);
  check_close(actual_N, expected_N, "scalar Rusanov mismatch", case_index);

  /* A late invalid-input failure must not publish a candidate. */
  double rejected_E = 91.0;
  double rejected_F[3] = { 92.0, 93.0, 94.0 };
  if(ghl_m1_compute_rusanov_flux(
           state_L, state_R, flux_E_L, flux_F_L, flux_E_R, flux_F_R, -1.0, &rejected_E,
           rejected_F)
           != ghl_error_m1_invalid_state
     || rejected_E != 91.0 || rejected_F[0] != 92.0 || rejected_F[1] != 93.0
     || rejected_F[2] != 94.0) {
    fail_case("invalid typed Rusanov input was not transactional", case_index);
  }

  double rejected_N = 95.0;
  if(ghl_m1_compute_number_rusanov_flux(N_L, N_R, nflux_L, nflux_R, -1.0, &rejected_N)
           != ghl_error_m1_invalid_state
     || rejected_N != 95.0) {
    fail_case("invalid scalar Rusanov input was not transactional", case_index);
  }
}

static void check_high_energy_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims) {

  const double energies[] = { 1.0e155, 1.0e308 };
  const double flux_factors[] = { 0.0, 0.1 };
  int case_index = 0;
  for(size_t energy_index = 0; energy_index < sizeof(energies) / sizeof(energies[0]);
      ++energy_index) {
    for(size_t flux_index = 0;
        flux_index < sizeof(flux_factors) / sizeof(flux_factors[0]);
        ++flux_index, ++case_index) {
      const double energy = energies[energy_index];
      const ghl_m1_rad_state state
            = { .E = energy, .F = { flux_factors[flux_index] * energy, 0.0, 0.0 } };
      ghl_m1_closure closure = { 0 };
      if(ghl_m1_compute_closure_with_primitives(
               m1_params, metric, prims, &state, &closure)
         != ghl_success) {
        fail_case("finite high-energy closure was rejected", 100 + case_index);
      }
      if(!isfinite(closure.chi) || !isfinite(closure.xi)
         || !isfinite(closure.root_residual)
         || closure.solve_status == ghl_m1_closure_solve_invalid) {
        fail_case(
              "finite high-energy closure published invalid metadata", 100 + case_index);
      }
      for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
          if(!isfinite(closure.P[i][j] / energy)) {
            fail_case(
                  "finite high-energy closure published a nonfinite tensor",
                  100 + case_index);
          }
        }
      }
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
  char fixture_error[256] = { 0 };
  if(!m1_thcm1_rusanov_check_generic_fixture(
           fixture_dir, fixture_error, sizeof(fixture_error))) {
    fail_case(
          fixture_error[0] != '\0' ? fixture_error
                                   : "generic Rusanov fixture evaluation failed",
          -2);
  }
  else {
    ghl_info("unit_test_rusanov_flux: 1024 strict generic Rusanov pairs passed\n");
  }

  ghl_metric_quantities flat_metric;
  m1_setup_flat_metric(&flat_metric);
  ghl_primitive_quantities flat_prims;
  setup_primitives(&flat_prims, 0.0, 0.0, 0.0);
  check_high_energy_closure(&m1_params, &flat_metric, &flat_prims);
  const ghl_m1_rad_state flat_L = { .E = 1.0, .F = { 0.19, -0.11, 0.07 } };
  const ghl_m1_rad_state flat_R = { .E = 1.4, .F = { -0.08, 0.15, 0.12 } };
  check_case(&m1_params, &flat_metric, &flat_prims, &flat_L, &flat_R, 0);

  ghl_metric_quantities curved_metric;
  m1_setup_metric_anchor_B(&curved_metric);
  ghl_primitive_quantities curved_prims;
  setup_primitives(&curved_prims, 0.06, -0.04, 0.03);
  check_high_energy_closure(&m1_params, &curved_metric, &curved_prims);
  const ghl_m1_rad_state curved_L = { .E = 0.83, .F = { 0.12, -0.09, 0.05 } };
  const ghl_m1_rad_state curved_R = { .E = 1.17, .F = { -0.11, 0.14, 0.08 } };
  check_case(&m1_params, &curved_metric, &curved_prims, &curved_L, &curved_R, 1);

  check_shared_rusanov_boundaries();
  check_m1_rusanov_boundaries(&m1_params, &flat_metric, &flat_prims);

  if(ghl_m1_compute_rusanov_flux(
           NULL, &flat_L, 0.0, (double[3]){ 0.0, 0.0, 0.0 }, 0.0,
           (double[3]){ 0.0, 0.0, 0.0 }, 1.0, NULL, NULL)
     != ghl_error_m1_null_pointer) {
    fail_case("NULL typed Rusanov input was not rejected", -1);
  }

  ghl_info("unit_test_rusanov_flux: typed physical/E-F/scalar Rusanov cases passed\n");
  return 0;
}
