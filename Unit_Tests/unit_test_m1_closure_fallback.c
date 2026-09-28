#define _POSIX_C_SOURCE 200809L
#include "../GRHayL/Radiation/ghl_m1_closure_private.h"
#include "../GRHayL/Radiation/ghl_m1_utils.h"
#include "ghl_m1.h"
#include "m1_test_utils.h"

#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Coverage for the Eulerian Minerbo admissibility fallback in
 * GRHayL/Radiation/ghl_m1_closure.c.
 *
 * Two distinct states publish four_point_compatibility == false: an exact-zero
 * Eulerian flux whose covariant thin dyad vanishes identically, and a
 * finite-flux candidate rejected by the physical PSD check. Both are separately
 * counted by ghl_m1_closure_counters. This test owns three properties:
 *
 *   1. Every published pressure tensor is exactly symmetric, on the fallback
 *      path as well as the primary path. The fallback builds its tensor from a
 *      direction dyad and the inverse metric; without explicit symmetrization
 *      those terms differ in their last bit under an index swap, and the shared
 *      tensor validator then rejects the tensor whenever an off-diagonal
 *      component is small compared with E.
 *   2. The fallback publishes rather than failing. A directional sweep through
 *      the PSD regime returns ghl_success for every state.
 *   3. The PSD regime stays where it is documented. Flux parallel or
 *      antiparallel to the fluid velocity does not reach it; transverse flux
 *      above roughly half light speed does. Pinning both directions makes a
 *      future closure change that widens the regime a test failure rather than
 *      a silent behavior change.
 */

static int failures = 0;

static void fail(const char *const what) {
  fprintf(stderr, "unit_test_m1_closure_fallback: %s\n", what);
  failures++;
}

static void setup_parameters(ghl_m1_parameters *restrict params) {
  if(ghl_m1_initialize(1.0e-8, 1.0e-12, 1.0e-12, 1.0e-6, 1.0e-8, 20, 1.0e-5, params)
     != ghl_success) {
    fail("initializer rejected the campaign control values");
    exit(1);
  }
}

typedef struct {
  int calls;
  int fail_call;
  ghl_error_codes_t failure;
  enum {
    root_script_default,
    root_script_same_sign,
    root_script_short_interpolation,
    root_script_right_guard_false
  } mode;
  double xi[8];
  double same_sign_normalized_residual[2];
} scripted_root_context;

static ghl_error_codes_t scripted_root_evaluate(
      void *const opaque,
      const double xi,
      ghl_m1_closure_evaluation *const evaluation) {
  scripted_root_context *const context = opaque;
  const int call = context->calls++;
  context->xi[call] = xi;
  if(call == context->fail_call) {
    return context->failure;
  }
  memset(evaluation, 0, sizeof(*evaluation));
  evaluation->chi = 10.0 + call;
  evaluation->physical_xi = 20.0 + call;
  evaluation->P[0][0] = 30.0 + call;
  evaluation->normalized_residual = 0.5;
  evaluation->residual = call == 0 ? -1.0 : 1.0;
  if(context->mode == root_script_same_sign) {
    evaluation->residual = 2.0 - call;
    evaluation->normalized_residual = context->same_sign_normalized_residual[call];
  }
  else if(context->mode == root_script_short_interpolation) {
    if(call == 0) {
      evaluation->residual = -1.0;
      evaluation->normalized_residual = 1.0;
    }
    else if(call == 1) {
      evaluation->residual = 1.0e-7;
      evaluation->normalized_residual = 1.0e-7;
    }
    else if(call == 2) {
      evaluation->residual = 0.5e-7;
      evaluation->normalized_residual = 0.5e-7;
    }
    else if(call == 3) {
      evaluation->residual = 0.25e-7;
      evaluation->normalized_residual = 0.25e-7;
    }
    else {
      evaluation->residual = -0.5;
      evaluation->normalized_residual = 0.5;
    }
  }
  else if(context->mode == root_script_right_guard_false) {
    if(call == 0) {
      evaluation->residual = -1.0;
      evaluation->normalized_residual = 1.0;
    }
    else if(call == 1) {
      evaluation->residual = 0.1;
      evaluation->normalized_residual = 0.1;
    }
    else if(call == 2) {
      evaluation->residual = 0.5;
      evaluation->normalized_residual = 0.5;
    }
    else {
      evaluation->residual = -0.2;
      evaluation->normalized_residual = 0.2;
    }
  }
  return ghl_success;
}

static int root_result_is_sentinel(const ghl_m1_closure_root_result *const root) {
  return root->xi == -7.0 && root->iterations == -8
         && root->status == ghl_m1_closure_solve_iteration_exhausted
         && root->evaluation.P[0][0] == -9.0 && root->evaluation.chi == -10.0
         && root->evaluation.normalized_residual == -11.0;
}

static ghl_m1_closure_root_result root_result_sentinel(void) {
  ghl_m1_closure_root_result root = { 0 };
  root.xi = -7.0;
  root.iterations = -8;
  root.status = ghl_m1_closure_solve_iteration_exhausted;
  root.evaluation.P[0][0] = -9.0;
  root.evaluation.chi = -10.0;
  root.evaluation.normalized_residual = -11.0;
  return root;
}

static void check_invalid_state_counter(const unsigned long long expected) {
  ghl_m1_closure_counters counters;
  ghl_m1_get_closure_counters(&counters);
  if(counters.invalid_state != expected) {
    fail("private root driver counted an evaluation failure incorrectly");
  }
}

static void check_private_root_driver(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  const int fail_calls[] = { 0, 1, 2 };
  const int expected_calls[] = { 1, 2, 3 };
  for(size_t i = 0; i < sizeof(fail_calls) / sizeof(fail_calls[0]); ++i) {
    scripted_root_context context
          = { .fail_call = fail_calls[i], .failure = ghl_error_m1_invalid_state };
    ghl_m1_closure_root_result root = root_result_sentinel();
    ghl_m1_reset_closure_counters();
    const ghl_error_codes_t error = ghl_m1_closure_private_solve_root(
          &params, scripted_root_evaluate, &context, &root);
    if(error != ghl_error_m1_invalid_state || context.calls != expected_calls[i]
       || !root_result_is_sentinel(&root)) {
      fail("private root driver did not propagate an endpoint/interior error "
           "atomically");
    }
    check_invalid_state_counter(1ULL);
  }

  scripted_root_context same_sign = { .fail_call = -1,
                                      .mode = root_script_same_sign,
                                      .same_sign_normalized_residual = { 0.2, 0.1 } };
  ghl_m1_closure_root_result root = { 0 };
  if(ghl_m1_closure_private_solve_root(
           &params, scripted_root_evaluate, &same_sign, &root)
           != ghl_success
     || same_sign.calls != 2 || root.status != ghl_m1_closure_solve_endpoint_fallback
     || root.xi != 1.0 || root.evaluation.P[0][0] != 31.0 || root.evaluation.chi != 11.0
     || root.evaluation.normalized_residual != 0.1) {
    fail("same-sign root endpoints did not select the lower-residual endpoint payload");
    return;
  }
  if(ghl_m1_closure_private_check_residual_gate(root.evaluation.normalized_residual, 0.1)
           != ghl_success
     || ghl_m1_closure_private_check_residual_gate(
              root.evaluation.normalized_residual, 0.05)
              != ghl_error_m1_closure_residual_too_large) {
    fail("selected endpoint did not retain the existing residual publication gate");
  }
  ghl_m1_closure_failure_stage_t stage;
  ghl_m1_get_last_closure_failure_stage(&stage);
  if(stage != ghl_m1_closure_failure_residual_gate) {
    fail("endpoint residual rejection did not retain its diagnostic stage");
  }

  scripted_root_context same_sign_first
        = { .fail_call = -1,
            .mode = root_script_same_sign,
            .same_sign_normalized_residual = { 0.1, 0.2 } };
  root = root_result_sentinel();
  if(ghl_m1_closure_private_solve_root(
           &params, scripted_root_evaluate, &same_sign_first, &root)
           != ghl_success
     || same_sign_first.calls != 2
     || root.status != ghl_m1_closure_solve_endpoint_fallback || root.xi != 0.0
     || root.iterations != 0 || root.evaluation.P[0][0] != 30.0
     || root.evaluation.chi != 10.0 || root.evaluation.residual != 2.0
     || root.evaluation.normalized_residual != 0.1) {
    fail("same-sign root endpoints did not select the first endpoint payload");
    return;
  }

  params.closure_root_tolerance = 1.0e-5;
  params.closure_root_max_iterations = 3;
  scripted_root_context short_step
        = { .fail_call = -1, .mode = root_script_short_interpolation };
  if(ghl_m1_closure_private_solve_root(
           &params, scripted_root_evaluate, &short_step, &root)
           != ghl_success
     || short_step.calls != 5 || root.status != ghl_m1_closure_solve_iteration_exhausted
     || root.iterations != 3 || root.xi != short_step.xi[4]) {
    fail("short interpolation case did not finish the scripted root schedule");
    return;
  }
  const double first_tolerance = 2.0 * DBL_EPSILON + 0.5 * params.closure_root_tolerance;
  const double raw_secant_step = 1.0e-7 / (1.0 + 1.0e-7);
  if(!(raw_secant_step < 0.5 * params.closure_root_tolerance)
     || fabs((1.0 - short_step.xi[2]) - first_tolerance) > 8.0 * DBL_EPSILON
     || fabs((short_step.xi[2] - short_step.xi[3]) - first_tolerance) > 8.0 * DBL_EPSILON
     || fabs(short_step.xi[4] - 0.5 * short_step.xi[3]) > DBL_EPSILON) {
    fail("sub-tolerance interpolation did not fall back to bracket bisection");
  }

  params.closure_root_max_iterations = 2;
  scripted_root_context right_guard
        = { .fail_call = -1, .mode = root_script_right_guard_false };
  if(ghl_m1_closure_private_solve_root(
           &params, scripted_root_evaluate, &right_guard, &root)
           != ghl_success
     || right_guard.calls != 4
     || fabs(right_guard.xi[3] - 0.5 * right_guard.xi[2]) > DBL_EPSILON) {
    fail("Brent interpolation did not bisect when the residual ordering rejected it");
  }
}

static void record_validation_reason_between_snapshot_load_and_cas(void *const opaque) {
  int *const calls = opaque;
  ++*calls;
  ghl_m1_record_closure_validation_failure(ghl_m1_closure_validation_symmetry);
}

static void check_forced_closure_snapshot_retry(void) {
  int interleave_calls = 0;
  ghl_m1_reset_closure_counters();
  ghl_m1_get_closure_counters(NULL);
  ghl_m1_get_last_closure_failure_stage(NULL);
  ghl_m1_get_last_closure_validation_reason(NULL);
  ghl_m1_closure_private_record_stage_with_interleave(
        ghl_m1_closure_failure_tensor_validation,
        record_validation_reason_between_snapshot_load_and_cas, &interleave_calls);
  ghl_m1_closure_failure_stage_t stage;
  int reason;
  ghl_m1_get_last_closure_failure_stage(&stage);
  ghl_m1_get_last_closure_validation_reason(&reason);
  if(interleave_calls != 1 || stage != ghl_m1_closure_failure_tensor_validation
     || reason != ghl_m1_closure_validation_symmetry) {
    fail("stale closure snapshot retry lost the competing stage or reason field");
  }
}

static void check_private_comoving_invariant(void) {
  double H2_clipped = -4.0, residual = -5.0, normalized_residual = -6.0;
  double physical_xi = -7.0;
  ghl_m1_reset_closure_counters();
  if(ghl_m1_closure_private_evaluate_invariant(
           2.0, 1.0, 4.0, 0.5, &H2_clipped, &residual, &normalized_residual,
           &physical_xi)
           != ghl_success
     || H2_clipped != 1.0 || residual != 0.0 || normalized_residual != 0.0
     || physical_xi != 0.5) {
    fail("private invariant evaluation changed the normalized comoving quantities");
  }

  if(ghl_m1_closure_private_evaluate_invariant(
           NAN, 0.0, 1.0, 0.0, &H2_clipped, &residual, &normalized_residual,
           &physical_xi)
     != ghl_error_m1_invalid_state) {
    fail("private invariant accepted a nonfinite comoving energy");
  }
  ghl_m1_closure_failure_stage_t stage;
  ghl_m1_get_last_closure_failure_stage(&stage);
  if(stage != ghl_m1_closure_failure_comoving_energy) {
    fail("nonfinite comoving energy recorded the wrong failure stage");
  }
  if(ghl_m1_closure_private_evaluate_invariant(
           0.0, 0.0, 1.0, 0.0, &H2_clipped, &residual, &normalized_residual,
           &physical_xi)
     != ghl_error_m1_invalid_state) {
    fail("private invariant accepted zero comoving energy");
  }

  const double tolerance = 1024.0 * DBL_EPSILON;
  if(ghl_m1_closure_private_evaluate_invariant(
           1.0, -tolerance, 1.0, 0.0, &H2_clipped, &residual, &normalized_residual,
           &physical_xi)
           != ghl_success
     || H2_clipped != 0.0 || physical_xi != 0.0) {
    fail("the exact tolerated negative comoving contraction was not clipped to zero");
  }
  if(ghl_m1_closure_private_evaluate_invariant(
           1.0, nextafter(-tolerance, -INFINITY), 1.0, 0.0, &H2_clipped, &residual,
           &normalized_residual, &physical_xi)
     != ghl_error_m1_invalid_state) {
    fail("private invariant accepted a contraction below its tolerance");
  }
  ghl_m1_get_last_closure_failure_stage(&stage);
  if(stage != ghl_m1_closure_failure_comoving_flux_norm) {
    fail("negative comoving contraction recorded the wrong failure stage");
  }
  if(ghl_m1_closure_private_evaluate_invariant(
           1.0, NAN, 1.0, 0.0, &H2_clipped, &residual, &normalized_residual,
           &physical_xi)
     != ghl_error_m1_invalid_state) {
    fail("private invariant accepted a nonfinite comoving contraction");
  }

  if(ghl_m1_closure_private_evaluate_invariant(
           1.0, 0.0, 0.0, 0.0, &H2_clipped, &residual, &normalized_residual,
           &physical_xi)
           != ghl_error_m1_invalid_state
     || ghl_m1_closure_private_evaluate_invariant(
              1.0, 0.0, INFINITY, 0.0, &H2_clipped, &residual, &normalized_residual,
              &physical_xi)
              != ghl_error_m1_invalid_state) {
    fail("underflowed or overflowed residual scale was accepted");
  }
  const double smallest_positive = nextafter(0.0, 1.0);
  H2_clipped = 11.0;
  residual = 12.0;
  normalized_residual = 13.0;
  physical_xi = 14.0;
  if(ghl_m1_closure_private_evaluate_invariant(
           smallest_positive, 1.0, 1.0, 0.0, &H2_clipped, &residual,
           &normalized_residual, &physical_xi)
           != ghl_error_m1_invalid_state
     || H2_clipped != 11.0 || residual != 12.0 || normalized_residual != 13.0
     || physical_xi != 14.0) {
    fail("tiny-J division did not reject atomically after finite normalization");
  }

  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  const double Pdd[4][4] = { { 0.0 }, { 0.0, 2.0, 0.0, 0.0 }, { 0.0 }, { 0.0 } };
  ghl_m1_closure_evaluation evaluation = { 0 };
  if(ghl_m1_closure_private_finish_evaluation(
           &metric, Pdd, 2.0, 1.0, 2.0, 1.0, 4.0, 0.5, 0.5, &evaluation)
           != ghl_success
     || evaluation.P[0][0] != 2.0 || evaluation.P[1][1] != 0.0
     || evaluation.P[2][2] != 0.0 || evaluation.chi != 0.5 || evaluation.residual != 0.0
     || evaluation.normalized_residual != 0.0 || evaluation.physical_xi != 0.5) {
    fail("production evaluation finalizer did not preserve its invariant payload");
  }
  evaluation = (ghl_m1_closure_evaluation){
    .P = { { 17.0, 18.0, 19.0 }, { 20.0, 21.0, 22.0 }, { 23.0, 24.0, 25.0 } },
    .chi = 26.0,
    .residual = 27.0,
    .normalized_residual = 28.0,
    .physical_xi = 29.0
  };
  if(ghl_m1_closure_private_finish_evaluation(
           &metric, Pdd, 2.0, 1.0, 2.0, 1.0, 0.0, 0.5, 0.5, &evaluation)
           != ghl_error_m1_invalid_state
     || evaluation.P[0][0] != 17.0 || evaluation.P[1][1] != 21.0
     || evaluation.chi != 26.0 || evaluation.residual != 27.0
     || evaluation.normalized_residual != 28.0 || evaluation.physical_xi != 29.0) {
    fail("failed production evaluation finalizer changed its output payload");
  }
}

static void check_private_fallback_invariants(void) {
  ghl_m1_closure candidate
        = { .P = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }, { 7.0, 8.0, 9.0 } },
            .chi = 0.25,
            .xi = 0.75,
            .root_residual = 0.125,
            .root_iterations = 4,
            .solve_status = ghl_m1_closure_solve_endpoint_fallback,
            .four_point_compatibility = false };
  const ghl_m1_closure sentinel = { .xi = -17.0, .chi = -18.0, .root_iterations = -19 };
  ghl_m1_closure closure = sentinel;
  ghl_m1_comoving comoving
        = { .J = 2.0, .HD = { 1.0, 0.0, 0.0 }, .HU = { 1.0, 0.0, 0.0 } };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_error_m1_invalid_state, NULL, NULL, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("failed comoving conversion changed the fallback output");
  }
  comoving.J = 0.0;
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("fallback finalizer accepted nonpositive comoving energy");
  }
  comoving.J = -1.0;
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("fallback finalizer accepted negative comoving energy");
  }

  comoving = (ghl_m1_comoving){ .J = 2.0,
                                .HD = { 1.0, 0.0, 0.0 },
                                .HU = { 1.0, 0.0, 0.0 } };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_success
     || closure.xi != 0.5 || closure.P[0][0] != candidate.P[0][0]
     || closure.chi != candidate.chi
     || closure.root_iterations != candidate.root_iterations) {
    fail("ordinary fallback finalizer did not publish its validated flux factor");
  }

  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = 1.0 };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_success
     || closure.xi != 0.0) {
    fail("ordinary fallback did not preserve a zero comoving contraction");
  }
  closure = sentinel;
  const double tolerated_xi = 1.0 + 1024.0 * DBL_EPSILON;
  comoving = (ghl_m1_comoving){ .J = 1.0,
                                .HD = { tolerated_xi * tolerated_xi, 0.0, 0.0 },
                                .HU = { 1.0, 0.0, 0.0 } };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_success
     || closure.xi != 1.0) {
    fail("ordinary fallback did not clip the exact tolerated streaming overshoot");
  }

  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = 1.0,
                                .HD = { 1.21, 0.0, 0.0 },
                                .HU = { 1.0, 0.0, 0.0 } };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("ordinary fallback accepted a flux contraction above J squared");
  }
  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = 1.0, .Hn = 1.0e-8 };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_success
     || closure.xi != 0.0) {
    fail("ordinary fallback failed to clip a tolerated negative contraction");
  }
  closure = sentinel;
  comoving.Hn = 1.0e-6;
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("ordinary fallback accepted a negative contraction beyond its tolerance");
  }
  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = 1.0, .HD = { INFINITY, 0.0, 0.0 }, .HU = { 1.0 } };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("ordinary fallback accepted a nonfinite contraction");
  }

  const double scaled_J = DBL_MIN * 0.5;
  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = scaled_J,
                                .HD = { scaled_J * 0.5, 0.0, 0.0 },
                                .HU = { scaled_J * 0.5, 0.0, 0.0 } };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_success
     || closure.xi != 0.5) {
    fail("scaled fallback finalizer did not preserve the normalized flux factor");
  }
  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = scaled_J, .Hn = scaled_J * 1.0e-8 };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_success
     || closure.xi != 0.0) {
    fail("scaled fallback failed to clip a tolerated negative contraction");
  }
  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = scaled_J, .Hn = scaled_J * 2.0 };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("scaled fallback accepted a negative contraction beyond its tolerance");
  }
  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = scaled_J,
                                .HD = { scaled_J * 1.1, 0.0, 0.0 },
                                .HU = { scaled_J * 1.1, 0.0, 0.0 } };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("scaled fallback accepted a flux contraction outside the tolerated cone");
  }
  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = scaled_J, .Hn = INFINITY };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("scaled fallback accepted a nonfinite normalized contraction");
  }
  closure = sentinel;
  comoving = (ghl_m1_comoving){ .J = scaled_J,
                                .HD = { INFINITY, 0.0, 0.0 },
                                .HU = { 1.0, 0.0, 0.0 } };
  if(ghl_m1_closure_private_finalize_fallback(
           ghl_success, &comoving, &candidate, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("scaled fallback accepted a nonfinite normalized contraction norm");
  }
}

static void check_private_eulerian_pressure_construction(void) {
  const double gammaUU[3][3]
        = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
  const double direction[3] = { 1.0, 0.0, 0.0 };
  const double factors[] = { 1.0, 1.0 + 64.0 * DBL_EPSILON };
  for(size_t k = 0; k < sizeof(factors) / sizeof(factors[0]); ++k) {
    double P[3][3], chi;
    if(ghl_m1_closure_private_construct_eulerian_pressure(
             2.0, gammaUU, factors[k], direction, P, &chi)
             != ghl_success
       || chi != 1.0 || P[0][0] != 2.0 || P[1][1] != 0.0 || P[2][2] != 0.0) {
      fail("Eulerian pressure did not clip exact or tolerated streaming to chi=1");
      return;
    }
  }
  const double overflowing_metric[3][3]
        = { { DBL_MAX, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
  const double zero_direction[3] = { 0.0, 0.0, 0.0 };
  double P[3][3], chi;
  if(ghl_m1_closure_private_construct_eulerian_pressure(
           8.0, overflowing_metric, 0.0, zero_direction, P, &chi)
     != ghl_error_m1_invalid_state) {
    fail("Eulerian pressure constructor accepted an overflowing component");
  }
}

typedef struct {
  pthread_mutex_t mutex;
  pthread_cond_t condition;
  int ready;
  bool start;
  ghl_m1_parameters params;
  ghl_metric_quantities metric;
  ghl_primitive_quantities prims;
  ghl_m1_rad_state state;
} closure_contention_context;

typedef struct {
  closure_contention_context *context;
  ghl_error_codes_t result;
} closure_contention_worker;

static void *record_concurrent_workspace_failure(void *const opaque) {
  closure_contention_worker *const worker = opaque;
  closure_contention_context *const context = worker->context;
  pthread_mutex_lock(&context->mutex);
  ++context->ready;
  pthread_cond_broadcast(&context->condition);
  while(!context->start) {
    pthread_cond_wait(&context->condition, &context->mutex);
  }
  pthread_mutex_unlock(&context->mutex);

  ghl_m1_closure closure = { 0 };
  worker->result = ghl_m1_compute_closure_with_primitives(
        &context->params, &context->metric, &context->prims, &context->state, &closure);
  return NULL;
}

/* Two simultaneous workspace failures publish the same stage through the
 * packed snapshot update. This keeps the stage/reason pair valid while
 * exercising the compare-exchange retry when both writers read one snapshot. */
static void check_concurrent_closure_failure_snapshot(void) {
  closure_contention_context context = { 0 };
  setup_parameters(&context.params);
  m1_setup_flat_metric(&context.metric);
  /* Keep the metric and state valid while making the lowered closure tensor
   * unrepresentable. Both threads then fail at the same workspace publication
   * site after competing to update the packed diagnostic snapshot. */
  context.metric.gammaDD[0][0] = 1.0e308;
  context.metric.gammaUU[0][0] = 1.0e-308;
  context.metric.detgamma = 1.0e308;
  context.metric.sqrt_detgamma = 1.0e154;
  context.state.E = 1.0e154;
  context.state.F[0] = 0.99e308;
  if(!ghl_m1_metric_is_symmetric_spd(&context.metric)
     || ghl_m1_validate_realizability_state(
              &context.params, &context.metric, &context.state, 128.0, NULL)
              != ghl_success) {
    fail("closure contention fixture is not a valid metric and radiation state");
    return;
  }

  if(pthread_mutex_init(&context.mutex, NULL) != 0) {
    fail("could not initialize closure contention synchronization");
    return;
  }
  if(pthread_cond_init(&context.condition, NULL) != 0) {
    pthread_mutex_destroy(&context.mutex);
    fail("could not initialize closure contention synchronization");
    return;
  }
  pthread_t threads[2];
  closure_contention_worker workers[2]
        = { { .context = &context }, { .context = &context } };
  int created = 0;
  for(; created < 2; ++created) {
    if(pthread_create(
             &threads[created], NULL, record_concurrent_workspace_failure,
             &workers[created])
       != 0) {
      fail("could not start closure contention worker");
      break;
    }
  }
  if(created == 2) {
    pthread_mutex_lock(&context.mutex);
    while(context.ready != 2) {
      pthread_cond_wait(&context.condition, &context.mutex);
    }
    context.start = true;
    pthread_cond_broadcast(&context.condition);
    pthread_mutex_unlock(&context.mutex);
  }
  else {
    pthread_mutex_lock(&context.mutex);
    context.start = true;
    pthread_cond_broadcast(&context.condition);
    pthread_mutex_unlock(&context.mutex);
  }
  for(int i = 0; i < created; ++i) {
    pthread_join(threads[i], NULL);
  }
  if(created == 2) {
    if(workers[0].result != ghl_error_m1_invalid_state
       || workers[1].result != ghl_error_m1_invalid_state) {
      fail("concurrent workspace failures returned an unexpected error");
    }
    ghl_m1_closure_failure_stage_t stage;
    int reason;
    ghl_m1_get_last_closure_failure_stage(&stage);
    ghl_m1_get_last_closure_validation_reason(&reason);
    if(stage != ghl_m1_closure_failure_workspace || reason != 0) {
      fail("concurrent closure writers published a torn stage/reason pair");
    }
  }
  pthread_cond_destroy(&context.condition);
  pthread_mutex_destroy(&context.mutex);
}

static void check_both_closure_endpoints(void) {
  ghl_m1_parameters params;
  if(ghl_m1_initialize(DBL_MIN, DBL_MIN, 1.0e-8, 1.0e-6, 1.0e-12, 20, 1.0e-10, &params)
     != ghl_success) {
    fail("endpoint closure parameters were rejected");
    return;
  }
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  const ghl_primitive_quantities prims = { .u0 = 1.0 };
  const ghl_m1_rad_state endpoint_states[]
        = { { .E = 1.0, .F = { 0.0, 0.0, 0.0 } }, { .E = 1.0, .F = { 1.0, 0.0, 0.0 } } };
  for(size_t i = 0; i < sizeof(endpoint_states) / sizeof(endpoint_states[0]); ++i) {
    ghl_m1_closure closure = { 0 };
    if(ghl_m1_compute_closure_with_primitives(
             &params, &metric, &prims, &endpoint_states[i], &closure)
             != ghl_success
       || closure.solve_status != ghl_m1_closure_solve_converged
       || closure.root_iterations != 0) {
      fail("exact Minerbo endpoint did not converge without root iterations");
      return;
    }
    const double expected_xi = i == 0 ? 0.0 : 1.0;
    if(fabs(closure.xi - expected_xi) > 8.0 * DBL_EPSILON) {
      fail("exact Minerbo endpoint published the wrong flux factor");
      return;
    }
  }
}

/* Flat metric with zero shift and unit lapse, so prims->vU is the Eulerian
 * three-velocity and the reduced flux factor is the Euclidean |F|/E. */
static void setup_velocity(ghl_primitive_quantities *restrict prims, const double V[3]) {
  memset(prims, 0, sizeof(*prims));
  prims->vU[0] = V[0];
  prims->vU[1] = V[1];
  prims->vU[2] = V[2];
}

static int
published_tensor_is_exactly_symmetric(const ghl_m1_closure *restrict closure) {
  for(int i = 0; i < 3; ++i) {
    for(int j = i + 1; j < 3; ++j) {
      if(closure->P[i][j] != closure->P[j][i]) {
        return 0;
      }
    }
  }
  return 1;
}

/* Sweep transverse flux directions at a speed inside the PSD regime. Every
 * state is admissible, so every call must publish, and every published tensor
 * must be exactly symmetric. */
static void check_psd_fallback_publishes_symmetric_tensors(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);

  const double V[3] = { 0.8, 0.0, 0.0 };
  ghl_primitive_quantities prims;
  setup_velocity(&prims, V);

  ghl_m1_reset_closure_counters();

  int sampled = 0;
  int fallbacks = 0;
  for(int step = 0; step < 720; ++step) {
    const double phi = 2.0 * M_PI * (double)step / 720.0;
    for(int k = 1; k <= 9; ++k) {
      const double flux_factor = 0.1 * (double)k;
      ghl_m1_rad_state rad_state;
      rad_state.E = 1.0;
      /* Transverse to V: no x component. */
      rad_state.F[0] = 0.0;
      rad_state.F[1] = flux_factor * cos(phi);
      rad_state.F[2] = flux_factor * sin(phi);

      ghl_m1_closure closure;
      memset(&closure, 0, sizeof(closure));
      const ghl_error_codes_t error = ghl_m1_compute_closure_minerbo(
            &params, &metric, &prims, &rad_state, &closure);
      if(error != ghl_success) {
        fail("admissible transverse-flux state failed to publish a closure");
        return;
      }
      if(!published_tensor_is_exactly_symmetric(&closure)) {
        fail("published pressure tensor is not exactly symmetric");
        return;
      }
      sampled++;
      if(!closure.four_point_compatibility) {
        fallbacks++;
      }
    }
  }

  if(sampled == 0) {
    fail("transverse sweep produced no samples");
    return;
  }
  if(fallbacks == 0) {
    fail("transverse sweep never reached the PSD admissibility fallback; this "
         "test no longer covers the path it owns");
    return;
  }

  ghl_m1_closure_counters counters;
  ghl_m1_get_closure_counters(&counters);
  if(counters.admissibility_fallback_psd != (unsigned long long)fallbacks) {
    fail("admissibility_fallback_psd does not match the observed substitutions");
  }
  if(counters.admissibility_fallback_zero_flux != 0ULL) {
    fail("finite-flux sweep incremented the zero-flux fallback counter");
  }
}

/* An exact-zero Eulerian flux with a moving fluid takes the same fallback for a
 * different, documented reason, and must be accounted separately. */
static void check_zero_flux_fallback_is_counted_separately(const double energy) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  params.E_floor = fmin(params.E_floor, energy);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);

  const double V[3] = { 0.3, 0.1, -0.2 };
  ghl_primitive_quantities prims;
  setup_velocity(&prims, V);

  ghl_m1_rad_state rad_state;
  rad_state.E = energy;
  rad_state.F[0] = 0.0;
  rad_state.F[1] = 0.0;
  rad_state.F[2] = 0.0;

  ghl_m1_reset_closure_counters();

  ghl_m1_closure closure;
  memset(&closure, 0, sizeof(closure));
  if(ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &rad_state, &closure)
     != ghl_success) {
    fail("exact-zero-flux state failed to publish a closure");
    return;
  }
  if(!published_tensor_is_exactly_symmetric(&closure)) {
    fail("zero-flux published pressure tensor is not exactly symmetric");
  }
  if(closure.four_point_compatibility) {
    fail("exact-zero-flux state with a moving fluid did not take the fallback; "
         "this test no longer covers the path it owns");
    return;
  }

  ghl_m1_closure_counters counters;
  ghl_m1_get_closure_counters(&counters);
  if(counters.admissibility_fallback_zero_flux != 1ULL) {
    fail("zero-flux fallback was not counted");
  }
  if(counters.admissibility_fallback_psd != 0ULL) {
    fail("zero-flux fallback incremented the PSD fallback counter");
  }
}

/* Pin the documented regime boundary in both directions. */
static void check_psd_regime_boundary(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);

  const double flux_factors[] = { 0.05, 0.2, 0.4, 0.6, 0.8, 0.95 };
  const int flux_count = (int)(sizeof(flux_factors) / sizeof(flux_factors[0]));

  for(int speed_step = 1; speed_step <= 9; ++speed_step) {
    const double speed = 0.1 * (double)speed_step;
    const double V[3] = { speed, 0.0, 0.0 };
    ghl_primitive_quantities prims;
    setup_velocity(&prims, V);

    for(int k = 0; k < flux_count; ++k) {
      const double flux_factor = flux_factors[k];

      /* Flux aligned with and opposed to the flow must stay on the primary
       * construction at every speed tested. */
      for(int sign = -1; sign <= 1; sign += 2) {
        ghl_m1_rad_state rad_state;
        rad_state.E = 1.0;
        rad_state.F[0] = (double)sign * flux_factor;
        rad_state.F[1] = 0.0;
        rad_state.F[2] = 0.0;
        ghl_m1_closure closure;
        memset(&closure, 0, sizeof(closure));
        if(ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &rad_state, &closure)
           != ghl_success) {
          fail("aligned-flux state failed to publish a closure");
          return;
        }
        if(!published_tensor_is_exactly_symmetric(&closure)) {
          fail("aligned-flux published pressure tensor is not exactly symmetric");
          return;
        }
        if(!closure.four_point_compatibility) {
          fail("flux aligned with the fluid velocity reached the PSD fallback; "
               "the documented regime has widened");
          return;
        }
      }

      /* Transverse flux: below the documented onset the primary construction
       * must still be published. */
      if(speed <= 0.4) {
        ghl_m1_rad_state rad_state;
        rad_state.E = 1.0;
        rad_state.F[0] = 0.0;
        rad_state.F[1] = flux_factor;
        rad_state.F[2] = 0.0;
        ghl_m1_closure closure;
        memset(&closure, 0, sizeof(closure));
        if(ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &rad_state, &closure)
           != ghl_success) {
          fail("transverse-flux state below the onset failed to publish a closure");
          return;
        }
        if(!closure.four_point_compatibility) {
          fail("transverse flux below half light speed reached the PSD fallback; "
               "the documented regime has widened");
          return;
        }
      }
    }
  }
}

/* Crossing the E^2 overflow boundary must not change the dimensionless
 * closure. A nonidentity metric distinguishes F_i V^i from F^i V^i. */
static void check_large_energy_scaling(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  metric.gammaDD[0][0] = 4.0;
  metric.gammaUU[0][0] = 0.25;
  metric.detgamma = 4.0;
  metric.sqrt_detgamma = 2.0;
  const double V[3] = { 0.1, 0.0, 0.0 };
  ghl_primitive_quantities prims;
  setup_velocity(&prims, V);

  /* 2^520 has a finite pressure tensor but its square exceeds DBL_MAX. */
  const double energies[] = { 1.0, 0x1p520 };
  ghl_m1_closure closures[2];
  for(int k = 0; k < 2; ++k) {
    const ghl_m1_rad_state state
          = { .E = energies[k], .F = { 0.8 * energies[k], 0.0, 0.0 } };
    if(ghl_m1_compute_closure_with_primitives(
             &params, &metric, &prims, &state, &closures[k])
       != ghl_success) {
      fail("energy-scaled closure failed to publish");
      return;
    }
  }
  const double tolerance = params.closure_root_residual_tolerance;
  if(fabs(closures[0].xi - closures[1].xi) > tolerance
     || fabs(closures[0].chi - closures[1].chi) > tolerance) {
    fail("dimensionless closure changed across the E^2 overflow boundary");
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      if(fabs(closures[0].P[i][j] - closures[1].P[i][j] / energies[1]) > tolerance) {
        fail("normalized pressure changed across the E^2 overflow boundary");
      }
    }
  }
}

/* An accepted interval tolerance below the double spacing near the root must
 * still converge rather than exhaust the iteration budget. */
static void check_sub_ulp_root_tolerance_converges(void) {
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  const double V[3] = { 0.45, 0.05, 0.0 };
  ghl_primitive_quantities prims;
  setup_velocity(&prims, V);
  const ghl_m1_rad_state state = { .E = 1.0, .F = { 0.3, 0.4, 0.0 } };

  ghl_m1_parameters params;
  setup_parameters(&params);
  if(ghl_m1_set_closure_solver_controls(1.0e-12, 100, &params) != ghl_success) {
    fail("closure solver rejected the reference tolerance");
    return;
  }
  ghl_m1_closure reference;
  if(ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &state, &reference)
           != ghl_success
     || reference.solve_status != ghl_m1_closure_solve_converged) {
    fail("reference closure root did not converge");
    return;
  }

  const double tolerances[] = { 1.0e-16, 1.0e-17 };
  for(int k = 0; k < 2; ++k) {
    setup_parameters(&params);
    if(ghl_m1_set_closure_solver_controls(tolerances[k], 100, &params) != ghl_success) {
      fail("closure solver rejected a sub-ulp tolerance");
      return;
    }
    ghl_m1_closure closure;
    if(ghl_m1_compute_closure_minerbo(&params, &metric, &prims, &state, &closure)
       != ghl_success) {
      fail("sub-ulp tolerance closure failed to publish");
      return;
    }
    if(closure.solve_status != ghl_m1_closure_solve_converged
       || closure.root_iterations >= 100) {
      fail("sub-ulp tolerance closure root exhausted its iterations");
    }
    if(fabs(closure.xi - reference.xi) > 1.0e-12) {
      fail("sub-ulp tolerance closure root moved");
    }
  }
}

/* Coordinate permutations of the same sheared metric must reject an
 * inaccurate finite-flux trace on every axis, without publishing a candidate. */
static void check_trace_failure_on_each_axis(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  const double shear = 1.0 - 0x1p-10;
  const double determinant = 1.0 - shear * shear;
  const ghl_primitive_quantities prims = { .u0 = 1.0 };
  for(int axis = 0; axis < 3; ++axis) {
    const int other = (axis + 1) % 3;
    ghl_metric_quantities metric;
    m1_setup_flat_metric(&metric);
    metric.gammaDD[axis][other] = metric.gammaDD[other][axis] = shear;
    metric.gammaUU[axis][axis] = metric.gammaUU[other][other] = 1.0 / determinant;
    metric.gammaUU[axis][other] = metric.gammaUU[other][axis] = -shear / determinant;
    metric.detgamma = determinant;
    metric.sqrt_detgamma = sqrt(determinant);
    ghl_m1_rad_state state = { .E = 1.0 };
    state.F[axis] = 0.3 * sqrt(determinant);
    const ghl_m1_closure sentinel = { .xi = -1.0 };
    ghl_m1_closure closure = sentinel;
    const ghl_error_codes_t error = ghl_m1_compute_closure_with_primitives(
          &params, &metric, &prims, &state, &closure);
    ghl_m1_closure_failure_stage_t stage;
    int reason;
    ghl_m1_get_last_closure_failure_stage(&stage);
    ghl_m1_get_last_closure_validation_reason(&reason);
    if(error != ghl_error_m1_invalid_state
       || stage != ghl_m1_closure_failure_tensor_validation
       || reason != ghl_m1_closure_validation_trace
       || !m1_closure_identical(&closure, &sentinel)) {
      fail("axis-permuted trace failure changed its diagnostic or published output");
    }
  }
}

static void check_large_energy_failure_paths(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  const struct {
    double energy, metric_scale, speed, flux_factor;
    ghl_error_codes_t error;
    ghl_m1_closure_failure_stage_t stage;
  } cases[] = {
    /* The thin dyad is zero; only the covariant thick pressure overflows. */
    { 0x1p520, 1.0e200, 0.0, 0.0, ghl_error_m1_invalid_state,
      ghl_m1_closure_failure_workspace },
    /* Homogeneous normalization must retain the endpoint residual rejection
     * already required for a unit-energy, nearly luminal flow. */
    { 0x1p520, 1.0, 0x1.fffffffffffffp-1, 0.3, ghl_error_m1_closure_residual_too_large,
      ghl_m1_closure_failure_residual_gate },
    /* The zero-flux fallback pressure is finite, but its comoving energy
     * exceeds DBL_MAX. It must fail atomically instead of publishing. */
    { DBL_MAX, 1.0, 0.3, 0.0, ghl_error_m1_invalid_state,
      ghl_m1_closure_failure_tensor_validation }
  };
  for(size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
    ghl_metric_quantities metric;
    m1_setup_flat_metric(&metric);
    metric.gammaDD[0][0] = cases[i].metric_scale;
    metric.gammaUU[0][0] = 1.0 / cases[i].metric_scale;
    metric.detgamma = cases[i].metric_scale;
    metric.sqrt_detgamma = sqrt(cases[i].metric_scale);
    const ghl_primitive_quantities prims = { .vU = { cases[i].speed, 0.0, 0.0 } };
    const ghl_m1_rad_state state
          = { .E = cases[i].energy,
              .F = { cases[i].energy * cases[i].flux_factor, 0.0, 0.0 } };
    const ghl_m1_closure sentinel = { .xi = -1.0 };
    ghl_m1_closure closure = sentinel;
    const ghl_error_codes_t error = ghl_m1_compute_closure_with_primitives(
          &params, &metric, &prims, &state, &closure);
    ghl_m1_closure_failure_stage_t stage;
    ghl_m1_get_last_closure_failure_stage(&stage);
    if(error != cases[i].error || stage != cases[i].stage
       || !m1_closure_identical(&closure, &sentinel)) {
      fail("large-energy closure failure changed its stage or published output");
    }
  }
}

static void check_workspace_failure_paths(void) {
  ghl_m1_parameters params;
  setup_parameters(&params);
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  const ghl_m1_rad_state zero_flux = { .E = 1.0 };
  const ghl_primitive_quantities nonfinite_velocity = { .vU = { NAN, 0.0, 0.0 } };
  const ghl_m1_closure sentinel = { .xi = -1.0 };
  ghl_m1_closure closure = sentinel;
  if(ghl_m1_compute_closure_with_primitives(
           &params, &metric, &nonfinite_velocity, &zero_flux, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("workspace propagated a nonfinite primitive velocity incorrectly");
  }

  const double large_shift = DBL_MAX * 0.5;
  metric.betaU[0] = large_shift;
  const ghl_primitive_quantities cancelling_velocity
        = { .vU = { -large_shift, 0.0, 0.0 } };
  const ghl_m1_rad_state overflowing_time_flux = { .E = 8.0, .F = { 4.0, 0.0, 0.0 } };
  closure = sentinel;
  if(ghl_m1_compute_closure_with_primitives(
           &params, &metric, &cancelling_velocity, &overflowing_time_flux, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("workspace accepted an overflowing time component of the flux covector");
  }
  ghl_m1_closure_failure_stage_t stage;
  ghl_m1_get_last_closure_failure_stage(&stage);
  if(stage != ghl_m1_closure_failure_workspace) {
    fail("overflowing time flux recorded the wrong failure stage");
  }

  metric.betaU[0] = 1.0e200;
  const ghl_primitive_quantities cancelling_large_shift
        = { .vU = { -1.0e200, 0.0, 0.0 } };
  closure = sentinel;
  if(ghl_m1_compute_closure_with_primitives(
           &params, &metric, &cancelling_large_shift, &zero_flux, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("workspace accepted an overflowing ordinary-range thick pressure");
  }
  ghl_m1_get_last_closure_failure_stage(&stage);
  if(stage != ghl_m1_closure_failure_workspace) {
    fail("ordinary thick-pressure overflow recorded the wrong failure stage");
  }

  const ghl_m1_rad_state large_energy_zero_flux = { .E = 0x1p520 };
  closure = sentinel;
  if(ghl_m1_compute_closure_with_primitives(
           &params, &metric, &cancelling_large_shift, &large_energy_zero_flux, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("workspace accepted an overflowing scaled thick pressure");
  }
  ghl_m1_get_last_closure_failure_stage(&stage);
  if(stage != ghl_m1_closure_failure_workspace) {
    fail("scaled thick-pressure overflow recorded the wrong failure stage");
  }

  m1_setup_flat_metric(&metric);
  metric.gammaDD[0][0] = 1.0e308;
  metric.gammaUU[0][0] = 1.0e-308;
  metric.detgamma = 1.0e308;
  metric.sqrt_detgamma = sqrt(metric.detgamma);
  const ghl_m1_rad_state ordinary_thin_overflow
        = { .E = 2.0, .F = { 1.0e154, 0.0, 0.0 } };
  const ghl_primitive_quantities rest = { .vU = { 0.0, 0.0, 0.0 } };
  closure = sentinel;
  if(ghl_m1_compute_closure_with_primitives(
           &params, &metric, &rest, &ordinary_thin_overflow, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("workspace accepted an overflowing ordinary-range thin pressure");
  }
  ghl_m1_get_last_closure_failure_stage(&stage);
  if(stage != ghl_m1_closure_failure_workspace) {
    fail("ordinary thin-pressure overflow recorded the wrong failure stage");
  }

  const ghl_m1_rad_state scaled_thin_overflow
        = { .E = 0x1p520, .F = { DBL_MAX, 0.0, 0.0 } };
  closure = sentinel;
  if(ghl_m1_compute_closure_with_primitives(
           &params, &metric, &rest, &scaled_thin_overflow, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("workspace accepted an overflowing scaled thin pressure");
  }
  ghl_m1_get_last_closure_failure_stage(&stage);
  if(stage != ghl_m1_closure_failure_workspace) {
    fail("scaled thin-pressure overflow recorded the wrong failure stage");
  }

  m1_setup_flat_metric(&metric);
  metric.lapse = 2.0;
  const ghl_m1_rad_state large_energy = { .E = DBL_MAX * 0.5 };
  closure = sentinel;
  if(ghl_m1_compute_closure_with_primitives(
           &params, &metric, &rest, &large_energy, &closure)
           != ghl_error_m1_invalid_state
     || !m1_closure_identical(&closure, &sentinel)) {
    fail("closure evaluation accepted an overflowing stress-energy trace term");
  }
  ghl_m1_get_last_closure_failure_stage(&stage);
  if(stage != ghl_m1_closure_failure_workspace) {
    fail("overflowing stress-energy trace recorded the wrong failure stage");
  }
}

int main(void) {
  check_forced_closure_snapshot_retry();
  check_private_comoving_invariant();
  check_private_fallback_invariants();
  check_private_eulerian_pressure_construction();
  check_private_root_driver();
  check_trace_failure_on_each_axis();
  check_large_energy_failure_paths();
  check_workspace_failure_paths();
  check_concurrent_closure_failure_snapshot();
  check_both_closure_endpoints();
  check_large_energy_scaling();
  check_sub_ulp_root_tolerance_converges();
  check_psd_fallback_publishes_symmetric_tensors();
  /* The fallback's invariant H^2/J^2 must survive both square-range limits. */
  check_zero_flux_fallback_is_counted_separately(1.0);
  check_zero_flux_fallback_is_counted_separately(0x1p520);
  check_zero_flux_fallback_is_counted_separately(0x1p-520);
  check_psd_regime_boundary();

  if(failures != 0) {
    fprintf(stderr, "unit_test_m1_closure_fallback: %d check(s) failed\n", failures);
    return 1;
  }

  ghl_info(
        "unit_test_m1_closure_fallback: admissibility fallback symmetry, "
        "accounting, and regime checks passed\n");
  return 0;
}
