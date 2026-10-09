#ifndef UNIT_TESTS_M1_THCM1_RUSANOV_FIXTURE_H_
#define UNIT_TESTS_M1_THCM1_RUSANOV_FIXTURE_H_

/*
 * Test-local evaluators for the retained paired Rusanov observations.
 *
 * The fixture inputs are the complete operands consumed by the corresponding
 * public GRHayL operation.  The neutrino fixture stores caller-prepared metric,
 * state, closure, number-current, velocity, and speed operands.  Its closure
 * and physical-flux preparation are therefore part of the serialized common
 * input, not an independent THC physical-flux validation.  The generic fixture
 * stores the public E/F state and physical-flux operands directly.
 */

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../m1_test_utils.h"
#include "m1_thcm1_transport_fixture.h"

#define M1_THCM1_RUSANOV_POLICY "strict_relative_2e-12_propagated_response_v1"
#define M1_THCM1_RUSANOV_CURRENT_POLICY \
  "strict_relative_2e-12_propagated_response_current_v2"
#define M1_THCM1_RUSANOV_CURRENT_OPERATION     "neutrino_rusanov_flux_current_v2"
#define M1_THCM1_RUSANOV_RECORD_COUNT          1024
#define M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT  54
#define M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT 5
#define M1_THCM1_RUSANOV_GENERIC_INPUT_COUNT   17
#define M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT  4

static inline int m1_thcm1_rusanov_pair_id_matches_operation(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict operation) {
  if(record == NULL || operation == NULL
     || !m1_thcm1_fixture_pair_metadata_valid(record)) {
    return 0;
  }
  char expected[M1_THCM1_FIXTURE_TEXT_MAX] = { 0 };
  const int written
        = snprintf(expected, sizeof(expected), "%s:%s", operation, record->case_id);
  return written >= 0 && (size_t)written < sizeof(expected)
         && strcmp(record->pair_id, expected) == 0;
}

static inline int m1_thcm1_rusanov_compare_baseline_response(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict policy,
      const char *restrict operation,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size) {
  const int supported_policy
        = policy != NULL
          && (strcmp(policy, M1_THCM1_RUSANOV_POLICY) == 0
              || strcmp(policy, M1_THCM1_RUSANOV_CURRENT_POLICY) == 0);
  if(record == NULL || !supported_policy
     || !m1_thcm1_rusanov_pair_id_matches_operation(record, operation)) {
    m1_thcm1_fixture_set_error(
          error, error_size, "invalid Rusanov baseline-response metadata");
    return 0;
  }
  return m1_thcm1_fixture_compare_baseline_response(
        record, policy, computed_normalization, computed_baseline, report, error,
        error_size);
}

/* Strict paired replay rule for the retained Rusanov corpora.  Both endpoint
 * roles and the paired response are compared against the retained values
 * with the declared strict_relative policy, reusing the shared comparator so
 * the endpoint scales propagate identically to the transport and
 * stress-energy owners.  The operation predicate mirrors the baseline
 * envelope entry. This propagated response rule does not independently test
 * derivative accuracy. */
static inline int m1_thcm1_rusanov_compare_paired_strict_relative(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict policy,
      const char *restrict operation,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      const double *restrict computed_perturbed,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size) {
  const int supported_policy
        = policy != NULL
          && (strcmp(policy, M1_THCM1_RUSANOV_POLICY) == 0
              || strcmp(policy, M1_THCM1_RUSANOV_CURRENT_POLICY) == 0);
  if(record == NULL || !supported_policy
     || !m1_thcm1_rusanov_pair_id_matches_operation(record, operation)
     || computed_normalization == NULL || computed_baseline == NULL
     || computed_perturbed == NULL) {
    m1_thcm1_fixture_set_error(
          error, error_size, "invalid Rusanov strict paired comparison argument");
    return 0;
  }
  return m1_thcm1_fixture_compare_paired_strict_relative(
        record, computed_normalization, computed_baseline, computed_perturbed, report,
        error, error_size);
}

static inline void m1_thcm1_rusanov_metric(
      const double input[M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT],
      ghl_metric_quantities *restrict metric) {
  *metric = (ghl_metric_quantities){ 0 };
  metric->lapse = input[1];
  metric->lapseinv = 1.0 / metric->lapse;
  metric->lapseinv2 = metric->lapseinv * metric->lapseinv;
  for(int i = 0; i < 3; ++i) {
    metric->betaU[i] = input[2 + i];
    metric->gammaDD[i][i] = input[5 + i];
    metric->gammaUU[i][i] = input[8 + i];
  }
  metric->detgamma = input[11];
  metric->sqrt_detgamma = input[12];
}

static inline void m1_thcm1_rusanov_state(
      const double *restrict values,
      ghl_m1_neutrino_state *restrict state) {
  state->N = values[0];
  state->E = values[1];
  for(int i = 0; i < 3; ++i) {
    state->F[i] = values[2 + i];
  }
}

static inline void m1_thcm1_rusanov_closure(
      const double *restrict values,
      ghl_m1_closure *restrict closure) {
  *closure = (ghl_m1_closure){ 0 };
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      closure->P[i][j] = values[3 * i + j];
    }
  }
}

int m1_thcm1_rusanov_check_neutrino_fixture(
      const char *restrict fixture_dir,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      char *restrict error,
      const size_t error_size);

/* The current shard deliberately has a separate operation identity.  It
 * retains the same 1024 face pairs but requires every number-current operand
 * to be nonzero and proves that radiation_N perturbations reach the number
 * output.  The public evaluator is still the only test-side computation. */
int m1_thcm1_rusanov_check_current_fixture(
      const char *restrict fixture_dir,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      char *restrict error,
      const size_t error_size);

int m1_thcm1_rusanov_check_generic_fixture(
      const char *restrict fixture_dir,
      char *restrict error,
      const size_t error_size);

/* Negative assertions use a stack-owned record copy: a generic perturbed
 * reference defect, reversed generic reference endpoints, or a corrupted
 * perturbed number-current operand. No payload bytes are rewritten. */
#define M1_THCM1_RUSANOV_REGRESSION_PERTURBED_OUTPUT 1
#define M1_THCM1_RUSANOV_REGRESSION_SWAPPED_ENDPOINTS 2
#define M1_THCM1_RUSANOV_REGRESSION_CORRUPTED_CURRENT 3

int m1_thcm1_rusanov_check_fixture_regression(
      const char *restrict fixture_dir,
      const int mutation_mode,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      char *restrict error,
      const size_t error_size);

#endif /* UNIT_TESTS_M1_THCM1_RUSANOV_FIXTURE_H_ */
