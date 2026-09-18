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

#include "m1_test_utils.h"
#include "m1_thcm1_transport_fixture.h"

#define M1_THCM1_RUSANOV_POLICY \
  "strict_relative_2e-12_propagated_response_v1"
#define M1_THCM1_RUSANOV_CURRENT_POLICY \
  "strict_relative_2e-12_propagated_response_current_v2"
#define M1_THCM1_RUSANOV_CURRENT_OPERATION "neutrino_rusanov_flux_current_v2"
#define M1_THCM1_RUSANOV_RECORD_COUNT 1024
#define M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT 54
#define M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT 5
#define M1_THCM1_RUSANOV_GENERIC_INPUT_COUNT 17
#define M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT 4

static inline int m1_thcm1_rusanov_pair_id_matches_operation(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict operation) {
  if(record == NULL || operation == NULL ||
     !m1_thcm1_fixture_pair_metadata_valid(record)) return 0;
  char expected[M1_THCM1_FIXTURE_TEXT_MAX] = {0};
  const int written = snprintf(expected, sizeof(expected), "%s:%s",
                               operation, record->case_id);
  return written >= 0 && (size_t)written < sizeof(expected) &&
         strcmp(record->pair_id, expected) == 0;
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
  const int supported_policy = policy != NULL &&
      (strcmp(policy, M1_THCM1_RUSANOV_POLICY) == 0 ||
       strcmp(policy, M1_THCM1_RUSANOV_CURRENT_POLICY) == 0);
  if(record == NULL || !supported_policy ||
     !m1_thcm1_rusanov_pair_id_matches_operation(record, operation)) {
    m1_thcm1_fixture_set_error(
        error, error_size, "invalid Rusanov baseline-response metadata");
    return 0;
  }
  return m1_thcm1_fixture_compare_baseline_response(
      record, policy, computed_normalization, computed_baseline,
      report, error, error_size);
}

static inline void m1_thcm1_rusanov_metric(
      const double input[M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT],
      ghl_metric_quantities *restrict metric) {
  *metric = (ghl_metric_quantities){0};
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
  for(int i = 0; i < 3; ++i) state->F[i] = values[2 + i];
}

static inline void m1_thcm1_rusanov_closure(
      const double *restrict values,
      ghl_m1_closure *restrict closure) {
  *closure = (ghl_m1_closure){0};
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      closure->P[i][j] = values[3 * i + j];
}

static inline int m1_thcm1_rusanov_eval_neutrino(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const m1_thcm1_fixture_record *restrict record,
      const int role,
      double output[M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT]) {
  if(m1_params == NULL || nu_params == NULL || record == NULL ||
     record->input_count != M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT ||
     record->output_count != M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT ||
     output == NULL)
    return 0;
  if(role < 0 || role > 1) return 0;
  const double *input = role == 0 ? record->baseline_input :
                                     record->perturbed_input;
  if(input == NULL || !isfinite(input[0]) || input[0] < 0.0 ||
     input[0] > 2.0 || floor(input[0]) != input[0])
    return 0;
  ghl_metric_quantities metric;
  m1_thcm1_rusanov_metric(input, &metric);
  ghl_m1_neutrino_state state_L, state_R;
  m1_thcm1_rusanov_state(&input[13], &state_L);
  m1_thcm1_rusanov_state(&input[18], &state_R);
  ghl_m1_closure closure_L, closure_R;
  m1_thcm1_rusanov_closure(&input[23], &closure_L);
  m1_thcm1_rusanov_closure(&input[32], &closure_R);
  double flux_F[3] = {0.0, 0.0, 0.0};
  const ghl_error_codes_t status = ghl_m1_compute_neutrino_rusanov_flux(
             m1_params, nu_params, &metric,
             (ghl_m1_direction_t)(int)input[0], &state_L, &state_R,
             &closure_L, &closure_R, &input[41], &input[44], &input[47],
             &input[50], input[53], &output[0], &output[1], flux_F);
  if(status != ghl_success) return 0;
  output[2] = flux_F[0];
  output[3] = flux_F[1];
  output[4] = flux_F[2];
  return 1;
}

static inline void m1_thcm1_rusanov_neutrino_normalization(
      const m1_thcm1_fixture_record *restrict record,
      double normalization[M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT]) {
  for(int component = 0; component < 5; ++component)
    normalization[component] = 0.0;
  for(int role = 0; role < 2; ++role) {
    const double *input = role == 0 ? record->baseline_input :
                                       record->perturbed_input;
    for(int component = 0; component < 5; ++component) {
      const double state_L = input[13 + component];
      const double state_R = input[18 + component];
      const double scale = fmax(
          fmax(fabs(state_L), fabs(state_R)),
          fmax(input[53] * fabs(state_L), input[53] * fabs(state_R)));
      normalization[component] = fmax(
          normalization[component],
          fmax(1.0e-300, input[12] * scale));
    }
  }
}

static inline int m1_thcm1_rusanov_eval_generic(
      const m1_thcm1_fixture_record *restrict record,
      const int role,
      double output[M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT]) {
  if(record == NULL ||
     record->input_count != M1_THCM1_RUSANOV_GENERIC_INPUT_COUNT ||
     record->output_count != M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT ||
     output == NULL)
    return 0;
  if(role < 0 || role > 1) return 0;
  const double *input = role == 0 ? record->baseline_input :
                                     record->perturbed_input;
  if(input == NULL) return 0;
  const ghl_m1_rad_state state_L = {
      .E = input[0], .F = {input[1], input[2], input[3]}};
  const ghl_m1_rad_state state_R = {
      .E = input[4], .F = {input[5], input[6], input[7]}};
  return ghl_m1_compute_rusanov_flux(
             &state_L, &state_R, input[8], &input[9], input[12],
             &input[13], input[16], &output[0], &output[1]) == ghl_success;
}

static inline void m1_thcm1_rusanov_generic_normalization(
      const m1_thcm1_fixture_record *restrict record,
      double normalization[M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT]) {
  for(int component = 0; component < 4; ++component)
    normalization[component] = 0.0;
  for(int role = 0; role < 2; ++role) {
    const double *input = role == 0 ? record->baseline_input :
                                       record->perturbed_input;
    for(int component = 0; component < 4; ++component) {
      const double state_L = input[component];
      const double state_R = input[4 + component];
      const double scale = fmax(
          fmax(fabs(state_L), fabs(state_R)),
          fmax(input[16] * fabs(state_L), input[16] * fabs(state_R)));
      normalization[component] = fmax(normalization[component],
                                      fmax(1.0e-300, scale));
    }
  }
}

static inline int m1_thcm1_rusanov_join_path(
      const char *restrict directory,
      const char *restrict name,
      char **restrict path) {
  if(directory == NULL || name == NULL || path == NULL) return 0;
  const size_t directory_length = strlen(directory);
  const size_t name_length = strlen(name);
  *path = (char *)malloc(directory_length + name_length + 2);
  if(*path == NULL) return 0;
  memcpy(*path, directory, directory_length);
  (*path)[directory_length] = '/';
  memcpy(*path + directory_length + 1, name, name_length + 1);
  return 1;
}

static inline int m1_thcm1_rusanov_check_neutrino_fixture(
      const char *restrict fixture_dir,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      char *restrict error,
      const size_t error_size) {
  char *path = NULL;
  if(!m1_thcm1_rusanov_join_path(fixture_dir, "rusanov_neutrino.dat", &path)) {
    m1_thcm1_fixture_set_error(error, error_size,
                               "neutrino Rusanov fixture path allocation failed");
    return 0;
  }
  m1_thcm1_fixture_collection collection = {0};
  const int loaded = m1_thcm1_fixture_load(
      path, "neutrino_rusanov_flux",
      M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT,
      M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT, &collection, error, error_size);
  free(path);
  if(!loaded) return 0;
  if(strcmp(collection.policy, M1_THCM1_RUSANOV_POLICY) != 0 ||
     collection.record_count != M1_THCM1_RUSANOV_RECORD_COUNT) {
    m1_thcm1_fixture_set_error(error, error_size,
                               "neutrino Rusanov fixture policy or count mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  size_t failures = 0;
  size_t direction_counts[3] = {0, 0, 0};
  for(size_t index = 0; index < collection.record_count; ++index) {
    const m1_thcm1_fixture_record *record = &collection.records[index];
    char expected_case_id[32];
    snprintf(expected_case_id, sizeof(expected_case_id), "face%04zu", index);
    if(strcmp(record->case_id, expected_case_id) != 0 ||
       !m1_thcm1_rusanov_pair_id_matches_operation(
           record, "neutrino_rusanov_flux") ||
       !isfinite(record->baseline_input[0]) ||
       floor(record->baseline_input[0]) != record->baseline_input[0] ||
       record->baseline_input[0] < 0.0 || record->baseline_input[0] > 2.0) {
      m1_thcm1_fixture_set_error(error, error_size,
                                 "neutrino Rusanov fixture face ordering mismatch");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    ++direction_counts[(int)record->baseline_input[0]];
    double normalization[5] = {0.0};
    double output[5] = {0.0};
    m1_thcm1_rusanov_neutrino_normalization(record, normalization);
    const int status_ok = m1_thcm1_rusanov_eval_neutrino(
        m1_params, nu_params, record, 0, output);
    m1_thcm1_fixture_comparison_report report = {0};
    char compare_error[256] = {0};
    if(!status_ok || !m1_thcm1_rusanov_compare_baseline_response(
           record, collection.policy, "neutrino_rusanov_flux", normalization,
           output, &report, compare_error, sizeof(compare_error))) {
      ++failures;
      if(failures == 1)
        m1_thcm1_fixture_set_error(
            error, error_size,
            compare_error[0] != '\0' ? compare_error :
                "neutrino Rusanov fixture evaluation failed");
    }
  }
  if(direction_counts[0] != 342 || direction_counts[1] != 341 ||
     direction_counts[2] != 341) {
    m1_thcm1_fixture_set_error(error, error_size,
                               "neutrino Rusanov fixture direction coverage mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  m1_thcm1_fixture_free(&collection);
  if(failures != 0) return 0;
  return 1;
}

/* The current shard deliberately has a separate operation identity.  It
 * retains the same 1024 face pairs but requires every number-current operand
 * to be nonzero and proves that radiation_N perturbations reach the number
 * output.  The public evaluator is still the only test-side computation. */
static inline int m1_thcm1_rusanov_check_current_fixture(
      const char *restrict fixture_dir,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      char *restrict error,
      const size_t error_size) {
  char *path = NULL;
  if(!m1_thcm1_rusanov_join_path(
         fixture_dir, "rusanov_neutrino_current.dat", &path)) {
    m1_thcm1_fixture_set_error(error, error_size,
                               "current Rusanov fixture path allocation failed");
    return 0;
  }
  m1_thcm1_fixture_collection collection = {0};
  const int loaded = m1_thcm1_fixture_load(
      path, M1_THCM1_RUSANOV_CURRENT_OPERATION,
      M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT,
      M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT, &collection, error, error_size);
  free(path);
  if(!loaded) return 0;
  if(strcmp(collection.policy, M1_THCM1_RUSANOV_CURRENT_POLICY) != 0 ||
     collection.record_count != M1_THCM1_RUSANOV_RECORD_COUNT) {
    m1_thcm1_fixture_set_error(error, error_size,
                               "current Rusanov fixture policy or count mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  size_t failures = 0;
  size_t direction_counts[3] = {0, 0, 0};
  size_t number_response_cases = 0;
  for(size_t index = 0; index < collection.record_count; ++index) {
    const m1_thcm1_fixture_record *record = &collection.records[index];
    char expected_case_id[32];
    snprintf(expected_case_id, sizeof(expected_case_id), "face%04zu", index);
    if(strcmp(record->case_id, expected_case_id) != 0 ||
       !m1_thcm1_rusanov_pair_id_matches_operation(
           record, M1_THCM1_RUSANOV_CURRENT_OPERATION) ||
       !isfinite(record->baseline_input[0]) ||
       floor(record->baseline_input[0]) != record->baseline_input[0] ||
       record->baseline_input[0] < 0.0 || record->baseline_input[0] > 2.0) {
      m1_thcm1_fixture_set_error(error, error_size,
                                 "current Rusanov fixture face ordering mismatch");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    ++direction_counts[(int)record->baseline_input[0]];
    for(int role = 0; role < 2; ++role) {
      const double *input = role == 0 ? record->baseline_input :
                                        record->perturbed_input;
      for(int component = 0; component < 12; ++component) {
        if(!isfinite(input[41 + component]) || input[41 + component] == 0.0)
          ++failures;
      }
      for(int component = 0; component < 3; ++component) {
        if(input[41 + component] != input[13] * input[47 + component] ||
           input[44 + component] != input[18] * input[50 + component])
          ++failures;
      }
    }
    if(strcmp(record->family, "radiation_N") == 0 &&
       record->baseline_output[0] != record->perturbed_output[0])
      ++number_response_cases;
    double normalization[5] = {0.0};
    double output[5] = {0.0};
    m1_thcm1_rusanov_neutrino_normalization(record, normalization);
    const int status_ok = m1_thcm1_rusanov_eval_neutrino(
        m1_params, nu_params, record, 0, output);
    m1_thcm1_fixture_comparison_report report = {0};
    char compare_error[256] = {0};
    if(!status_ok || !m1_thcm1_rusanov_compare_baseline_response(
           record, collection.policy, M1_THCM1_RUSANOV_CURRENT_OPERATION,
           normalization, output, &report, compare_error, sizeof(compare_error))) {
      ++failures;
      if(failures == 1)
        m1_thcm1_fixture_set_error(
            error, error_size,
            compare_error[0] != '\0' ? compare_error :
                "current Rusanov fixture evaluation failed");
    }
  }
  if(direction_counts[0] != 342 || direction_counts[1] != 341 ||
     direction_counts[2] != 341) {
    m1_thcm1_fixture_set_error(error, error_size,
                               "current Rusanov fixture direction coverage mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  if(number_response_cases == 0) {
    m1_thcm1_fixture_set_error(
        error, error_size,
        "current Rusanov fixture has no radiation_N output response");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  m1_thcm1_fixture_free(&collection);
  if(failures != 0) {
    if(error == NULL || error[0] == '\0')
      m1_thcm1_fixture_set_error(error, error_size,
                                 "current Rusanov fixture operand validation failed");
    return 0;
  }
  return 1;
}

static inline int m1_thcm1_rusanov_check_generic_fixture(
      const char *restrict fixture_dir,
      char *restrict error,
      const size_t error_size) {
  char *path = NULL;
  if(!m1_thcm1_rusanov_join_path(fixture_dir, "rusanov_generic.dat", &path)) {
    m1_thcm1_fixture_set_error(error, error_size,
                               "generic Rusanov fixture path allocation failed");
    return 0;
  }
  m1_thcm1_fixture_collection collection = {0};
  const int loaded = m1_thcm1_fixture_load(
      path, "rusanov_flux", M1_THCM1_RUSANOV_GENERIC_INPUT_COUNT,
      M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT, &collection, error, error_size);
  free(path);
  if(!loaded) return 0;
  if(strcmp(collection.policy, M1_THCM1_RUSANOV_POLICY) != 0 ||
     collection.record_count != M1_THCM1_RUSANOV_RECORD_COUNT) {
    m1_thcm1_fixture_set_error(error, error_size,
                               "generic Rusanov fixture policy or count mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  size_t failures = 0;
  for(size_t index = 0; index < collection.record_count; ++index) {
    const m1_thcm1_fixture_record *record = &collection.records[index];
    char expected_case_id[32];
    snprintf(expected_case_id, sizeof(expected_case_id), "face%04zu", index);
    if(strcmp(record->case_id, expected_case_id) != 0 ||
       !m1_thcm1_rusanov_pair_id_matches_operation(record, "rusanov_flux")) {
      m1_thcm1_fixture_set_error(error, error_size,
                                 "generic Rusanov fixture face ordering mismatch");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    double normalization[4] = {0.0};
    double output[4] = {0.0};
    m1_thcm1_rusanov_generic_normalization(record, normalization);
    const int status_ok = m1_thcm1_rusanov_eval_generic(record, 0, output);
    m1_thcm1_fixture_comparison_report report = {0};
    char compare_error[256] = {0};
    if(!status_ok || !m1_thcm1_rusanov_compare_baseline_response(
           record, collection.policy, "rusanov_flux", normalization, output,
           &report, compare_error, sizeof(compare_error))) {
      ++failures;
      if(failures == 1)
        m1_thcm1_fixture_set_error(
            error, error_size,
            compare_error[0] != '\0' ? compare_error :
                "generic Rusanov fixture evaluation failed");
    }
  }
  m1_thcm1_fixture_free(&collection);
  return failures == 0;
}

#endif /* UNIT_TESTS_M1_THCM1_RUSANOV_FIXTURE_H_ */
