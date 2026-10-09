#include "m1_thcm1_rusanov_fixture.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int m1_thcm1_rusanov_eval_neutrino(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const m1_thcm1_fixture_record *restrict record,
      const int role,
      double output[M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT]) {
  if(m1_params == NULL || nu_params == NULL || record == NULL
     || record->input_count != M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT
     || record->output_count != M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT
     || output == NULL) {
    return 0;
  }
  if(role < 0 || role > 1) {
    return 0;
  }
  const double *input = role == 0 ? record->baseline_input : record->perturbed_input;
  if(input == NULL || !isfinite(input[0]) || input[0] < 0.0 || input[0] > 2.0
     || floor(input[0]) != input[0]) {
    return 0;
  }
  ghl_metric_quantities metric;
  m1_thcm1_rusanov_metric(input, &metric);
  ghl_m1_neutrino_state state_L, state_R;
  m1_thcm1_rusanov_state(&input[13], &state_L);
  m1_thcm1_rusanov_state(&input[18], &state_R);
  ghl_m1_closure closure_L, closure_R;
  m1_thcm1_rusanov_closure(&input[23], &closure_L);
  m1_thcm1_rusanov_closure(&input[32], &closure_R);
  double flux_F[3] = { 0.0, 0.0, 0.0 };
  const ghl_error_codes_t status = ghl_m1_compute_neutrino_rusanov_flux(
        m1_params, nu_params, &metric, (ghl_m1_direction_t)(int)input[0], &state_L,
        &state_R, &closure_L, &closure_R, &input[41], &input[44], &input[47], &input[50],
        input[53], &output[0], &output[1], flux_F);
  if(status != ghl_success) {
    return 0;
  }
  output[2] = flux_F[0];
  output[3] = flux_F[1];
  output[4] = flux_F[2];
  return 1;
}

static void m1_thcm1_rusanov_neutrino_normalization(
      const m1_thcm1_fixture_record *restrict record,
      double normalization[M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT]) {
  for(int component = 0; component < 5; ++component) {
    normalization[component] = 0.0;
  }
  for(int role = 0; role < 2; ++role) {
    const double *input = role == 0 ? record->baseline_input : record->perturbed_input;
    for(int component = 0; component < 5; ++component) {
      const double state_L = input[13 + component];
      const double state_R = input[18 + component];
      const double scale = fmax(
            fmax(fabs(state_L), fabs(state_R)),
            fmax(input[53] * fabs(state_L), input[53] * fabs(state_R)));
      normalization[component]
            = fmax(normalization[component], fmax(1.0e-300, input[12] * scale));
    }
  }
}

static int m1_thcm1_rusanov_eval_generic(
      const m1_thcm1_fixture_record *restrict record,
      const int role,
      double output[M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT]) {
  if(record == NULL || record->input_count != M1_THCM1_RUSANOV_GENERIC_INPUT_COUNT
     || record->output_count != M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT
     || output == NULL) {
    return 0;
  }
  if(role < 0 || role > 1) {
    return 0;
  }
  const double *input = role == 0 ? record->baseline_input : record->perturbed_input;
  if(input == NULL) {
    return 0;
  }
  const ghl_m1_rad_state state_L
        = { .E = input[0], .F = { input[1], input[2], input[3] } };
  const ghl_m1_rad_state state_R
        = { .E = input[4], .F = { input[5], input[6], input[7] } };
  return ghl_m1_compute_rusanov_flux(
               &state_L, &state_R, input[8], &input[9], input[12], &input[13], input[16],
               &output[0], &output[1])
         == ghl_success;
}

static void m1_thcm1_rusanov_generic_normalization(
      const m1_thcm1_fixture_record *restrict record,
      double normalization[M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT]) {
  for(int component = 0; component < 4; ++component) {
    normalization[component] = 0.0;
  }
  for(int role = 0; role < 2; ++role) {
    const double *input = role == 0 ? record->baseline_input : record->perturbed_input;
    for(int component = 0; component < 4; ++component) {
      const double state_L = input[component];
      const double state_R = input[4 + component];
      const double scale = fmax(
            fmax(fabs(state_L), fabs(state_R)),
            fmax(input[16] * fabs(state_L), input[16] * fabs(state_R)));
      normalization[component] = fmax(normalization[component], fmax(1.0e-300, scale));
    }
  }
}

/* Evaluate both endpoint roles and run the mandatory baseline-envelope gate
 * plus the strict paired gate for the operation.  Only the two supported
 * Rusanov policies reach this point, so the retained corpus always pays for
 * the strict endpoint/response comparisons after the envelope. */
static int m1_thcm1_rusanov_eval_and_compare(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const int is_neutrino,
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict policy,
      const char *restrict operation,
      const double *restrict normalization,
      char *restrict error,
      const size_t error_size) {
  double output_baseline[M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT] = { 0.0 };
  double output_perturbed[M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT] = { 0.0 };
  if(is_neutrino) {
    if(!m1_thcm1_rusanov_eval_neutrino(m1_params, nu_params, record, 0, output_baseline)
       || !m1_thcm1_rusanov_eval_neutrino(
             m1_params, nu_params, record, 1, output_perturbed)) {
      m1_thcm1_fixture_set_error(
            error, error_size, "Rusanov fixture evaluator rejected a valid record");
      return 0;
    }
  }
  else {
    if(!m1_thcm1_rusanov_eval_generic(record, 0, output_baseline)
       || !m1_thcm1_rusanov_eval_generic(record, 1, output_perturbed)) {
      m1_thcm1_fixture_set_error(
            error, error_size, "Rusanov fixture evaluator rejected a valid record");
      return 0;
    }
  }
  m1_thcm1_fixture_comparison_report report = { 0 };
  if(!m1_thcm1_rusanov_compare_baseline_response(
        record, policy, operation, normalization, output_baseline, &report, error,
        error_size)
     || !m1_thcm1_rusanov_compare_paired_strict_relative(
           record, policy, operation, normalization, output_baseline, output_perturbed,
           &report, error, error_size)) {
    return 0;
  }
  return 1;
}

static int m1_thcm1_rusanov_join_path(
      const char *restrict directory,
      const char *restrict name,
      char **restrict path) {
  if(directory == NULL || name == NULL || path == NULL) {
    return 0;
  }
  const size_t directory_length = strlen(directory);
  const size_t name_length = strlen(name);
  *path = (char *)malloc(directory_length + name_length + 2);
  if(*path == NULL) {
    return 0;
  }
  memcpy(*path, directory, directory_length);
  (*path)[directory_length] = '/';
  memcpy(*path + directory_length + 1, name, name_length + 1);
  return 1;
}

int m1_thcm1_rusanov_check_neutrino_fixture(
      const char *restrict fixture_dir,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      char *restrict error,
      const size_t error_size) {
  char *path = NULL;
  if(!m1_thcm1_rusanov_join_path(fixture_dir, "rusanov_neutrino.bin", &path)) {
    m1_thcm1_fixture_set_error(
          error, error_size, "neutrino Rusanov fixture path allocation failed");
    return 0;
  }
  m1_thcm1_fixture_collection collection = { 0 };
  const int loaded = m1_thcm1_fixture_load(
        path, "neutrino_rusanov_flux", M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT,
        M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT, &collection, error, error_size);
  free(path);
  if(!loaded) {
    return 0;
  }
  if(strcmp(collection.policy, M1_THCM1_RUSANOV_POLICY) != 0
     || collection.record_count != M1_THCM1_RUSANOV_RECORD_COUNT) {
    m1_thcm1_fixture_set_error(
          error, error_size, "neutrino Rusanov fixture policy or count mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  size_t failures = 0;
  size_t direction_counts[3] = { 0, 0, 0 };
  for(size_t index = 0; index < collection.record_count; ++index) {
    const m1_thcm1_fixture_record *record = &collection.records[index];
    char expected_case_id[32];
    snprintf(expected_case_id, sizeof(expected_case_id), "face%04zu", index);
    if(strcmp(record->case_id, expected_case_id) != 0
       || !m1_thcm1_rusanov_pair_id_matches_operation(record, "neutrino_rusanov_flux")
       || !isfinite(record->baseline_input[0])
       || floor(record->baseline_input[0]) != record->baseline_input[0]
       || record->baseline_input[0] < 0.0 || record->baseline_input[0] > 2.0) {
      m1_thcm1_fixture_set_error(
            error, error_size, "neutrino Rusanov fixture face ordering mismatch");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    ++direction_counts[(int)record->baseline_input[0]];
    double normalization[5] = { 0.0 };
    m1_thcm1_rusanov_neutrino_normalization(record, normalization);
    char compare_error[256] = { 0 };
    if(!m1_thcm1_rusanov_eval_and_compare(
          m1_params, nu_params, 1, record, collection.policy,
          "neutrino_rusanov_flux", normalization, compare_error,
          sizeof(compare_error))) {
      ++failures;
      if(failures == 1) {
        m1_thcm1_fixture_set_error(
              error, error_size,
              compare_error[0] != '\0' ? compare_error
                                       : "neutrino Rusanov fixture evaluation failed");
      }
    }
  }
  if(direction_counts[0] != 342 || direction_counts[1] != 341
     || direction_counts[2] != 341) {
    m1_thcm1_fixture_set_error(
          error, error_size, "neutrino Rusanov fixture direction coverage mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  m1_thcm1_fixture_free(&collection);
  if(failures != 0) {
    return 0;
  }
  return 1;
}

int m1_thcm1_rusanov_check_current_fixture(
      const char *restrict fixture_dir,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      char *restrict error,
      const size_t error_size) {
  char *path = NULL;
  if(!m1_thcm1_rusanov_join_path(fixture_dir, "rusanov_neutrino_current.bin", &path)) {
    m1_thcm1_fixture_set_error(
          error, error_size, "current Rusanov fixture path allocation failed");
    return 0;
  }
  m1_thcm1_fixture_collection collection = { 0 };
  const int loaded = m1_thcm1_fixture_load(
        path, M1_THCM1_RUSANOV_CURRENT_OPERATION, M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT,
        M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT, &collection, error, error_size);
  free(path);
  if(!loaded) {
    return 0;
  }
  if(strcmp(collection.policy, M1_THCM1_RUSANOV_CURRENT_POLICY) != 0
     || collection.record_count != M1_THCM1_RUSANOV_RECORD_COUNT) {
    m1_thcm1_fixture_set_error(
          error, error_size, "current Rusanov fixture policy or count mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  size_t failures = 0;
  size_t direction_counts[3] = { 0, 0, 0 };
  size_t number_response_cases = 0;
  for(size_t index = 0; index < collection.record_count; ++index) {
    const m1_thcm1_fixture_record *record = &collection.records[index];
    char expected_case_id[32];
    snprintf(expected_case_id, sizeof(expected_case_id), "face%04zu", index);
    if(strcmp(record->case_id, expected_case_id) != 0
       || !m1_thcm1_rusanov_pair_id_matches_operation(
             record, M1_THCM1_RUSANOV_CURRENT_OPERATION)
       || !isfinite(record->baseline_input[0])
       || floor(record->baseline_input[0]) != record->baseline_input[0]
       || record->baseline_input[0] < 0.0 || record->baseline_input[0] > 2.0) {
      m1_thcm1_fixture_set_error(
            error, error_size, "current Rusanov fixture face ordering mismatch");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    ++direction_counts[(int)record->baseline_input[0]];
    for(int role = 0; role < 2; ++role) {
      const double *input = role == 0 ? record->baseline_input : record->perturbed_input;
      for(int component = 0; component < 12; ++component) {
        if(!isfinite(input[41 + component]) || input[41 + component] == 0.0) {
          ++failures;
        }
      }
      for(int component = 0; component < 3; ++component) {
        if(input[41 + component] != input[13] * input[47 + component]
           || input[44 + component] != input[18] * input[50 + component]) {
          ++failures;
        }
      }
    }
    if(strcmp(record->family, "radiation_N") == 0
       && record->baseline_output[0] != record->perturbed_output[0]) {
      ++number_response_cases;
    }
    double normalization[5] = { 0.0 };
    m1_thcm1_rusanov_neutrino_normalization(record, normalization);
    char compare_error[256] = { 0 };
    if(!m1_thcm1_rusanov_eval_and_compare(
          m1_params, nu_params, 1, record, collection.policy,
          M1_THCM1_RUSANOV_CURRENT_OPERATION, normalization, compare_error,
          sizeof(compare_error))) {
      ++failures;
      if(failures == 1) {
        m1_thcm1_fixture_set_error(
              error, error_size,
              compare_error[0] != '\0' ? compare_error
                                       : "current Rusanov fixture evaluation failed");
      }
    }
  }
  if(direction_counts[0] != 342 || direction_counts[1] != 341
     || direction_counts[2] != 341) {
    m1_thcm1_fixture_set_error(
          error, error_size, "current Rusanov fixture direction coverage mismatch");
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
    if(error == NULL || error[0] == '\0') {
      m1_thcm1_fixture_set_error(
            error, error_size, "current Rusanov fixture operand validation failed");
    }
    return 0;
  }
  return 1;
}

int m1_thcm1_rusanov_check_generic_fixture(
      const char *restrict fixture_dir,
      char *restrict error,
      const size_t error_size) {
  char *path = NULL;
  if(!m1_thcm1_rusanov_join_path(fixture_dir, "rusanov_generic.bin", &path)) {
    m1_thcm1_fixture_set_error(
          error, error_size, "generic Rusanov fixture path allocation failed");
    return 0;
  }
  m1_thcm1_fixture_collection collection = { 0 };
  const int loaded = m1_thcm1_fixture_load(
        path, "rusanov_flux", M1_THCM1_RUSANOV_GENERIC_INPUT_COUNT,
        M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT, &collection, error, error_size);
  free(path);
  if(!loaded) {
    return 0;
  }
  if(strcmp(collection.policy, M1_THCM1_RUSANOV_POLICY) != 0
     || collection.record_count != M1_THCM1_RUSANOV_RECORD_COUNT) {
    m1_thcm1_fixture_set_error(
          error, error_size, "generic Rusanov fixture policy or count mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }
  size_t failures = 0;
  for(size_t index = 0; index < collection.record_count; ++index) {
    const m1_thcm1_fixture_record *record = &collection.records[index];
    char expected_case_id[32];
    snprintf(expected_case_id, sizeof(expected_case_id), "face%04zu", index);
    if(strcmp(record->case_id, expected_case_id) != 0
       || !m1_thcm1_rusanov_pair_id_matches_operation(record, "rusanov_flux")) {
      m1_thcm1_fixture_set_error(
            error, error_size, "generic Rusanov fixture face ordering mismatch");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    double normalization[4] = { 0.0 };
    m1_thcm1_rusanov_generic_normalization(record, normalization);
    char compare_error[256] = { 0 };
    if(!m1_thcm1_rusanov_eval_and_compare(
          NULL, NULL, 0, record, collection.policy, "rusanov_flux", normalization,
          compare_error, sizeof(compare_error))) {
      ++failures;
      if(failures == 1) {
        m1_thcm1_fixture_set_error(
              error, error_size,
              compare_error[0] != '\0' ? compare_error
                                       : "generic Rusanov fixture evaluation failed");
      }
    }
  }
  m1_thcm1_fixture_free(&collection);
  return failures == 0;
}

/* Corrupt only stack-owned copies; the retained package stays untouched.
 * Return success only when the intended negative assertion rejects its copy. */
int m1_thcm1_rusanov_check_fixture_regression(
      const char *restrict fixture_dir,
      const int mutation_mode,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      char *restrict error,
      const size_t error_size) {
  const int current = mutation_mode == M1_THCM1_RUSANOV_REGRESSION_CORRUPTED_CURRENT;
  if(!current && mutation_mode != M1_THCM1_RUSANOV_REGRESSION_PERTURBED_OUTPUT
     && mutation_mode != M1_THCM1_RUSANOV_REGRESSION_SWAPPED_ENDPOINTS) {
    m1_thcm1_fixture_set_error(error, error_size, "invalid Rusanov regression mode");
    return 0;
  }
  const char *operation = current ? M1_THCM1_RUSANOV_CURRENT_OPERATION : "rusanov_flux";
  const char *policy = current ? M1_THCM1_RUSANOV_CURRENT_POLICY : M1_THCM1_RUSANOV_POLICY;
  char *path = NULL;
  if(!m1_thcm1_rusanov_join_path(
           fixture_dir, current ? "rusanov_neutrino_current.bin" : "rusanov_generic.bin",
           &path)) {
    m1_thcm1_fixture_set_error(error, error_size, "Rusanov regression path allocation failed");
    return 0;
  }
  m1_thcm1_fixture_collection collection = { 0 };
  const int loaded = m1_thcm1_fixture_load(
        path, operation,
        current ? M1_THCM1_RUSANOV_NEUTRINO_INPUT_COUNT : M1_THCM1_RUSANOV_GENERIC_INPUT_COUNT,
        current ? M1_THCM1_RUSANOV_NEUTRINO_OUTPUT_COUNT : M1_THCM1_RUSANOV_GENERIC_OUTPUT_COUNT,
        &collection, error, error_size);
  free(path);
  if(!loaded) {
    return 0;
  }
  int rejected = 0;
  for(size_t i = 0; i < collection.record_count; ++i) {
    m1_thcm1_fixture_record record = collection.records[i];
    if(record.baseline_output[0] == record.perturbed_output[0]) {
      continue;
    }
    double normalization[5] = { 0.0 };
    double perturbed_input[54];
    double perturbed_output[5];
    if(current) {
      m1_thcm1_rusanov_neutrino_normalization(&record, normalization);
      memcpy(perturbed_input, record.perturbed_input, sizeof(perturbed_input));
      perturbed_input[41] *= 2.0;
      record.perturbed_input = perturbed_input;
    }
    else {
      m1_thcm1_rusanov_generic_normalization(&record, normalization);
      if(mutation_mode == M1_THCM1_RUSANOV_REGRESSION_PERTURBED_OUTPUT) {
        memcpy(perturbed_output, record.perturbed_output, 4 * sizeof(double));
        perturbed_output[0] += fmax(1.0, fabs(perturbed_output[0]));
        record.perturbed_output = perturbed_output;
      }
      else {
        double *baseline_output = record.baseline_output;
        record.baseline_output = record.perturbed_output;
        record.perturbed_output = baseline_output;
      }
    }
    char rejection[256] = { 0 };
    rejected = !m1_thcm1_rusanov_eval_and_compare(
          m1_params, nu_params, current, &record, policy, operation, normalization,
          rejection, sizeof(rejection));
    if(current) {
      rejected &= strstr(rejection, "evaluator rejected") != NULL;
    }
    else if(mutation_mode == M1_THCM1_RUSANOV_REGRESSION_PERTURBED_OUTPUT) {
      rejected &= strstr(rejection, "perturbed comparison failed") != NULL;
    }
    break;
  }
  m1_thcm1_fixture_free(&collection);
  if(!rejected) {
    m1_thcm1_fixture_set_error(error, error_size, "Rusanov regression corruption was not rejected");
  }
  return rejected;
}
