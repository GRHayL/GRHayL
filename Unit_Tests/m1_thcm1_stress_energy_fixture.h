#ifndef UNIT_TESTS_M1_THCM1_STRESS_ENERGY_FIXTURE_H_
#define UNIT_TESTS_M1_THCM1_STRESS_ENERGY_FIXTURE_H_

/*
 * Replay the portable stress_energy pairs exported from the Verification/
 * CL-04 campaign.  The fixture contains the complete consumed radiation
 * state, ADM metric, and fluid velocity for each pair.  The current public
 * implementation is evaluated at the baseline input; the trusted perturbed
 * output is retained as the local numerical-response envelope, exactly as
 * in the other frozen THC_M1 replay lanes.
 */

#include "ghl_m1.h"
#include "m1_thcm1_fixture_utils.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#define M1_THCM1_STRESS_ENERGY_FIXTURE_OPERATION "stress_energy"
#define M1_THCM1_STRESS_ENERGY_FIXTURE_POLICY \
  "strict_relative_2e-12_propagated_response_v1"
#define M1_THCM1_STRESS_ENERGY_FIXTURE_INPUT_COUNT 21
#define M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT 10
/* This is the registry-backed CL-04 corpus size exported into this package. */
#define M1_THCM1_STRESS_ENERGY_FIXTURE_RECORD_COUNT 1024

enum {
  M1_THCM1_STRESS_ENERGY_N = 0,
  M1_THCM1_STRESS_ENERGY_E = 1,
  M1_THCM1_STRESS_ENERGY_F0 = 2,
  M1_THCM1_STRESS_ENERGY_F1 = 3,
  M1_THCM1_STRESS_ENERGY_F2 = 4,
  M1_THCM1_STRESS_ENERGY_LAPSE = 5,
  M1_THCM1_STRESS_ENERGY_SHIFT0 = 6,
  M1_THCM1_STRESS_ENERGY_SHIFT1 = 7,
  M1_THCM1_STRESS_ENERGY_SHIFT2 = 8,
  M1_THCM1_STRESS_ENERGY_GAMMA00 = 9,
  M1_THCM1_STRESS_ENERGY_V0 = 18
};

#define M1_THCM1_STRESS_ENERGY_EPSILON_C 1.0e-10
#define M1_THCM1_STRESS_ENERGY_ENERGY_FLOOR 1.0e-12
#define M1_THCM1_STRESS_ENERGY_ZETA_MIN 1.0e-8
#define M1_THCM1_STRESS_ENERGY_FD_EPSILON_REL 1.0e-6
#define M1_THCM1_STRESS_ENERGY_FD_EPSILON_ABS 1.0e-12
#define M1_THCM1_STRESS_ENERGY_NEWTON_MAX_ITERATIONS 20
#define M1_THCM1_STRESS_ENERGY_NEWTON_TOLERANCE 1.0e-10

static inline void m1_thcm1_stress_energy_fixture_error(
      char *restrict error,
      const size_t error_size,
      const char *restrict message) {
  if(error != NULL && error_size > 0)
    snprintf(error, error_size, "%s", message != NULL ? message :
             "stress-energy fixture error");
}

static inline int m1_thcm1_stress_energy_fixture_build_metric(
      const double *restrict input,
      ghl_metric_quantities *restrict metric) {
  if(input == NULL || metric == NULL) return 0;
  ghl_initialize_metric(
      input[M1_THCM1_STRESS_ENERGY_LAPSE],
      input[M1_THCM1_STRESS_ENERGY_SHIFT0],
      input[M1_THCM1_STRESS_ENERGY_SHIFT1],
      input[M1_THCM1_STRESS_ENERGY_SHIFT2],
      input[M1_THCM1_STRESS_ENERGY_GAMMA00 + 0],
      input[M1_THCM1_STRESS_ENERGY_GAMMA00 + 1],
      input[M1_THCM1_STRESS_ENERGY_GAMMA00 + 2],
      input[M1_THCM1_STRESS_ENERGY_GAMMA00 + 4],
      input[M1_THCM1_STRESS_ENERGY_GAMMA00 + 5],
      input[M1_THCM1_STRESS_ENERGY_GAMMA00 + 8], metric);
  return metric->lapse > 0.0 && isfinite(metric->lapse) &&
         metric->detgamma > 0.0 && isfinite(metric->detgamma) &&
         metric->sqrt_detgamma > 0.0 && isfinite(metric->sqrt_detgamma);
}

static inline int m1_thcm1_stress_energy_fixture_build_primitives(
      const double *restrict input,
      const ghl_metric_quantities *restrict metric,
      ghl_primitive_quantities *restrict prims) {
  if(input == NULL || metric == NULL || prims == NULL) return 0;
  *prims = (ghl_primitive_quantities){0};
  double v_con[3] = {0.0, 0.0, 0.0};
  for(int i = 0; i < 3; ++i)
    v_con[i] = (input[M1_THCM1_STRESS_ENERGY_V0 + i] + metric->betaU[i]) /
               metric->lapse;
  double v2 = 0.0;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      v2 += metric->gammaDD[i][j] * v_con[i] * v_con[j];
  if(!isfinite(v2) || v2 < 0.0 || v2 >= 1.0) return 0;
  for(int i = 0; i < 3; ++i)
    prims->vU[i] = input[M1_THCM1_STRESS_ENERGY_V0 + i];
  prims->u0 = 1.0 / (metric->lapse * sqrt(1.0 - v2));
  return isfinite(prims->u0);
}

static inline void m1_thcm1_stress_energy_fixture_four_metric(
      const ghl_metric_quantities *restrict metric,
      double g_dd[4][4]) {
  memset(g_dd, 0, sizeof(double[4][4]));
  double beta_d[3] = {0.0, 0.0, 0.0};
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      beta_d[i] += metric->gammaDD[i][j] * metric->betaU[j];
  g_dd[0][0] = -metric->lapse * metric->lapse;
  for(int i = 0; i < 3; ++i) {
    g_dd[0][0] += beta_d[i] * metric->betaU[i];
    g_dd[0][i+1] = beta_d[i];
    g_dd[i+1][0] = beta_d[i];
    for(int j = 0; j < 3; ++j)
      g_dd[i+1][j+1] = metric->gammaDD[i][j];
  }
}

static inline int m1_thcm1_stress_energy_fixture_evaluate(
      const double input[M1_THCM1_STRESS_ENERGY_FIXTURE_INPUT_COUNT],
      double output[M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT],
      char *restrict error,
      const size_t error_size) {
  if(input == NULL || output == NULL) return 0;
  ghl_metric_quantities metric = {0};
  ghl_primitive_quantities prims = {0};
  ghl_m1_parameters params = {0};
  ghl_m1_neutrino_state state = {
      .N = input[M1_THCM1_STRESS_ENERGY_N],
      .E = input[M1_THCM1_STRESS_ENERGY_E],
      .F = {input[M1_THCM1_STRESS_ENERGY_F0],
            input[M1_THCM1_STRESS_ENERGY_F1],
            input[M1_THCM1_STRESS_ENERGY_F2]}};
  ghl_m1_rad_state radiation = {
      .E = state.E,
      .F = {state.F[0], state.F[1], state.F[2]}};
  ghl_m1_closure closure = {0};
  ghl_stress_energy stress = {0};
  if(!m1_thcm1_stress_energy_fixture_build_metric(input, &metric) ||
     !m1_thcm1_stress_energy_fixture_build_primitives(input, &metric, &prims) ||
     ghl_m1_initialize(
         M1_THCM1_STRESS_ENERGY_EPSILON_C,
         M1_THCM1_STRESS_ENERGY_ENERGY_FLOOR,
         M1_THCM1_STRESS_ENERGY_ZETA_MIN,
         M1_THCM1_STRESS_ENERGY_FD_EPSILON_REL,
         M1_THCM1_STRESS_ENERGY_FD_EPSILON_ABS,
         M1_THCM1_STRESS_ENERGY_NEWTON_MAX_ITERATIONS,
         M1_THCM1_STRESS_ENERGY_NEWTON_TOLERANCE, &params) != ghl_success) {
    m1_thcm1_stress_energy_fixture_error(
        error, error_size, "stress-energy fixture initialization failed");
    return 0;
  }
  if(ghl_m1_compute_closure_with_primitives(
         &params, &metric, &prims, &radiation, &closure) != ghl_success ||
     ghl_m1_compute_neutrino_stress_energy(
         &params, &metric, &state, &closure, &stress) != ghl_success) {
    m1_thcm1_stress_energy_fixture_error(
        error, error_size, "stress-energy fixture public evaluation failed");
    return 0;
  }

  double g_dd[4][4];
  double lowered[4][4] = {{0.0}};
  m1_thcm1_stress_energy_fixture_four_metric(&metric, g_dd);
  for(int mu = 0; mu < 4; ++mu)
    for(int nu = 0; nu < 4; ++nu)
      for(int a = 0; a < 4; ++a)
        for(int b = 0; b < 4; ++b)
          lowered[mu][nu] += g_dd[mu][a] * g_dd[nu][b] * stress.T4[a][b];
  static const int components[M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT][2] = {
      {0,0}, {0,1}, {0,2}, {0,3}, {1,1},
      {1,2}, {1,3}, {2,2}, {2,3}, {3,3}};
  for(size_t component = 0;
      component < M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT; ++component) {
    output[component] = lowered[components[component][0]][components[component][1]];
    if(!isfinite(output[component])) {
      m1_thcm1_stress_energy_fixture_error(
          error, error_size, "stress-energy fixture produced a nonfinite output");
      return 0;
    }
  }
  return 1;
}

static inline int m1_thcm1_stress_energy_fixture_normalization(
      const double baseline[M1_THCM1_STRESS_ENERGY_FIXTURE_INPUT_COUNT],
      const double perturbed[M1_THCM1_STRESS_ENERGY_FIXTURE_INPUT_COUNT],
      double normalization[M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT]) {
  if(baseline == NULL || perturbed == NULL || normalization == NULL) return 0;
  const double scale = fmax(
      1.0e-14,
      fmax(fabs(baseline[M1_THCM1_STRESS_ENERGY_E]),
           fabs(perturbed[M1_THCM1_STRESS_ENERGY_E])));
  if(!isfinite(scale) || !(scale > 0.0)) return 0;
  for(size_t component = 0;
      component < M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT; ++component)
    normalization[component] = scale;
  return 1;
}

static inline int m1_thcm1_run_stress_energy_fixtures(
      const char *restrict fixture_path,
      size_t *restrict records_run,
      char *restrict error,
      const size_t error_size) {
  if(records_run != NULL) *records_run = 0;
  m1_thcm1_fixture_collection collection = {0};
  if(!m1_thcm1_fixture_load(
         fixture_path, M1_THCM1_STRESS_ENERGY_FIXTURE_OPERATION,
         M1_THCM1_STRESS_ENERGY_FIXTURE_INPUT_COUNT,
         M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT,
         &collection, error, error_size))
    return 0;
  if(collection.policy == NULL ||
     strcmp(collection.policy, M1_THCM1_STRESS_ENERGY_FIXTURE_POLICY) != 0 ||
     collection.record_count != M1_THCM1_STRESS_ENERGY_FIXTURE_RECORD_COUNT) {
    m1_thcm1_stress_energy_fixture_error(
        error, error_size, "stress-energy fixture inventory or policy mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }

  for(size_t index = 0; index < collection.record_count; ++index) {
    m1_thcm1_fixture_record *record = &collection.records[index];
    char expected_case[64];
    char expected_pair[72];
    char expected_seed[64];
    snprintf(expected_case, sizeof(expected_case), "stress-energy-%04zu", index);
    snprintf(expected_pair, sizeof(expected_pair), "%s-pair", expected_case);
    snprintf(expected_seed, sizeof(expected_seed), "stress-energy-seed-v1-%04zu", index);
    if(strcmp(record->case_id, expected_case) != 0 ||
       strcmp(record->pair_id, expected_pair) != 0 ||
       strcmp(record->origin, "external_thcm1_direct_grid") != 0 ||
       strcmp(record->seed_id, expected_seed) != 0 ||
       strcmp(record->family, "radiation") != 0 ||
       strcmp(record->perturbation, "radiation.E") != 0 ||
       record->sensitivity_start != M1_THCM1_STRESS_ENERGY_E ||
       record->sensitivity_count != 1) {
      m1_thcm1_stress_energy_fixture_error(
          error, error_size, "stress-energy fixture pair metadata mismatch");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    double normalization[M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT] = {0};
    if(!m1_thcm1_stress_energy_fixture_normalization(
           record->baseline_input, record->perturbed_input, normalization)) {
      m1_thcm1_stress_energy_fixture_error(
          error, error_size, "stress-energy normalization could not be recomputed");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    for(size_t component = 0;
        component < M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT; ++component)
      if(normalization[component] != record->normalization[component]) {
        m1_thcm1_stress_energy_fixture_error(
            error, error_size, "stress-energy normalization differs from its input");
        m1_thcm1_fixture_free(&collection);
        return 0;
      }
    double baseline[M1_THCM1_STRESS_ENERGY_FIXTURE_OUTPUT_COUNT] = {0};
    if(!m1_thcm1_stress_energy_fixture_evaluate(
           record->baseline_input, baseline, error, error_size)) {
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    m1_thcm1_fixture_comparison_report report = {0};
    if(!m1_thcm1_fixture_compare_baseline_response(
           record, collection.policy, normalization, baseline,
           &report, error, error_size)) {
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    if(records_run != NULL) ++(*records_run);
  }
  const size_t completed = records_run != NULL ? *records_run : collection.record_count;
  m1_thcm1_fixture_free(&collection);
  ghl_info("THC_M1 stress-energy fixtures: %zu trusted/perturbed agreement pairs passed\n",
           completed);
  return 1;
}

#endif /* UNIT_TESTS_M1_THCM1_STRESS_ENERGY_FIXTURE_H_ */
