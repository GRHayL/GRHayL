#ifndef UNIT_TESTS_M1_THCM1_SOURCE_FIXTURE_H_
#define UNIT_TESTS_M1_THCM1_SOURCE_FIXTURE_H_

/*
 * Test-local replay of the retained instantaneous frozen-rate source lane.
 * The common fixture reader owns parsing, duplicate-ID rejection,
 * finite/status validation, and the baseline/paired-response comparators.
 * This adapter
 * only maps the source operation's complete input vector to the public GRHayL
 * API and recomputes its input-defined normalization.  The retained
 * source_a1_a2_rate_normalized_v1 policy retains the trusted source baseline
 * and perturbed outputs for paired-response comparison.  Ordinary replay
 * evaluates the current API at both endpoints and compares the paired
 * response.
 * The six named zero-flux policy records remain
 * explicit local two-state checks.
 */

#include "ghl_m1.h"
#include "m1_thcm1_fixture_utils.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#define M1_THCM1_SOURCE_FIXTURE_OPERATION     "instantaneous_interaction_sources"
#define M1_THCM1_SOURCE_FIXTURE_POLICY        "source_a1_a2_rate_normalized_v1"
#define M1_THCM1_SOURCE_FIXTURE_INPUT_COUNT   40
#define M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT  9
/* Derived from the retained source_campaign.json corpus: 56 pairs times the
 * three retained species, not an acceptance limit invented by this adapter. */
#define M1_THCM1_SOURCE_FIXTURE_RECORD_COUNT  168

/* These are the exact controls passed by current_official/source_campaign.cc
 * to ghl_m1_initialize and the retained source packet contract. */
#define M1_THCM1_SOURCE_EPSILON_C             1.0e-8
#define M1_THCM1_SOURCE_ENERGY_FLOOR          1.0e-12
#define M1_THCM1_SOURCE_ZETA_MIN              1.0e-12
#define M1_THCM1_SOURCE_FD_EPSILON_REL        1.0e-6
#define M1_THCM1_SOURCE_FD_EPSILON_ABS        1.0e-8
#define M1_THCM1_SOURCE_NEWTON_MAX_ITERATIONS 20
#define M1_THCM1_SOURCE_NEWTON_TOLERANCE      1.0e-10
#define M1_THCM1_SOURCE_NUMBER_FLOOR          1.0e-14

/* Flattened source input layout.  The six pair-emissivity slots are zero in
 * the retained canonical sidecar, but remain explicit so the complete public
 * rates struct is represented rather than silently defaulted by the replay.
 * Baryon density, dt, and source thermal limit are endpoint/host controls and
 * are deliberately outside this instantaneous operation's input vector. */
enum {
  M1_THCM1_SOURCE_N = 0,
  M1_THCM1_SOURCE_E = 1,
  M1_THCM1_SOURCE_F0 = 2,
  M1_THCM1_SOURCE_F1 = 3,
  M1_THCM1_SOURCE_F2 = 4,
  M1_THCM1_SOURCE_LAPSE = 5,
  M1_THCM1_SOURCE_SHIFT0 = 6,
  M1_THCM1_SOURCE_SHIFT1 = 7,
  M1_THCM1_SOURCE_SHIFT2 = 8,
  M1_THCM1_SOURCE_GAMMA00 = 9,
  M1_THCM1_SOURCE_V0 = 18,
  M1_THCM1_SOURCE_ETA_N = 21,
  M1_THCM1_SOURCE_ETA_E = 22,
  M1_THCM1_SOURCE_KAPPA_A_N = 23,
  M1_THCM1_SOURCE_KAPPA_A_E = 24,
  M1_THCM1_SOURCE_KAPPA_S = 25,
  M1_THCM1_SOURCE_KAPPA_TR = 26,
  M1_THCM1_SOURCE_N_EQ = 27,
  M1_THCM1_SOURCE_J_EQ = 28,
  M1_THCM1_SOURCE_MEAN_ENERGY = 29,
  M1_THCM1_SOURCE_LEPTON_WEIGHT = 30,
  M1_THCM1_SOURCE_ETA_N_CC = 31,
  M1_THCM1_SOURCE_KAPPA_A_N_CC = 32,
  M1_THCM1_SOURCE_ETA_N_PAIR0 = 33,
  M1_THCM1_SOURCE_ETA_E_PAIR0 = 36,
  M1_THCM1_SOURCE_SQRT_DETGAMMA = 39
};

static inline void m1_thcm1_source_fixture_error(
      char *restrict error,
      const size_t error_size,
      const char *restrict message) {
  if(error != NULL && error_size > 0) {
    snprintf(
          error, error_size, "%s", message != NULL ? message : "source fixture error");
  }
}

static inline int m1_thcm1_source_fixture_species(const char *restrict name) {
  if(name == NULL) {
    return -1;
  }
  if(strcmp(name, "nue") == 0) {
    return ghl_m1_neutrino_nue;
  }
  if(strcmp(name, "anue") == 0) {
    return ghl_m1_neutrino_anue;
  }
  if(strcmp(name, "nux") == 0) {
    return ghl_m1_neutrino_nux;
  }
  return -1;
}

static inline int m1_thcm1_source_fixture_same_input(
      const double *restrict baseline,
      const double *restrict perturbed) {
  for(size_t i = 0; i < M1_THCM1_SOURCE_FIXTURE_INPUT_COUNT; ++i) {
    if(baseline[i] != perturbed[i]) {
      return 0;
    }
  }
  return 1;
}

static inline int m1_thcm1_source_fixture_build_metric(
      const double *restrict input,
      ghl_metric_quantities *restrict metric) {
  if(input == NULL || metric == NULL) {
    return 0;
  }
  ghl_initialize_metric(
        input[M1_THCM1_SOURCE_LAPSE], input[M1_THCM1_SOURCE_SHIFT0],
        input[M1_THCM1_SOURCE_SHIFT1], input[M1_THCM1_SOURCE_SHIFT2],
        input[M1_THCM1_SOURCE_GAMMA00 + 0 * 3 + 0],
        input[M1_THCM1_SOURCE_GAMMA00 + 0 * 3 + 1],
        input[M1_THCM1_SOURCE_GAMMA00 + 0 * 3 + 2],
        input[M1_THCM1_SOURCE_GAMMA00 + 1 * 3 + 1],
        input[M1_THCM1_SOURCE_GAMMA00 + 1 * 3 + 2],
        input[M1_THCM1_SOURCE_GAMMA00 + 2 * 3 + 2], metric);
  /* The source campaign consumed this retained sidecar volume.  Preserve it
   * exactly for the matter-coupling units; do not substitute a fixture value. */
  metric->sqrt_detgamma = input[M1_THCM1_SOURCE_SQRT_DETGAMMA];
  metric->detgamma = metric->sqrt_detgamma * metric->sqrt_detgamma;
  return isfinite(metric->sqrt_detgamma) && metric->sqrt_detgamma > 0.0;
}

static inline int m1_thcm1_source_fixture_build_rates(
      const double *restrict input,
      const int species,
      ghl_m1_neutrino_rates *restrict rates) {
  if(input == NULL || rates == NULL || species < 0 || species >= 3) {
    return 0;
  }
  *rates = (ghl_m1_neutrino_rates){ 0 };
  rates->species = (ghl_m1_neutrino_species_t)species;
  rates->eta_N = input[M1_THCM1_SOURCE_ETA_N];
  rates->eta_E = input[M1_THCM1_SOURCE_ETA_E];
  rates->kappa_a_N = input[M1_THCM1_SOURCE_KAPPA_A_N];
  rates->kappa_a_E = input[M1_THCM1_SOURCE_KAPPA_A_E];
  rates->kappa_s = input[M1_THCM1_SOURCE_KAPPA_S];
  rates->kappa_tr = input[M1_THCM1_SOURCE_KAPPA_TR];
  rates->n_eq = input[M1_THCM1_SOURCE_N_EQ];
  rates->J_eq = input[M1_THCM1_SOURCE_J_EQ];
  rates->mean_energy = input[M1_THCM1_SOURCE_MEAN_ENERGY];
  rates->lepton_weight = input[M1_THCM1_SOURCE_LEPTON_WEIGHT];
  rates->eta_N_cc = input[M1_THCM1_SOURCE_ETA_N_CC];
  rates->kappa_a_N_cc = input[M1_THCM1_SOURCE_KAPPA_A_N_CC];
  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    rates->eta_N_pair[process] = input[M1_THCM1_SOURCE_ETA_N_PAIR0 + process];
    rates->eta_E_pair[process] = input[M1_THCM1_SOURCE_ETA_E_PAIR0 + process];
  }
  return 1;
}

static inline int m1_thcm1_source_fixture_local_policy_case(const char *case_id) {
  static const char *const prefixes[]
        = { "rngpkt-v2-rd-radiation-anchor-energy-floor-a01:baseline:",
            "rngpkt-v2-rd-metric-anchor-offdiagonal-spd-a01:baseline:" };
  for(size_t i = 0; i < sizeof(prefixes) / sizeof(prefixes[0]); ++i) {
    const size_t length = strlen(prefixes[i]);
    if(strncmp(case_id, prefixes[i], length) == 0
       && m1_thcm1_source_fixture_species(case_id + length) >= 0) {
      return 1;
    }
  }
  return 0;
}

static inline int m1_thcm1_source_fixture_evaluate_one(
      const m1_thcm1_fixture_record *restrict record,
      const double *restrict input,
      const int species,
      double output[M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT],
      char *restrict error,
      const size_t error_size) {
  if(record == NULL || input == NULL || output == NULL || species < 0 || species >= 3) {
    return 0;
  }
  ghl_metric_quantities metric = { 0 };
  ghl_primitive_quantities prims = { 0 };
  ghl_m1_neutrino_state state
        = { .N = input[M1_THCM1_SOURCE_N],
            .E = input[M1_THCM1_SOURCE_E],
            .F = { input[M1_THCM1_SOURCE_F0], input[M1_THCM1_SOURCE_F1],
                   input[M1_THCM1_SOURCE_F2] } };
  ghl_m1_neutrino_rates rates = { 0 };
  ghl_m1_parameters m1_params = { 0 };
  ghl_m1_neutrino_parameters nu_params = { 0 };
  if(!m1_thcm1_source_fixture_build_metric(input, &metric)
     || !m1_thcm1_source_fixture_build_rates(input, species, &rates)
     || ghl_m1_initialize(
              M1_THCM1_SOURCE_EPSILON_C, M1_THCM1_SOURCE_ENERGY_FLOOR,
              M1_THCM1_SOURCE_ZETA_MIN, M1_THCM1_SOURCE_FD_EPSILON_REL,
              M1_THCM1_SOURCE_FD_EPSILON_ABS, M1_THCM1_SOURCE_NEWTON_MAX_ITERATIONS,
              M1_THCM1_SOURCE_NEWTON_TOLERANCE, &m1_params)
              != ghl_success) {
    m1_thcm1_source_fixture_error(
          error, error_size, "source fixture input initialization failed");
    return 0;
  }
  for(int i = 0; i < 3; ++i) {
    prims.vU[i] = input[M1_THCM1_SOURCE_V0 + i];
  }
  /* Match source_campaign.cc::make_packet exactly.  The current source API
   * recomputes W from vU, but the producer packet also materializes u0 and the
   * replay must verify that complete primitive requirement. */
  if(!(metric.lapse > 0.0) || !isfinite(metric.lapse)) {
    m1_thcm1_source_fixture_error(error, error_size, "source fixture lapse is invalid");
    return 0;
  }
  double v_con[3] = { 0.0, 0.0, 0.0 };
  double v2 = 0.0;
  for(int i = 0; i < 3; ++i) {
    v_con[i] = (prims.vU[i] + metric.betaU[i]) / metric.lapse;
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      v2 += metric.gammaDD[i][j] * v_con[i] * v_con[j];
    }
  }
  if(!isfinite(v2) || v2 < 0.0 || v2 >= 1.0) {
    m1_thcm1_source_fixture_error(
          error, error_size, "source fixture velocity is invalid");
    return 0;
  }
  const double W = 1.0 / sqrt(1.0 - v2);
  prims.u0 = W / metric.lapse;
  if(!isfinite(W) || !isfinite(prims.u0)) {
    m1_thcm1_source_fixture_error(error, error_size, "source fixture u0 is invalid");
    return 0;
  }
  nu_params.N_floor = M1_THCM1_SOURCE_NUMBER_FLOOR;
  nu_params.terminal_fallback_policy = ghl_m1_neutrino_terminal_fallback_no_update_all;

  if(m1_thcm1_source_fixture_local_policy_case(record->case_id)) {
    const ghl_m1_rad_state radiation
          = { .E = state.E, .F = { state.F[0], state.F[1], state.F[2] } };
    ghl_m1_closure closure = { 0 };
    if(state.F[0] != 0.0 || state.F[1] != 0.0 || state.F[2] != 0.0
       || ghl_m1_compute_closure_with_primitives(
                &m1_params, &metric, &prims, &radiation, &closure)
                != ghl_success
       || closure.solve_status != ghl_m1_closure_solve_endpoint_fallback
       || closure.four_point_compatibility) {
      m1_thcm1_source_fixture_error(
            error, error_size,
            "named source policy case no longer uses the zero-flux fallback");
      return 0;
    }
    long double trace = 0.0L;
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        trace += (long double)metric.gammaDD[i][j] * closure.P[i][j];
      }
    }
    if(!isfinite(trace)
       || fabsl(trace - state.E)
                > 128.0L * DBL_EPSILON * fmaxl(fabsl(trace), fabsl(state.E))
       || fabs(closure.chi - 1.0 / 3.0) > 128.0 * DBL_EPSILON) {
      m1_thcm1_source_fixture_error(
            error, error_size,
            "source zero-flux fallback violates its published-tensor contract");
      return 0;
    }
  }

  ghl_m1_sources ef_sources = { 0 };
  double number_source = NAN;
  const ghl_error_codes_t source_error = ghl_m1_compute_neutrino_interaction_sources(
        &m1_params, &nu_params, &metric, &prims, &state, &rates, &ef_sources,
        &number_source);
  if(source_error != ghl_success) {
    if(error != NULL && error_size > 0) {
      snprintf(
            error, error_size, "source fixture %s runtime source failed (%d)",
            record->case_id != NULL ? record->case_id : "<unknown>", (int)source_error);
    }
    return 0;
  }
  double matter_energy = NAN;
  double matter_momentum[3] = { NAN, NAN, NAN };
  const ghl_error_codes_t matter_error = ghl_m1_compute_neutrino_matter_coupling_sources(
        &metric, &ef_sources, &matter_energy, matter_momentum);
  if(matter_error != ghl_success) {
    if(error != NULL && error_size > 0) {
      snprintf(
            error, error_size, "source fixture %s runtime matter coupling failed (%d)",
            record->case_id != NULL ? record->case_id : "<unknown>", (int)matter_error);
    }
    return 0;
  }
  output[0] = number_source;
  output[1] = ef_sources.S_E;
  output[2] = ef_sources.S[0];
  output[3] = ef_sources.S[1];
  output[4] = ef_sources.S[2];
  output[5] = matter_energy;
  output[6] = matter_momentum[0];
  output[7] = matter_momentum[1];
  output[8] = matter_momentum[2];
  if(!m1_thcm1_fixture_finite_vector(output, M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT)) {
    return 0;
  }
  if(m1_thcm1_source_fixture_local_policy_case(record->case_id)) {
    const double volume = metric.lapse * metric.sqrt_detgamma;
    for(size_t i = 0; i < 4; ++i) {
      if(output[5 + i] != -volume * output[1 + i]) {
        m1_thcm1_source_fixture_error(
              error, error_size,
              "source policy case violates radiation/matter exchange balance");
        return 0;
      }
    }
  }
  return 1;
}

static inline int m1_thcm1_source_fixture_normalization(
      const double *restrict baseline,
      const double *restrict perturbed,
      double normalization[M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT]) {
  const double n_scale
        = fmax(fabs(baseline[M1_THCM1_SOURCE_N]), fabs(perturbed[M1_THCM1_SOURCE_N]));
  const double e_scale
        = fmax(fabs(baseline[M1_THCM1_SOURCE_E]), fabs(perturbed[M1_THCM1_SOURCE_E]));
  /* The exporter evaluates multiplication and addition separately in
   * binary64. Materialize those products so compiler FMA contraction cannot
   * change the fixture-metadata normalization by one rounding unit. */
  const volatile double baseline_n_absorption
        = baseline[M1_THCM1_SOURCE_KAPPA_A_N] * baseline[M1_THCM1_SOURCE_N];
  const volatile double perturbed_n_absorption
        = perturbed[M1_THCM1_SOURCE_KAPPA_A_N] * perturbed[M1_THCM1_SOURCE_N];
  const volatile double baseline_e_absorption
        = baseline[M1_THCM1_SOURCE_KAPPA_TR] * baseline[M1_THCM1_SOURCE_E];
  const volatile double perturbed_e_absorption
        = perturbed[M1_THCM1_SOURCE_KAPPA_TR] * perturbed[M1_THCM1_SOURCE_E];
  const double rate_n = fmax(
        fabs(baseline[M1_THCM1_SOURCE_ETA_N]) + baseline_n_absorption,
        fabs(perturbed[M1_THCM1_SOURCE_ETA_N]) + perturbed_n_absorption);
  const double rate_e = fmax(
        fabs(baseline[M1_THCM1_SOURCE_ETA_E]) + baseline_e_absorption,
        fabs(perturbed[M1_THCM1_SOURCE_ETA_E]) + perturbed_e_absorption);
  const double n_normalization = rate_n != 0.0 ? rate_n : n_scale;
  const double e_normalization = rate_e != 0.0 ? rate_e : e_scale;
  const double volume = fmax(
        baseline[M1_THCM1_SOURCE_LAPSE] * baseline[M1_THCM1_SOURCE_SQRT_DETGAMMA],
        perturbed[M1_THCM1_SOURCE_LAPSE] * perturbed[M1_THCM1_SOURCE_SQRT_DETGAMMA]);
  normalization[0] = n_normalization;
  for(size_t i = 1; i < 5; ++i) {
    normalization[i] = e_normalization;
  }
  for(size_t i = 5; i < 9; ++i) {
    normalization[i] = volume * e_normalization;
  }
  for(size_t i = 0; i < M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT; ++i) {
    if(!isfinite(normalization[i]) || !(normalization[i] > 0.0)) {
      return 0;
    }
  }
  return 1;
}

/* Ordinary records reach the current source API at both paired endpoints and
 * use the common baseline/perturbed/response comparator.  The six named
 * zero-flux policy cases have local admissibility/exchange checks, not
 * reference agreement.  records_run counts only paired-agreement records. */
static inline int m1_thcm1_run_instantaneous_source_fixtures(
      const char *restrict fixture_path,
      size_t *restrict records_run,
      char *restrict error,
      const size_t error_size) {
  if(records_run != NULL) {
    *records_run = 0;
  }
  m1_thcm1_fixture_collection collection = { 0 };
  if(!m1_thcm1_fixture_load(
           fixture_path, M1_THCM1_SOURCE_FIXTURE_OPERATION,
           M1_THCM1_SOURCE_FIXTURE_INPUT_COUNT, M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT,
           &collection, error, error_size)) {
    return 0;
  }

  if(collection.policy == NULL
     || strcmp(collection.policy, M1_THCM1_SOURCE_FIXTURE_POLICY) != 0
     || collection.record_count != M1_THCM1_SOURCE_FIXTURE_RECORD_COUNT) {
    m1_thcm1_source_fixture_error(
          error, error_size,
          "instantaneous source fixture policy or retained record count mismatch");
    m1_thcm1_fixture_free(&collection);
    return 0;
  }

  size_t local_policy_cases = 0;
  for(size_t index = 0; index < collection.record_count; ++index) {
    m1_thcm1_fixture_record *record = &collection.records[index];
    const int species = m1_thcm1_source_fixture_species(record->seed_id);
    char expected_pair_id[M1_THCM1_FIXTURE_TEXT_MAX] = { 0 };
    snprintf(
          expected_pair_id, sizeof(expected_pair_id), "pair%zu:%s", index / 3,
          record->seed_id != NULL ? record->seed_id : "<missing>");
    if(species < 0 || strcmp(record->origin, "source_campaign") != 0
       || record->input_count != M1_THCM1_SOURCE_FIXTURE_INPUT_COUNT
       || record->output_count != M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT
       || strcmp(record->pair_id, expected_pair_id) != 0
       || (record->sensitivity_count == 0
           && !m1_thcm1_source_fixture_same_input(
                 record->baseline_input, record->perturbed_input))) {
      m1_thcm1_source_fixture_error(
            error, error_size, "invalid instantaneous source fixture metadata");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    double normalization[M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT] = { 0 };
    if(!m1_thcm1_source_fixture_normalization(
             record->baseline_input, record->perturbed_input, normalization)) {
      m1_thcm1_source_fixture_error(
            error, error_size, "instantaneous source normalization failed");
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    for(size_t component = 0; component < M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT;
        ++component) {
      if(normalization[component] != record->normalization[component]) {
        m1_thcm1_source_fixture_error(
              error, error_size,
              "source normalization differs from its consumed inputs");
        m1_thcm1_fixture_free(&collection);
        return 0;
      }
    }
    double baseline[M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT] = { 0 };
    if(!m1_thcm1_source_fixture_evaluate_one(
             record, record->baseline_input, species, baseline, error, error_size)
       || !m1_thcm1_fixture_finite_vector(
             baseline, M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT)) {
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    m1_thcm1_fixture_comparison_report report = { 0 };
    if(m1_thcm1_source_fixture_local_policy_case(record->case_id)) {
      double perturbed[M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT] = { 0 };
      if(!m1_thcm1_source_fixture_evaluate_one(
               record, record->perturbed_input, species, perturbed, error, error_size)
         || !m1_thcm1_fixture_finite_vector(
               perturbed, M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT)) {
        m1_thcm1_fixture_free(&collection);
        return 0;
      }
      ++local_policy_cases;
      /* These named records are local two-state admissibility/exchange
       * checks.  They are deliberately excluded from source agreement counts
       * because their published closure policy differs from THC_M1. */
      continue;
    }
    double perturbed[M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT] = { 0 };
    if(!m1_thcm1_source_fixture_evaluate_one(
             record, record->perturbed_input, species, perturbed, error, error_size)
       || !m1_thcm1_fixture_finite_vector(
             perturbed, M1_THCM1_SOURCE_FIXTURE_OUTPUT_COUNT)) {
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    if(!m1_thcm1_fixture_compare_paired(
             record, normalization, baseline, perturbed, &report, error, error_size)) {
      m1_thcm1_fixture_free(&collection);
      return 0;
    }
    if(records_run != NULL) {
      ++(*records_run);
    }
  }
  m1_thcm1_fixture_free(&collection);
  if(local_policy_cases != 6) {
    m1_thcm1_source_fixture_error(
          error, error_size,
          "source corpus changed its named zero-flux policy inventory");
    return 0;
  }
  ghl_info(
        "THC_M1 source fixtures: %zu named zero-flux policy pairs passed local "
        "admissibility/exchange checks\n",
        local_policy_cases);
  return 1;
}

#endif /* UNIT_TESTS_M1_THCM1_SOURCE_FIXTURE_H_ */
