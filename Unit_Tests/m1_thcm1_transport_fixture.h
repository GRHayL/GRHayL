#ifndef UNIT_TESTS_M1_THCM1_TRANSPORT_FIXTURE_H_
#define UNIT_TESTS_M1_THCM1_TRANSPORT_FIXTURE_H_

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "m1_thcm1_fixture_utils.h"

/* The transport fixture wire format is intentionally shared by the
 * constant-volume and prepared variable-volume families.  The first six
 * fields are controls, fields 6..9 are the four stencil-cell volumes, fields
 * 10..29 are four five-component states, and fields 30..39 are the two
 * five-component physical face fluxes. */
#define M1_THCM1_TRANSPORT_INPUT_COUNT           50
#define M1_THCM1_TRANSPORT_OUTPUT_COUNT          5
#define M1_THCM1_TRANSPORT_FACE_VOLUME_INDEX     0
#define M1_THCM1_TRANSPORT_SPEED_INDEX           1
#define M1_THCM1_TRANSPORT_KAPPA_INDEX           2
#define M1_THCM1_TRANSPORT_DELTA_INDEX           3
#define M1_THCM1_TRANSPORT_THETA_INDEX           4
#define M1_THCM1_TRANSPORT_MINDISS_INDEX         5
#define M1_THCM1_TRANSPORT_CELL_VOLUME_START     6
#define M1_THCM1_TRANSPORT_STATE_START           10
#define M1_THCM1_TRANSPORT_FLUX_L_START          30
#define M1_THCM1_TRANSPORT_FLUX_R_START          35
#define M1_THCM1_TRANSPORT_PRODUCER_FLUX_L_START 40
#define M1_THCM1_TRANSPORT_PRODUCER_FLUX_R_START 45

/*
 * Host-boundary preparation for the THC variable-volume convention.
 *
 * This adapter deliberately stages into local candidates and publishes only
 * after every volume and operand has been checked.  Thus a bad volume cannot
 * leave a partially weighted stencil in a caller-owned output buffer.
 */
static inline int m1_thcm1_transport_prepare_volume_weighted(
      const double cell_volumes[4],
      const double face_volume,
      const double state_stencil[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      const double physical_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      const double physical_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      double weighted_state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      double weighted_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      double weighted_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      char *restrict error,
      const size_t error_size) {
  if(cell_volumes == NULL || state_stencil == NULL || physical_flux_L == NULL
     || physical_flux_R == NULL || weighted_state == NULL || weighted_flux_L == NULL
     || weighted_flux_R == NULL) {
    m1_thcm1_fixture_set_error(error, error_size, "null volume-preparation argument");
    return 0;
  }
  if(!isfinite(face_volume) || !(face_volume > 0.0)) {
    m1_thcm1_fixture_set_error(
          error, error_size, "face volume must be finite and positive");
    return 0;
  }
  for(int cell = 0; cell < 4; ++cell) {
    if(!isfinite(cell_volumes[cell]) || !(cell_volumes[cell] > 0.0)) {
      m1_thcm1_fixture_set_error(
            error, error_size, "cell volume must be finite and positive");
      return 0;
    }
    for(int component = 0; component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
      if(!isfinite(state_stencil[cell][component])) {
        m1_thcm1_fixture_set_error(error, error_size, "transport state must be finite");
        return 0;
      }
    }
  }
  for(int component = 0; component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
    if(!isfinite(physical_flux_L[component]) || !isfinite(physical_flux_R[component])) {
      m1_thcm1_fixture_set_error(
            error, error_size, "transport physical flux must be finite");
      return 0;
    }
  }

  double candidate_state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double candidate_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double candidate_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  for(int cell = 0; cell < 4; ++cell) {
    for(int component = 0; component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
      candidate_state[cell][component]
            = cell_volumes[cell] * state_stencil[cell][component];
    }
  }
  for(int cell = 0; cell < 4; ++cell) {
    for(int component = 0; component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
      if(!isfinite(candidate_state[cell][component])) {
        m1_thcm1_fixture_set_error(
              error, error_size, "weighted stencil state is not finite");
        return 0;
      }
    }
  }
  for(int component = 0; component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
    candidate_flux_L[component] = face_volume * physical_flux_L[component];
    candidate_flux_R[component] = face_volume * physical_flux_R[component];
    if(!isfinite(candidate_flux_L[component])
       || !isfinite(candidate_flux_R[component])) {
      m1_thcm1_fixture_set_error(
            error, error_size, "weighted physical flux is not finite");
      return 0;
    }
  }
  memcpy(weighted_state, candidate_state, sizeof(candidate_state));
  memcpy(weighted_flux_L, candidate_flux_L, sizeof(candidate_flux_L));
  memcpy(weighted_flux_R, candidate_flux_R, sizeof(candidate_flux_R));
  return 1;
}

/* Match the existing fixture producer's scale convention after the same
 * volume preparation: face physical operands and speed times the first two
 * stencil rows.  The stored normalization is metadata, not a tolerance. */
static inline int m1_thcm1_transport_compute_weighted_normalization(
      const double input[M1_THCM1_TRANSPORT_INPUT_COUNT],
      double normalization[M1_THCM1_TRANSPORT_OUTPUT_COUNT]) {
  if(input == NULL || normalization == NULL) {
    return 0;
  }
  const double face_volume = input[M1_THCM1_TRANSPORT_FACE_VOLUME_INDEX];
  const double speed = fabs(input[M1_THCM1_TRANSPORT_SPEED_INDEX]);
  if(!isfinite(face_volume) || !(face_volume > 0.0) || !isfinite(speed)) {
    return 0;
  }
  for(int cell = 0; cell < 4; ++cell) {
    const double volume = input[M1_THCM1_TRANSPORT_CELL_VOLUME_START + cell];
    if(!isfinite(volume) || !(volume > 0.0)) {
      return 0;
    }
  }
  for(int component = 0; component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
    const int state_L_start = M1_THCM1_TRANSPORT_STATE_START;
    const int state_C_start
          = M1_THCM1_TRANSPORT_STATE_START + M1_THCM1_TRANSPORT_OUTPUT_COUNT;
    const double flux_L = input[M1_THCM1_TRANSPORT_FLUX_L_START + component];
    const double flux_R = input[M1_THCM1_TRANSPORT_FLUX_R_START + component];
    const double state_L = input[state_L_start + component];
    const double state_C = input[state_C_start + component];
    if(!isfinite(flux_L) || !isfinite(flux_R) || !isfinite(state_L)
       || !isfinite(state_C)) {
      return 0;
    }
    const double scale = fmax(
          fmax(fabs(face_volume * flux_L), fabs(face_volume * flux_R)),
          fmax(speed * fabs(input[M1_THCM1_TRANSPORT_CELL_VOLUME_START] * state_L),
               speed * fabs(input[M1_THCM1_TRANSPORT_CELL_VOLUME_START + 1] * state_C)));
    normalization[component] = fmax(1.0e-300, scale);
  }
  return 1;
}

/* Explicit two-state transport comparator.  Stored transport replay and
 * focused comparator tests use this stronger
 * current-baseline/current-perturbed campaign contract.  The
 * transport producer uses abs(a-b)/max(1e-300,abs(a),abs(b)) <= 2e-12.
 * Calling the shared scalar primitive with that denominator, absolute
 * tolerance 2e-12, zero relative tolerance, and zero floor preserves the
 * rule including sub-floor and exact-zero cases.
 * The fixture's input-derived normalization is still checked separately by
 * the transport consumer; it is not used as a tolerance knob here. */
static inline int m1_thcm1_transport_compare_paired(
      const m1_thcm1_fixture_record *restrict record,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      const double *restrict computed_perturbed,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size) {
  if(record == NULL || computed_normalization == NULL || computed_baseline == NULL
     || computed_perturbed == NULL || record->normalization == NULL
     || record->baseline_output == NULL || record->perturbed_output == NULL
     || record->output_count == 0) {
    m1_thcm1_fixture_set_error(
          error, error_size, "invalid transport paired comparison argument");
    return 0;
  }
  if(report != NULL) {
    *report = (m1_thcm1_fixture_comparison_report){ .max_scaled_error = 0.0,
                                                    .component = 0,
                                                    .gate = "baseline" };
  }
  const double *expected[2] = { record->baseline_output, record->perturbed_output };
  const double *computed[2] = { computed_baseline, computed_perturbed };
  const char *gate_names[3] = { "baseline", "perturbed", "response" };
  for(size_t component = 0; component < record->output_count; ++component) {
    if(!isfinite(computed_normalization[component])
       || !(computed_normalization[component] > 0.0)
       || computed_normalization[component] != record->normalization[component]) {
      if(error != NULL && error_size > 0) {
        snprintf(
              error, error_size,
              "transport %s component %zu has invalid input normalization",
              record->case_id, component);
      }
      return 0;
    }
    for(int gate = 0; gate < 3; ++gate) {
      const double actual = gate < 2 ? computed[gate][component]
                                     : computed[1][component] - computed[0][component];
      const double reference = gate < 2
                                     ? expected[gate][component]
                                     : expected[1][component] - expected[0][component];
      double denominator = fmax(1.0e-300, fmax(fabs(actual), fabs(reference)));
      double bound = 2.0e-12;
      if(gate == 2) {
        /* A difference of nearby outputs must not be tested relative to its
         * nearly zero value. Propagate the two retained endpoint bounds.
         * The three double subtractions contribute at most 4*DBL_EPSILON
         * times the summed endpoint scales (using epsilon conservatively
         * rather than half-epsilon). This is not independent sensitivity
         * evidence beyond the endpoint comparisons. */
        const double baseline_scale = fmax(
              1.0e-300,
              fmax(fabs(computed[0][component]), fabs(expected[0][component])));
        const double perturbed_scale = fmax(
              1.0e-300,
              fmax(fabs(computed[1][component]), fabs(expected[1][component])));
        denominator = baseline_scale + perturbed_scale;
        bound += 4.0 * DBL_EPSILON;
      }
      double scaled_error = 0.0;
      if(!m1_thcm1_fixture_compare_one(
               actual, reference, denominator, bound, 0.0, 0.0, &scaled_error)) {
        if(report != NULL) {
          report->max_scaled_error = fmax(report->max_scaled_error, scaled_error);
          report->component = component;
          report->gate = gate_names[gate];
        }
        if(error != NULL && error_size > 0) {
          snprintf(
                error, error_size, "transport %s component %zu %s comparison failed",
                record->case_id, component, gate_names[gate]);
        }
        return 0;
      }
      if(report != NULL) {
        report->max_scaled_error = fmax(report->max_scaled_error, scaled_error);
      }
    }
  }
  return 1;
}

/* Standard offline replay: evaluate the current library only on the baseline
 * input and use the retained THC perturbation as the response envelope. */
static inline int m1_thcm1_transport_compare_baseline_response(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict policy,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size) {
  const int supported_policy
        = policy != NULL
          && (strcmp(policy, "strict_relative_2e-12_propagated_response_v1") == 0
              || strcmp(policy, "strict_relative_2e-12_propagated_response_current_v2")
                       == 0);
  if(record == NULL || !supported_policy
     || !m1_thcm1_fixture_pair_id_matches_case_suffix(record, "__pair")) {
    m1_thcm1_fixture_set_error(
          error, error_size, "invalid transport baseline-response metadata");
    return 0;
  }
  return m1_thcm1_fixture_compare_baseline_response(
        record, policy, computed_normalization, computed_baseline, report, error,
        error_size);
}

#endif /* UNIT_TESTS_M1_THCM1_TRANSPORT_FIXTURE_H_ */
