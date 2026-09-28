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
int m1_thcm1_transport_prepare_volume_weighted(
      const double cell_volumes[4],
      const double face_volume,
      const double state_stencil[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      const double physical_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      const double physical_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      double weighted_state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      double weighted_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      double weighted_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      char *restrict error,
      const size_t error_size);

/* Match the existing fixture producer's scale convention after the same
 * volume preparation: face physical operands and speed times the first two
 * stencil rows.  The stored normalization is metadata, not a tolerance. */
int m1_thcm1_transport_compute_weighted_normalization(
      const double input[M1_THCM1_TRANSPORT_INPUT_COUNT],
      double normalization[M1_THCM1_TRANSPORT_OUTPUT_COUNT]);

/* Stored transport replay and focused comparator tests use the shared
 * strict_relative_2e-12_propagated_response_v1 two-state comparator. */
static inline int m1_thcm1_transport_compare_paired(
      const m1_thcm1_fixture_record *restrict record,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      const double *restrict computed_perturbed,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size) {
  return m1_thcm1_fixture_compare_paired_strict_relative(
        record, computed_normalization, computed_baseline, computed_perturbed, report,
        error, error_size);
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
