#ifndef UNIT_TESTS_M1_THCM1_FIXTURE_UTILS_H_
#define UNIT_TESTS_M1_THCM1_FIXTURE_UTILS_H_

/*
 * Test-local reader and paired comparator for frozen THC_M1 observations.
 *
 * External TestData fixtures use a versioned little-endian binary stream.
 * The original token stream remains readable for local fixtures and parser
 * rejection checks. Neither format depends on native struct layout or an
 * installed GRHayL test API. Consumers validate operation-specific vector
 * counts before comparing the recorded values.
 */

#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define M1_THCM1_FIXTURE_MAGIC    "M1_THCM1_FIXTURE"
#define M1_THCM1_FIXTURE_VERSION  1
/* Caller-owned diagnostic storage; this is not a parser/input-text cap. */
#define M1_THCM1_FIXTURE_TEXT_MAX 256
typedef struct {
  char *case_id;
  char *pair_id;
  char *origin;
  char *seed_id;
  char *family;
  char *perturbation;
  size_t sensitivity_start;
  size_t sensitivity_count;
  size_t input_count;
  size_t output_count;
  size_t normalization_count;
  double *baseline_input;
  double *perturbed_input;
  double *baseline_output;
  double *perturbed_output;
  double *normalization;
  int baseline_available;
  int perturbed_available;
  int baseline_status;
  int perturbed_status;
  int baseline_reference_status;
  int perturbed_reference_status;
} m1_thcm1_fixture_record;

typedef struct {
  char *operation;
  char *policy;
  size_t record_count;
  m1_thcm1_fixture_record *records;
} m1_thcm1_fixture_collection;

typedef struct {
  double max_scaled_error;
  size_t component;
  const char *gate;
} m1_thcm1_fixture_comparison_report;

static inline void m1_thcm1_fixture_set_error(
      char *restrict error,
      const size_t error_size,
      const char *restrict message) {
  if(error != NULL && error_size > 0) {
    snprintf(error, error_size, "%s", message != NULL ? message : "fixture error");
  }
}

void m1_thcm1_fixture_free(m1_thcm1_fixture_collection *restrict collection);

static inline int
m1_thcm1_fixture_finite_vector(const double *restrict values, const size_t count) {
  if(values == NULL || count == 0) {
    return 0;
  }
  for(size_t i = 0; i < count; ++i) {
    if(!isfinite(values[i])) {
      return 0;
    }
  }
  return 1;
}

/* A pair identifier is a separate provenance key, not a second spelling of
 * the case identifier.  Operation owners apply their own producer-specific
 * relationship below; this common check prevents absent or collapsed pair
 * metadata from entering any comparator. */
static inline int
m1_thcm1_fixture_pair_metadata_valid(const m1_thcm1_fixture_record *restrict record) {
  return record != NULL && record->case_id != NULL && record->pair_id != NULL
         && record->case_id[0] != '\0' && record->pair_id[0] != '\0'
         && strcmp(record->case_id, record->pair_id) != 0;
}

static inline int m1_thcm1_fixture_pair_id_matches_case_suffix(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict suffix) {
  if(!m1_thcm1_fixture_pair_metadata_valid(record) || suffix == NULL) {
    return 0;
  }
  const size_t case_length = strlen(record->case_id);
  const size_t suffix_length = strlen(suffix);
  const size_t pair_length = strlen(record->pair_id);
  return pair_length == case_length + suffix_length
         && strncmp(record->pair_id, record->case_id, case_length) == 0
         && strcmp(record->pair_id + case_length, suffix) == 0;
}

int m1_thcm1_fixture_record_valid(
      const m1_thcm1_fixture_record *restrict record,
      const size_t expected_input_count,
      const size_t expected_output_count);

/* Load one operation.  A zero expected count means "operation-defined" and
 * is useful to the harness rejection tests; normal consumers pass their exact
 * schema counts. */
int m1_thcm1_fixture_load(
      const char *restrict path,
      const char *restrict operation,
      const size_t expected_input_count,
      const size_t expected_output_count,
      m1_thcm1_fixture_collection *restrict collection,
      char *restrict error,
      const size_t error_size);

/* Compare one current baseline evaluation with the retained THC baseline and
 * its paired perturbation response.  The operation policy is deliberately
 * enumerated here: a caller cannot silently obtain a generic fallback by
 * changing a fixture header. */
int m1_thcm1_fixture_compare_baseline_response(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict policy,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size);

/* The constants are copied from the recorded current_official A1/A2 policy.
 * They are a fixture policy, not a new production tolerance. */
int m1_thcm1_fixture_compare_paired(
      const m1_thcm1_fixture_record *restrict record,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      const double *restrict computed_perturbed,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size);

/* Explicit two-state comparator for the
 * strict_relative_2e-12_propagated_response_v1 policy. Stored transport and
 * stress-energy replay use this symmetric campaign rule: each endpoint uses
 * max(1e-300, |actual|, |reference|), and the paired response propagates the
 * two endpoint scales. The common trusted-baseline helper is asymmetric and
 * uses the perturbed endpoint as an envelope, so it cannot express this
 * contract without changing acceptance decisions. The fixture's
 * input-derived normalization is still checked separately by the consumer;
 * it is not used as a tolerance knob here. */
int m1_thcm1_fixture_compare_paired_strict_relative(
      const m1_thcm1_fixture_record *restrict record,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      const double *restrict computed_perturbed,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size);

/* Compare one current implementation result with a trusted/perturbed THC
 * pair.  This is the minimal offline replay shape used by the pointwise
 * campaign: the current implementation is evaluated only at the baseline
 * input, and the retained response participates in the target-local bound. */
int m1_thcm1_fixture_compare_envelope(
      const m1_thcm1_fixture_record *restrict record,
      const char *restrict policy,
      const double *restrict computed_normalization,
      const double *restrict computed_baseline,
      m1_thcm1_fixture_comparison_report *restrict report,
      char *restrict error,
      const size_t error_size);

#endif /* UNIT_TESTS_M1_THCM1_FIXTURE_UTILS_H_ */
