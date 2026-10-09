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

/* Ordinary records reach the current source API at both paired endpoints and
 * use the common baseline/perturbed/response comparator.  The six named
 * zero-flux policy cases have local admissibility/exchange checks, not
 * reference agreement.  records_run counts only paired-agreement records. */
int m1_thcm1_run_instantaneous_source_fixtures(
      const char *restrict fixture_path,
      size_t *restrict records_run,
      char *restrict error,
      const size_t error_size);

#endif /* UNIT_TESTS_M1_THCM1_SOURCE_FIXTURE_H_ */
