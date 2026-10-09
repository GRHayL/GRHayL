#ifndef UNIT_TESTS_M1_THCM1_STRESS_ENERGY_FIXTURE_H_
#define UNIT_TESTS_M1_THCM1_STRESS_ENERGY_FIXTURE_H_

/*
 * Replay the portable stress_energy pairs exported from the Verification/
 * CL-04 campaign.  The fixture contains the complete consumed radiation
 * state, ADM metric, and fluid velocity for each pair.  The current public
 * implementation is evaluated at both endpoints and compared with the stored
 * THC_M1 baseline, perturbed endpoint, and paired response.
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
#define M1_THCM1_STRESS_ENERGY_FIXTURE_INPUT_COUNT  21
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

#define M1_THCM1_STRESS_ENERGY_EPSILON_C             1.0e-10
#define M1_THCM1_STRESS_ENERGY_ENERGY_FLOOR          1.0e-12
#define M1_THCM1_STRESS_ENERGY_ZETA_MIN              1.0e-8
#define M1_THCM1_STRESS_ENERGY_FD_EPSILON_REL        1.0e-6
#define M1_THCM1_STRESS_ENERGY_FD_EPSILON_ABS        1.0e-12
#define M1_THCM1_STRESS_ENERGY_NEWTON_MAX_ITERATIONS 20
#define M1_THCM1_STRESS_ENERGY_NEWTON_TOLERANCE      1.0e-10

int m1_thcm1_run_stress_energy_fixtures(
      const char *restrict fixture_path,
      size_t *restrict records_run,
      char *restrict error,
      const size_t error_size);

#endif /* UNIT_TESTS_M1_THCM1_STRESS_ENERGY_FIXTURE_H_ */
