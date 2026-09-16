#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "m1_test_utils.h"
#include "m1_thcm1_transport_fixture.h"

/*
 * Local and retained-reference coverage for the canonical four-point blended
 * Rusanov transport operation. The test checks limiter/opacity policies,
 * metric and volume preparation, validation, and both transport fixture modes.
 */

static void fail_test(const char *message) {
  ghl_error("unit_test_m1_thcm1_blended_rusanov: %s\n", message);
}

static ghl_error_codes_t m1_thcm1_call_volume_weighted_transport(
      const ghl_m1_parameters *restrict m1_params,
      const double state_stencil[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      const double physical_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      const double physical_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      const double speed_L,
      const double speed_R,
      const double kappa_face,
      const double delta_x,
      const bool diffusion_correction_enabled,
      double flux_tilde[M1_THCM1_TRANSPORT_OUTPUT_COUNT],
      ghl_m1_four_point_transport_diagnostics *restrict diagnostics) {
  return ghl_m1_compute_neutrino_four_point_volume_weighted_transport_flux(
      m1_params, state_stencil, physical_flux_L, physical_flux_R, speed_L,
      speed_R, kappa_face, delta_x, diffusion_correction_enabled, flux_tilde,
      diagnostics);
}

static void check_close(const double actual, const double expected,
                        const char *label) {
  if(!m1_nearly_equal(actual, expected, 2.0e-13, 2.0e-14))
    fail_test(label);
}

typedef struct {
  uint64_t state;
} four_point_rng;

static uint64_t four_point_rng_next(four_point_rng *restrict rng) {
  uint64_t z = (rng->state += UINT64_C(0x9e3779b97f4a7c15));
  z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
  return z ^ (z >> 31);
}

static double four_point_rng_unit(four_point_rng *restrict rng) {
  return (double)(four_point_rng_next(rng) >> 11) * 0x1.0p-53;
}

static double four_point_rng_between(
      four_point_rng *restrict rng, const double lower, const double upper) {
  return lower + (upper - lower) * four_point_rng_unit(rng);
}

static bool same_nonzero_sign(const double left, const double right) {
  return (left > 0.0 && right > 0.0) ||
         (left < 0.0 && right < 0.0);
}

static bool opposite_nonzero_sign(const double left, const double right) {
  return (left > 0.0 && right < 0.0) ||
         (left < 0.0 && right > 0.0);
}

static double expected_limiter(
      const double dum, const double duc, const double dup,
      const double theta, bool *restrict sawtooth) {
  *sawtooth = false;
  if(same_nonzero_sign(dup, duc) && same_nonzero_sign(dum, duc)) {
    return fmin(1.0, fmin(theta * dum / duc, theta * dup / duc));
  }
  if(opposite_nonzero_sign(dup, duc) &&
     opposite_nonzero_sign(dum, duc))
    *sawtooth = true;
  return 0.0;
}

static double expected_opacity(
      const double kappa_face, const double delta,
      const double mindiss) {
  const double optical_width = kappa_face * delta;
  if(optical_width <= 1.0)
    return 1.0;
  return fmax(mindiss, fmin(1.0, 1.0 / optical_width));
}

static void check_full_four_point_operator(
      ghl_m1_parameters *restrict params) {
  four_point_rng rng = {.state = UINT64_C(0x4d315f464f4f5552)};

  for(int case_index = 0; case_index < 128; ++case_index) {
    ghl_metric_quantities metric;
    if(case_index & 1)
      m1_setup_metric_anchor_B(&metric);
    else
      m1_setup_flat_metric(&metric);

    params->minmod_theta = four_point_rng_between(&rng, 0.35, 1.85);
    params->mindiss = (case_index % 4 == 0) ? 0.0 :
                     four_point_rng_between(&rng, 0.05, 0.35);
    const double kappa_face = case_index % 3 == 0 ? 0.2 :
                              (case_index % 3 == 1 ? 4.0 : 100.0);
    const double delta = case_index % 3 == 0 ? 0.5 : 1.0;
    const double speed_L = four_point_rng_between(&rng, 0.05, 0.9);
    const double speed_R = four_point_rng_between(&rng, 0.05, 0.9);
    const double face_speed = fmax(speed_L, speed_R);
    double state_stencil[4][ghl_m1_neutrino_transport_component_count];
    double physical_flux_L[ghl_m1_neutrino_transport_component_count];
    double physical_flux_R[ghl_m1_neutrino_transport_component_count];
    double increments[3];

    for(int component = 0;
        component < ghl_m1_neutrino_transport_component_count; ++component) {
      const double base = four_point_rng_between(&rng, -1.0, 1.0);
      const double scale = four_point_rng_between(&rng, 0.1, 0.9);
      switch((case_index + component) % 4) {
        case 0:
          increments[0] = scale;
          increments[1] = 0.7 * scale;
          increments[2] = 1.3 * scale;
          break;
        case 1:
          increments[0] = -scale;
          increments[1] = -0.7 * scale;
          increments[2] = -1.3 * scale;
          break;
        case 2:
          increments[0] = -scale;
          increments[1] = 0.7 * scale;
          increments[2] = -1.3 * scale;
          break;
        default:
          increments[0] = scale;
          increments[1] = 0.0;
          increments[2] = -1.3 * scale;
          break;
      }
      state_stencil[0][component] = base;
      state_stencil[1][component] = base + increments[0];
      state_stencil[2][component] = state_stencil[1][component] + increments[1];
      state_stencil[3][component] = state_stencil[2][component] + increments[2];
      physical_flux_L[component] = four_point_rng_between(&rng, -2.0, 2.0);
      physical_flux_R[component] = four_point_rng_between(&rng, -2.0, 2.0);
    }

    double flux_tilde[ghl_m1_neutrino_transport_component_count];
    for(int component = 0;
        component < ghl_m1_neutrino_transport_component_count; ++component)
      flux_tilde[component] = 91.0 + component;
    ghl_m1_four_point_transport_diagnostics diagnostics = {0};
    if(ghl_m1_compute_neutrino_four_point_transport_flux(
           params, &metric, state_stencil, physical_flux_L, physical_flux_R,
           speed_L, speed_R, kappa_face, delta, false, flux_tilde,
           &diagnostics) != ghl_success)
      fail_test("randomized four-point transport operation failed");

    const double opacity = expected_opacity(kappa_face, delta, params->mindiss);
    check_close(diagnostics.opacity_suppression, opacity,
                "four-point opacity diagnostic mismatch");
    check_close(diagnostics.face_speed, face_speed,
                "four-point face speed diagnostic mismatch");
    for(int component = 0;
        component < ghl_m1_neutrino_transport_component_count; ++component) {
      bool expected_sawtooth = false;
      const double phi = expected_limiter(
          state_stencil[1][component] - state_stencil[0][component],
          state_stencil[2][component] - state_stencil[1][component],
          state_stencil[3][component] - state_stencil[2][component],
          params->minmod_theta, &expected_sawtooth);
      const double flux_high = 0.5 * (physical_flux_L[component] +
                                      physical_flux_R[component]);
      const double flux_low = flux_high - 0.5 * face_speed *
          (state_stencil[2][component] - state_stencil[1][component]);
      const double dissipation = expected_sawtooth ? 1.0 : opacity;
      const double expected = metric.sqrt_detgamma *
          (flux_high - dissipation * (1.0 - phi) *
           (flux_high - flux_low));
      check_close(flux_tilde[component], expected,
                  "four-point transport flux mismatch");
      check_close(diagnostics.phi[component], phi,
                  "four-point limiter diagnostic mismatch");
      if(diagnostics.sawtooth[component] != expected_sawtooth)
        fail_test("four-point sawtooth diagnostic mismatch");
    }

    /* The optional diagnostics packet must not be required by the operator. */
    double no_diagnostic_flux[ghl_m1_neutrino_transport_component_count];
    if(ghl_m1_compute_neutrino_four_point_transport_flux(
           params, &metric, state_stencil, physical_flux_L, physical_flux_R,
           speed_L, speed_R, kappa_face, delta, false, no_diagnostic_flux,
           NULL) != ghl_success)
      fail_test("four-point transport rejected NULL diagnostics");
    for(int component = 0;
        component < ghl_m1_neutrino_transport_component_count; ++component)
      check_close(no_diagnostic_flux[component], flux_tilde[component],
                  "four-point NULL-diagnostics result mismatch");
  }

  /* Invalid policy and prepared operands are rejected transactionally. */
  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  double state_stencil[4][ghl_m1_neutrino_transport_component_count] = {{0.0}};
  double physical_flux_L[ghl_m1_neutrino_transport_component_count] = {0.0};
  double physical_flux_R[ghl_m1_neutrino_transport_component_count] = {0.0};
  double flux_tilde[ghl_m1_neutrino_transport_component_count];
  for(int component = 0;
      component < ghl_m1_neutrino_transport_component_count; ++component)
    flux_tilde[component] = 101.0 + component;
  const double flux_before[ghl_m1_neutrino_transport_component_count] = {
      101.0, 102.0, 103.0, 104.0, 105.0};
  ghl_m1_four_point_transport_diagnostics diagnostics = {
      .opacity_suppression = 201.0, .face_speed = 202.0};
  const ghl_m1_four_point_transport_diagnostics diagnostics_before = diagnostics;
  if(ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &metric, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.7, 1.0, 1.0, true, flux_tilde, &diagnostics) !=
         ghl_error_m1_incompatible_transport_policy)
    fail_test("incompatible four-point transport policy was accepted");
  for(int component = 0;
      component < ghl_m1_neutrino_transport_component_count; ++component)
    if(flux_tilde[component] != flux_before[component])
      fail_test("policy rejection changed four-point flux");
  if(diagnostics.opacity_suppression != diagnostics_before.opacity_suppression ||
     diagnostics.face_speed != diagnostics_before.face_speed)
    fail_test("policy rejection changed four-point diagnostics");

  state_stencil[0][0] = NAN;
  for(int component = 0;
      component < ghl_m1_neutrino_transport_component_count; ++component)
    flux_tilde[component] = 111.0 + component;
  if(ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &metric, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.7, 1.0, 1.0, false, flux_tilde, NULL) !=
         ghl_error_m1_invalid_state)
    fail_test("invalid four-point stencil was accepted");
  for(int component = 0;
      component < ghl_m1_neutrino_transport_component_count; ++component)
    if(flux_tilde[component] != 111.0 + component)
      fail_test("invalid four-point stencil changed flux");

  metric.sqrt_detgamma = 0.0;
  state_stencil[0][0] = 0.0;
  for(int component = 0;
      component < ghl_m1_neutrino_transport_component_count; ++component)
    flux_tilde[component] = 121.0 + component;
  if(ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &metric, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.7, 1.0, 1.0, false, flux_tilde, NULL) !=
         ghl_error_m1_invalid_metric)
    fail_test("invalid four-point metric was accepted");
  for(int component = 0;
      component < ghl_m1_neutrino_transport_component_count; ++component)
    if(flux_tilde[component] != 121.0 + component)
      fail_test("invalid four-point metric changed flux");
}

static void check_volume_preparation_adapter(void) {
  double state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double physical_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double physical_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  for(int cell = 0; cell < 4; ++cell)
    for(int component = 0;
        component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
      state[cell][component] = 0.25 + 0.3 * cell + 0.02 * component;
  for(int component = 0;
      component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
    physical_flux_L[component] = 0.7 + 0.03 * component;
    physical_flux_R[component] = -0.2 + 0.04 * component;
  }
  const double cell_volumes[4] = {0.75, 1.10, 1.40, 1.85};
  const double face_volume = 0.90;
  double weighted_state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double weighted_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double weighted_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  char error[128] = {0};
  if(!m1_thcm1_transport_prepare_volume_weighted(
         cell_volumes, face_volume, state, physical_flux_L, physical_flux_R,
         weighted_state, weighted_flux_L, weighted_flux_R, error,
         sizeof(error)))
    fail_test("valid volume preparation was rejected");
  for(int cell = 0; cell < 4; ++cell)
    for(int component = 0;
        component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
      check_close(weighted_state[cell][component],
                  cell_volumes[cell] * state[cell][component],
                  "weighted stencil state mismatch");
  for(int component = 0;
      component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
    check_close(weighted_flux_L[component],
                face_volume * physical_flux_L[component],
                "weighted left face flux mismatch");
    check_close(weighted_flux_R[component],
                face_volume * physical_flux_R[component],
                "weighted right face flux mismatch");
  }

  const double state_before[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {
      {91.0, 92.0, 93.0, 94.0, 95.0},
      {96.0, 97.0, 98.0, 99.0, 100.0},
      {101.0, 102.0, 103.0, 104.0, 105.0},
      {106.0, 107.0, 108.0, 109.0, 110.0}};
  const double flux_L_before[M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {
      111.0, 112.0, 113.0, 114.0, 115.0};
  const double flux_R_before[M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {
      116.0, 117.0, 118.0, 119.0, 120.0};
  memcpy(weighted_state, state_before, sizeof(weighted_state));
  memcpy(weighted_flux_L, flux_L_before, sizeof(weighted_flux_L));
  memcpy(weighted_flux_R, flux_R_before, sizeof(weighted_flux_R));
  double invalid_volumes[4] = {0.75, 0.0, 1.40, 1.85};
  if(m1_thcm1_transport_prepare_volume_weighted(
         invalid_volumes, face_volume, state, physical_flux_L, physical_flux_R,
         weighted_state, weighted_flux_L, weighted_flux_R, error,
         sizeof(error)))
    fail_test("nonpositive cell volume was accepted");
  if(memcmp(weighted_state, state_before, sizeof(weighted_state)) != 0 ||
     memcmp(weighted_flux_L, flux_L_before, sizeof(weighted_flux_L)) != 0 ||
     memcmp(weighted_flux_R, flux_R_before, sizeof(weighted_flux_R)) != 0)
    fail_test("invalid volume preparation changed outputs");
  if(m1_thcm1_transport_prepare_volume_weighted(
         cell_volumes, NAN, state, physical_flux_L, physical_flux_R,
         weighted_state, weighted_flux_L, weighted_flux_R, error,
         sizeof(error)))
    fail_test("nonfinite face volume was accepted");

  double overflow_state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  memcpy(overflow_state, state, sizeof(overflow_state));
  overflow_state[0][0] = 2.0;
  const double overflow_volumes[4] = {
      DBL_MAX, 1.0, 1.0, 1.0};
  if(m1_thcm1_transport_prepare_volume_weighted(
         overflow_volumes, 1.0, overflow_state, physical_flux_L,
         physical_flux_R, weighted_state, weighted_flux_L, weighted_flux_R,
         error, sizeof(error)))
    fail_test("overflowing weighted state was accepted");
  if(memcmp(weighted_state, state_before, sizeof(weighted_state)) != 0 ||
     memcmp(weighted_flux_L, flux_L_before, sizeof(weighted_flux_L)) != 0 ||
     memcmp(weighted_flux_R, flux_R_before, sizeof(weighted_flux_R)) != 0)
    fail_test("overflowing volume preparation changed outputs");
}

static void check_prepared_transport_local_contract(void) {
  check_volume_preparation_adapter();

  ghl_m1_parameters params = {0};
  if(ghl_m1_initialize(1.0e-8, 1.0e-12, 1.0e-8, 1.0e-6,
                       1.0e-10, 100, 1.0e-10, &params) != ghl_success) {
    fail_test("prepared transport M1 initialization failed");
    return;
  }
  params.minmod_theta = 1.4;
  params.mindiss = 0.2;
  double state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double physical_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double physical_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  for(int cell = 0; cell < 4; ++cell)
    for(int component = 0;
        component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
      const double cell_value[4] = {0.10, 0.30, 0.70, 0.65};
      state[cell][component] = cell_value[cell] + 0.03 * component;
    }
  for(int component = 0;
      component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
    physical_flux_L[component] = 0.75 + 0.04 * component;
    physical_flux_R[component] = -0.25 + 0.06 * component;
  }

  const double volume = 1.75;
  const double constant_volumes[4] = {volume, volume, volume, volume};
  double weighted_state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double weighted_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double weighted_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  char error[128] = {0};
  if(!m1_thcm1_transport_prepare_volume_weighted(
         constant_volumes, volume, state, physical_flux_L, physical_flux_R,
         weighted_state, weighted_flux_L, weighted_flux_R, error,
         sizeof(error))) {
    fail_test("constant volume preparation failed");
    return;
  }
  ghl_metric_quantities metric;
  const double spatial_metric = pow(volume, 2.0 / 3.0);
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0,
                        spatial_metric, 0.0, 0.0,
                        spatial_metric, 0.0, spatial_metric, &metric);
  double pointwise_flux[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  double prepared_flux[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  ghl_m1_four_point_transport_diagnostics pointwise_diagnostics = {0};
  ghl_m1_four_point_transport_diagnostics prepared_diagnostics = {0};
  const ghl_error_codes_t pointwise_status =
      ghl_m1_compute_neutrino_four_point_transport_flux(
         &params, &metric, state, physical_flux_L, physical_flux_R,
         0.35, 0.65, 3.0, 0.7, false, pointwise_flux,
         &pointwise_diagnostics);
  const ghl_error_codes_t prepared_status = m1_thcm1_call_volume_weighted_transport(
         &params, weighted_state, weighted_flux_L, weighted_flux_R,
         0.35, 0.65, 3.0, 0.7, false, prepared_flux,
         &prepared_diagnostics);
  if(pointwise_status != ghl_success || prepared_status != ghl_success) {
    ghl_error("unit_test_m1_thcm1_blended_rusanov: constant-volume statuses "
              "pointwise=%d prepared=%d\n", pointwise_status, prepared_status);
    fail_test("constant volume reduction operation failed");
    return;
  }
  for(int component = 0;
      component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
    check_close(prepared_flux[component], pointwise_flux[component],
                "prepared constant-volume reduction mismatch");

  /* A weighted call through the old API applies the face factor twice.  The
   * deliberately wrong result must differ from the prepared result. */
  double double_densitized[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  if(ghl_m1_compute_neutrino_four_point_transport_flux(
         &params, &metric, weighted_state, weighted_flux_L, weighted_flux_R,
         0.35, 0.65, 3.0, 0.7, false, double_densitized, NULL) !=
         ghl_success)
    fail_test("double-densitization mutation could not be evaluated");
  bool detected_double_densitization = false;
  for(int component = 0;
      component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
    detected_double_densitization |= !m1_nearly_equal(
        double_densitized[component], prepared_flux[component],
        1.0e-12, 1.0e-14);
  if(!detected_double_densitization)
    fail_test("double-densitization mutation was not detected");

  /* A failed prepared call must leave both outputs and every diagnostic field
   * untouched.  Initialize every member so this checks the complete public
   * diagnostic packet, not only opacity and face speed. */
  double invalid_state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  memcpy(invalid_state, weighted_state, sizeof(invalid_state));
  invalid_state[0][0] = NAN;
  double invalid_flux[M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {
      131.0, 132.0, 133.0, 134.0, 135.0};
  const double invalid_flux_before[M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {
      131.0, 132.0, 133.0, 134.0, 135.0};
  ghl_m1_four_point_transport_diagnostics invalid_diagnostics = {
      .phi = {141.0, 142.0, 143.0, 144.0, 145.0},
      .sawtooth = {true, false, true, false, true},
      .opacity_suppression = 146.0,
      .face_speed = 147.0};
  const ghl_m1_four_point_transport_diagnostics diagnostics_before =
      invalid_diagnostics;
  if(m1_thcm1_call_volume_weighted_transport(
         &params, invalid_state, weighted_flux_L, weighted_flux_R,
         0.35, 0.65, 3.0, 0.7, false, invalid_flux,
         &invalid_diagnostics) != ghl_error_m1_invalid_state)
    fail_test("invalid prepared state was accepted");
  if(memcmp(invalid_flux, invalid_flux_before, sizeof(invalid_flux)) != 0)
    fail_test("invalid prepared state changed the output flux");
  bool diagnostics_unchanged =
      invalid_diagnostics.opacity_suppression ==
          diagnostics_before.opacity_suppression &&
      invalid_diagnostics.face_speed == diagnostics_before.face_speed;
  for(int component = 0;
      component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
    diagnostics_unchanged &=
        invalid_diagnostics.phi[component] ==
            diagnostics_before.phi[component] &&
        invalid_diagnostics.sawtooth[component] ==
            diagnostics_before.sawtooth[component];
  if(!diagnostics_unchanged)
    fail_test("invalid prepared state changed diagnostics");

  const double varying_volumes[4] = {0.75, 1.10, 1.40, 1.85};
  if(!m1_thcm1_transport_prepare_volume_weighted(
         varying_volumes, 0.90, state, physical_flux_L, physical_flux_R,
         weighted_state, weighted_flux_L, weighted_flux_R, error,
         sizeof(error))) {
    fail_test("varying volume preparation failed");
    return;
  }
  double varying_flux[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  if(m1_thcm1_call_volume_weighted_transport(
         &params, weighted_state, weighted_flux_L, weighted_flux_R,
         0.35, 0.65, 3.0, 0.7, false, varying_flux, NULL) != ghl_success) {
    fail_test("varying volume prepared operation failed");
    return;
  }
  for(int changed_cell = 0; changed_cell < 4; ++changed_cell) {
    double mutated_volumes[4];
    memcpy(mutated_volumes, varying_volumes, sizeof(mutated_volumes));
    mutated_volumes[changed_cell] *= 1.13;
    if(!m1_thcm1_transport_prepare_volume_weighted(
           mutated_volumes, 0.90, state, physical_flux_L, physical_flux_R,
           weighted_state, weighted_flux_L, weighted_flux_R, error,
           sizeof(error))) {
      fail_test("volume mutation preparation failed");
      return;
    }
    double mutated_flux[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
    if(m1_thcm1_call_volume_weighted_transport(
           &params, weighted_state, weighted_flux_L, weighted_flux_R,
           0.35, 0.65, 3.0, 0.7, false, mutated_flux, NULL) != ghl_success) {
      fail_test("volume mutation prepared operation failed");
      return;
    }
    bool changed = false;
    for(int component = 0;
        component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
      changed |= !m1_nearly_equal(mutated_flux[component], varying_flux[component],
                                  1.0e-12, 1.0e-14);
    if(!changed) {
      ghl_error("unit_test_m1_thcm1_blended_rusanov: unchanged cell volume "
                "%d in local mutation\n", changed_cell);
      fail_test("prepared transport ignored a stencil-cell volume");
    }
  }
  if(!m1_thcm1_transport_prepare_volume_weighted(
         varying_volumes, 1.07, state, physical_flux_L, physical_flux_R,
         weighted_state, weighted_flux_L, weighted_flux_R, error,
         sizeof(error)))
    fail_test("face-volume mutation preparation failed");
  double mutated_face_flux[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
  if(m1_thcm1_call_volume_weighted_transport(
         &params, weighted_state, weighted_flux_L, weighted_flux_R,
         0.35, 0.65, 3.0, 0.7, false, mutated_face_flux, NULL) != ghl_success)
    fail_test("face-volume mutation prepared operation failed");
  bool changed_face = false;
  for(int component = 0;
      component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
    changed_face |= !m1_nearly_equal(mutated_face_flux[component], varying_flux[component],
                                     1.0e-12, 1.0e-14);
  if(!changed_face)
    fail_test("prepared transport ignored the face volume");
}

static const char *m1_thcm1_transport_direction(const char *case_id) {
  if(case_id != NULL && strstr(case_id, "__d0__") != NULL) return "d0";
  if(case_id != NULL && strstr(case_id, "__d1__") != NULL) return "d1";
  if(case_id != NULL && strstr(case_id, "__d2__") != NULL) return "d2";
  return "unknown";
}

static char *m1_thcm1_transport_make_path(
      const char *restrict fixture_dir,
      const char *restrict fixture_name) {
  const size_t directory_length = strlen(fixture_dir);
  const size_t fixture_name_length = strlen(fixture_name);
  char *path = (char *)malloc(directory_length + fixture_name_length + 1);
  if(path == NULL) return NULL;
  memcpy(path, fixture_dir, directory_length);
  memcpy(path + directory_length, fixture_name, fixture_name_length + 1);
  return path;
}

static int m1_thcm1_transport_file_exists(const char *path) {
  FILE *file = fopen(path, "r");
  if(file == NULL) return 0;
  fclose(file);
  return 1;
}

static int m1_thcm1_load_variable_transport_fixture(
      const char *restrict path,
      m1_thcm1_fixture_collection *restrict collection,
      char *restrict error,
      const size_t error_size) {
  /* Keep the variable corpus distinct from the pointwise corpus.  A file
   * carrying the old operation label must not silently enter the prepared
   * volume-weighted comparison. */
  return m1_thcm1_fixture_load(
      path, "neutrino_four_point_prepared_transport_flux",
      M1_THCM1_TRANSPORT_INPUT_COUNT, M1_THCM1_TRANSPORT_OUTPUT_COUNT,
      collection, error, error_size);
}

static void check_variable_transport_fixtures(const char *restrict fixture_dir) {
  static const char *const fixture_names[] = {
      "/transport_four_point_varying_d0.dat",
      "/transport_four_point_varying_d1.dat",
      "/transport_four_point_varying_d2.dat",
      "/transport_four_point_varying_controls.dat"};
  const size_t shard_count = sizeof(fixture_names) / sizeof(fixture_names[0]);
  m1_thcm1_fixture_collection collections[4] = {{0}};
  size_t present_count = 0;
  for(size_t shard = 0; shard < shard_count; ++shard) {
    char *path = m1_thcm1_transport_make_path(fixture_dir, fixture_names[shard]);
    if(path == NULL) {
      fail_test("variable-volume fixture path allocation failed");
      return;
    }
    if(m1_thcm1_transport_file_exists(path))
      ++present_count;
    else
      ghl_error("unit_test_m1_thcm1_blended_rusanov: missing required "
                "variable-volume fixture %s\n", path);
    free(path);
  }
  if(present_count != shard_count) {
    fail_test("variable-volume fixture family is incomplete; no trusted "
              "outputs were fabricated");
    return;
  }

  size_t total_records = 0;
  for(size_t shard = 0; shard < shard_count; ++shard) {
    char *path = m1_thcm1_transport_make_path(fixture_dir, fixture_names[shard]);
    if(path == NULL) {
      fail_test("variable-volume fixture path allocation failed");
      for(size_t cleanup = 0; cleanup < shard_count; ++cleanup)
        m1_thcm1_fixture_free(&collections[cleanup]);
      return;
    }
    char error[256] = {0};
    if(!m1_thcm1_load_variable_transport_fixture(
           path, &collections[shard], error, sizeof(error))) {
      fail_test(error[0] != '\0' ? error :
                "variable-volume fixture load failed");
      free(path);
      for(size_t cleanup = 0; cleanup < shard_count; ++cleanup)
        m1_thcm1_fixture_free(&collections[cleanup]);
      return;
    }
    free(path);
    if(strcmp(collections[shard].policy,
              "strict_relative_2e-12_propagated_response_v1") != 0) {
      fail_test("variable-volume fixture uses an unsupported comparison policy");
      for(size_t cleanup = 0; cleanup < shard_count; ++cleanup)
        m1_thcm1_fixture_free(&collections[cleanup]);
      return;
    }
    total_records += collections[shard].record_count;
  }

  size_t failures = 0;
  size_t changed_records = 0;
  size_t control_records = 0;
  bool saw_nonconstant_volumes = false;
  for(size_t shard = 0; shard < shard_count; ++shard)
    for(size_t index = 0; index < collections[shard].record_count; ++index) {
      const m1_thcm1_fixture_record *record =
          &collections[shard].records[index];
      const double *inputs[2] = {
          record->baseline_input, record->perturbed_input};
      const bool input_invariant_control = record->sensitivity_count == 0;
      double output[M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {0.0};
      double normalizations[2][M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {{0.0}};
      double weighted_states[2][4][M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {{{0.0}}};
      double weighted_fluxes_L[2][M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {{0.0}};
      double weighted_fluxes_R[2][M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {{0.0}};
      char preparation_error[256] = {0};
      int status_ok = 1;
      bool changed_consumed_input = false;
      bool record_has_nonconstant_volumes = false;
      /* Compare the complete 50-field common wire input. This includes the
       * right physical face flux at fields 35..39 and the producer-side
       * physical operands at fields 40..49, so no consumed perturbation can
       * be hidden by a truncated range. */
      for(size_t field = 0;
          field < M1_THCM1_TRANSPORT_INPUT_COUNT;
          ++field)
        changed_consumed_input |= inputs[0][field] != inputs[1][field];
      for(int role = 0; role < 2; ++role) {
        const double *input = inputs[role];
        for(int component = 0;
            component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
          if(input[M1_THCM1_TRANSPORT_PRODUCER_FLUX_L_START + component] !=
                 input[M1_THCM1_TRANSPORT_FLUX_L_START + component] ||
             input[M1_THCM1_TRANSPORT_PRODUCER_FLUX_R_START + component] !=
                 input[M1_THCM1_TRANSPORT_FLUX_R_START + component]) {
            status_ok = 0;
            if(failures == 0)
              fail_test("variable-volume fixture broke the common physical-"
                        "flux transform");
            break;
          }
        if(!status_ok) break;
        const double face_volume =
            input[M1_THCM1_TRANSPORT_FACE_VOLUME_INDEX];
        double cell_volumes[4];
        double state[4][M1_THCM1_TRANSPORT_OUTPUT_COUNT] = {{0.0}};
        double physical_flux_L[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
        double physical_flux_R[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
        for(int cell = 0; cell < 4; ++cell)
          cell_volumes[cell] =
              input[M1_THCM1_TRANSPORT_CELL_VOLUME_START + cell];
        for(int cell = 0; cell < 4; ++cell)
          for(int component = 0;
              component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
            state[cell][component] = input[
                M1_THCM1_TRANSPORT_STATE_START +
                cell * M1_THCM1_TRANSPORT_OUTPUT_COUNT + component];
        for(int component = 0;
            component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component) {
          physical_flux_L[component] =
              input[M1_THCM1_TRANSPORT_FLUX_L_START + component];
          physical_flux_R[component] =
              input[M1_THCM1_TRANSPORT_FLUX_R_START + component];
        }
        for(int cell = 0; cell < 4; ++cell)
          record_has_nonconstant_volumes |= cell_volumes[cell] != face_volume;
        char *error = preparation_error;
        if(!m1_thcm1_transport_prepare_volume_weighted(
               cell_volumes, face_volume, state, physical_flux_L,
               physical_flux_R, weighted_states[role], weighted_fluxes_L[role],
               weighted_fluxes_R[role], error, sizeof(preparation_error)) ||
           !m1_thcm1_transport_compute_weighted_normalization(
               input, normalizations[role])) {
          status_ok = 0;
          if(failures == 0)
            fail_test(error[0] != '\0' ? error :
                      "variable-volume input preparation failed");
          break;
        }
      }
      if(status_ok) {
        /* Standard replay evaluates the current implementation once, using
         * only the baseline prepared operands.  The retained perturbed
         * operands are validated and used only to derive the response bound. */
        const double *input = inputs[0];
        ghl_m1_parameters params = {0};
        if(ghl_m1_initialize(1.0e-8, 1.0e-12, 1.0e-8, 1.0e-6,
                             1.0e-10, 100, 1.0e-10, &params) != ghl_success) {
          status_ok = 0;
        } else {
          /* transport.cc assigns the campaign theta field directly to the
           * limiter control after initialization. */
          params.minmod_theta = input[M1_THCM1_TRANSPORT_THETA_INDEX];
          params.mindiss = input[M1_THCM1_TRANSPORT_MINDISS_INDEX];
          ghl_m1_four_point_transport_diagnostics diagnostics = {0};
          if(m1_thcm1_call_volume_weighted_transport(
                 &params, weighted_states[0], weighted_fluxes_L[0],
                 weighted_fluxes_R[0],
                 input[M1_THCM1_TRANSPORT_SPEED_INDEX],
                 input[M1_THCM1_TRANSPORT_SPEED_INDEX],
                 input[M1_THCM1_TRANSPORT_KAPPA_INDEX],
                 input[M1_THCM1_TRANSPORT_DELTA_INDEX], false, output,
                 &diagnostics) != ghl_success ||
             !isfinite(diagnostics.opacity_suppression) ||
             !isfinite(diagnostics.face_speed))
            status_ok = 0;
        }
      }
      saw_nonconstant_volumes |= record_has_nonconstant_volumes;
      if(input_invariant_control) {
        ++control_records;
        if(changed_consumed_input) {
          status_ok = 0;
          if(failures == 0)
            fail_test("input-invariant variable control changed a consumed "
                      "input");
        }
      } else {
        ++changed_records;
      }
      if(!input_invariant_control && !changed_consumed_input) {
        status_ok = 0;
        if(failures == 0)
          fail_test("variable-volume pair changed no consumed input");
      }
      if(!record_has_nonconstant_volumes) {
        status_ok = 0;
        if(failures == 0)
          fail_test("variable-volume family contains a constant-volume pair");
      }
      m1_thcm1_fixture_comparison_report report = {0};
      char compare_error[256] = {0};
      double computed_normalization[M1_THCM1_TRANSPORT_OUTPUT_COUNT];
      for(int component = 0;
          component < M1_THCM1_TRANSPORT_OUTPUT_COUNT; ++component)
        computed_normalization[component] =
            fmax(normalizations[0][component], normalizations[1][component]);
      if(!status_ok || !m1_thcm1_transport_compare_baseline_response(
             record, collections[shard].policy, computed_normalization, output,
             &report, compare_error, sizeof(compare_error))) {
        ++failures;
        if(failures == 1)
          ghl_error("unit_test_m1_thcm1_blended_rusanov: variable-volume "
                    "comparison failed direction=%s case=%s pair=%s "
                    "component=%zu mode=%s: %s\n",
                    m1_thcm1_transport_direction(record->case_id),
                    record->case_id, record->pair_id, report.component,
                    report.gate != NULL ? report.gate : "unknown",
                    compare_error[0] != '\0' ? compare_error :
                    "prepared evaluation failed");
      }
    }
  if(failures != 0)
    fail_test("one or more variable-volume THC_M1 pairs failed");
  else if(!saw_nonconstant_volumes)
    fail_test("variable-volume fixture family has no varying volumes");
  else
    ghl_info("unit_test_m1_thcm1_blended_rusanov: %zu variable-volume "
             "THC_M1 pairs passed baseline/response-envelope gates "
             "(%zu changed-input, %zu input-invariant control)\n",
             total_records, changed_records, control_records);
  for(size_t shard = 0; shard < shard_count; ++shard)
    m1_thcm1_fixture_free(&collections[shard]);
}

static void check_transport_fixtures(const char *restrict fixture_dir) {
  static const char *const fixture_names[] = {
      "/transport_four_point_d0.dat",
      "/transport_four_point_d1.dat",
      "/transport_four_point_d2.dat"};
  /* Derived from the retained exporter corpus: all 63,504 constant-volume
   * pairs minus 6,120 exact duplicate input pairs. */
  static const size_t expected_shard_counts[] = {19128, 19128, 19128};
  static const size_t expected_record_count = 57384;
  static const size_t shard_count = sizeof(fixture_names) / sizeof(fixture_names[0]);
  if(fixture_dir == NULL) {
    fail_test("transport fixture directory is null");
    return;
  }
  m1_thcm1_fixture_collection collections[3] = {{0}};
  for(size_t shard = 0; shard < shard_count; ++shard) {
    const size_t directory_length = strlen(fixture_dir);
    const size_t fixture_name_length = strlen(fixture_names[shard]);
    char *path = (char *)malloc(directory_length + fixture_name_length + 1);
    if(path == NULL) {
      fail_test("transport fixture path allocation failed");
      for(size_t cleanup = 0; cleanup < shard_count; ++cleanup)
        m1_thcm1_fixture_free(&collections[cleanup]);
      return;
    }
    memcpy(path, fixture_dir, directory_length);
    memcpy(path + directory_length, fixture_names[shard], fixture_name_length + 1);
    char error[256] = {0};
    if(!m1_thcm1_fixture_load(path, "neutrino_four_point_transport_flux", 50, 5,
                              &collections[shard], error, sizeof(error))) {
      fail_test(error[0] != '\0' ? error : "transport fixture load failed");
      free(path);
      for(size_t cleanup = 0; cleanup < shard_count; ++cleanup)
        m1_thcm1_fixture_free(&collections[cleanup]);
      return;
    }
    free(path);
    if(strcmp(collections[shard].policy,
             "strict_relative_2e-12_propagated_response_v1") != 0 ||
       collections[shard].record_count != expected_shard_counts[shard]) {
      fail_test("transport fixture shard policy or exact record count mismatch");
      for(size_t cleanup = 0; cleanup < shard_count; ++cleanup)
        m1_thcm1_fixture_free(&collections[cleanup]);
      return;
    }
  }
  size_t total_records = 0;
  for(size_t shard = 0; shard < shard_count; ++shard)
    total_records += collections[shard].record_count;
  if(total_records != expected_record_count) {
    fail_test("transport fixture exact admitted-corpus count mismatch");
    for(size_t cleanup = 0; cleanup < shard_count; ++cleanup)
      m1_thcm1_fixture_free(&collections[cleanup]);
    return;
  }
  const char *directions[] = {"__d0__", "__d1__", "__d2__"};
  const char *profiles[] = {"__constant__", "__monotone__", "__sawtooth__"};
  const char *thetas[] = {"__th0.0__", "__th1.0__", "__th2.0__"};
  const char *mindiss[] = {"__md0.0__", "__md1.0__"};
  const char *resolutions[] = {"__r32__", "__r64__", "__r128__", "__r256__"};
  const char *opacities[] = {"__zero__", "__transition__", "__packet__"};
  const char *const *categories[] = {
      directions, profiles, thetas, mindiss, resolutions, opacities};
  const size_t category_counts[] = {3, 3, 3, 2, 4, 3};
  static const size_t expected_category_counts[6][4] = {
      {19128, 19128, 19128, 0},
      {16056, 20664, 20664, 0},
      {19128, 19128, 19128, 0},
      {28692, 28692, 0, 0},
      {14346, 14346, 14346, 14346},
      {18792, 18792, 19800, 0}};
  size_t observed_category_counts[6][4] = {{0}};
  for(size_t shard = 0; shard < shard_count; ++shard)
    for(size_t index = 0; index < collections[shard].record_count; ++index) {
      const char *case_id = collections[shard].records[index].case_id;
      for(size_t category = 0; category < 6; ++category)
        for(size_t value = 0; value < category_counts[category]; ++value)
          if(strstr(case_id, categories[category][value]) != NULL)
            ++observed_category_counts[category][value];
    }
  for(size_t category = 0; category < 6; ++category)
    for(size_t value = 0; value < category_counts[category]; ++value)
      if(observed_category_counts[category][value] !=
         expected_category_counts[category][value]) {
        fail_test("transport fixture exact branch inventory mismatch");
        for(size_t cleanup = 0; cleanup < shard_count; ++cleanup)
          m1_thcm1_fixture_free(&collections[cleanup]);
        return;
      }

  size_t failures = 0;
  bool saw_nonunit_constant_volume = false;
  for(size_t shard = 0; shard < shard_count; ++shard)
    for(size_t index = 0; index < collections[shard].record_count; ++index) {
    const m1_thcm1_fixture_record *record = &collections[shard].records[index];
    const double *base = record->baseline_input;
    const double *pert = record->perturbed_input;
    for(int role = 0; role < 2; ++role) {
      const double *input = role == 0 ? base : pert;
      if(!(input[0] > 0.0) ||
         input[6] != input[0] || input[7] != input[0] ||
         input[8] != input[0] || input[9] != input[0]) {
        fail_test("constant-volume transport record has inconsistent cell volumes");
        for(size_t cleanup = 0; cleanup < shard_count; ++cleanup)
          m1_thcm1_fixture_free(&collections[cleanup]);
        return;
      }
      if(input[0] != 1.0) saw_nonunit_constant_volume = true;
    }
    double computed_normalization[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    for(int role = 0; role < 2; ++role) {
      const double *input = role == 0 ? base : pert;
      const double speed = fabs(input[1]);
      for(int component = 0; component < 5; ++component) {
        const double scale = fmax(
            fmax(fabs(input[30 + component]), fabs(input[35 + component])),
            fmax(speed * fabs(input[10 + component]),
                 speed * fabs(input[15 + component])));
        computed_normalization[component] = fmax(
            computed_normalization[component], fmax(1.0e-300, scale));
      }
    }

    double output[5] = {0.0};
    int status_ok = 1;
    const double *input = base;
    ghl_m1_parameters params = {0};
    if(ghl_m1_initialize(1.0e-8, 1.0e-12, 1.0e-8, 1.0e-6,
                         1.0e-10, 100, 1.0e-10, &params) != ghl_success) {
      status_ok = 0;
    } else {
      /* transport.cc assigns the campaign theta field directly to the
       * limiter control after initialization. */
      params.minmod_theta = input[4];
      params.mindiss = input[5];
      ghl_metric_quantities metric = {0};
      const double spatial_metric = pow(input[0], 2.0 / 3.0);
      ghl_initialize_metric(1.0, 0.0, 0.0, 0.0,
                            spatial_metric, 0.0, 0.0,
                            spatial_metric, 0.0, spatial_metric, &metric);
      double state[4][ghl_m1_neutrino_transport_component_count] = {{0.0}};
      for(int cell = 0; cell < 4; ++cell)
        for(int component = 0; component < 5; ++component)
          state[cell][component] = input[10 + 5 * cell + component];
      const double *physical_left = &input[30];
      const double *physical_right = &input[35];
      ghl_m1_four_point_transport_diagnostics diagnostics = {0};
      if(ghl_m1_compute_neutrino_four_point_transport_flux(
             &params, &metric, state, physical_left, physical_right,
             input[1], input[1], input[2], input[3], false, output,
             &diagnostics) != ghl_success) {
        status_ok = 0;
      }
      if(!isfinite(diagnostics.opacity_suppression) ||
         !isfinite(diagnostics.face_speed)) {
        status_ok = 0;
      }
    }
    m1_thcm1_fixture_comparison_report report = {0};
    char compare_error[256] = {0};
    if(!status_ok || !m1_thcm1_transport_compare_baseline_response(
           record, collections[shard].policy, computed_normalization, output,
           &report, compare_error, sizeof(compare_error))) {
      ++failures;
      if(failures == 1)
        fail_test(compare_error[0] != '\0' ? compare_error :
                  "transport fixture GRHayL evaluation failed");
    }
    }
  if(failures != 0)
    fail_test("one or more stored transport pairs failed");
  else if(!saw_nonunit_constant_volume)
    fail_test("transport fixture lacks a nonunit constant face volume");
  else
    ghl_info("unit_test_m1_thcm1_blended_rusanov: %zu strict transport pairs passed\n",
             total_records);
  for(size_t shard = 0; shard < shard_count; ++shard)
    m1_thcm1_fixture_free(&collections[shard]);
}

static void check_four_point_branch_boundaries(
      ghl_m1_parameters *restrict params) {
  double phi = 41.0;
  bool sawtooth = true;
  if(ghl_m1_compute_four_point_flux_limiter(
         NULL, 1.0, 1.0, 1.0, &phi, &sawtooth) != ghl_error_m1_null_pointer ||
     ghl_m1_compute_four_point_flux_limiter(
         params, 1.0, 1.0, 1.0, NULL, &sawtooth) !=
         ghl_error_m1_null_pointer ||
     ghl_m1_compute_four_point_flux_limiter(
         params, 1.0, 1.0, 1.0, &phi, NULL) != ghl_error_m1_null_pointer)
    fail_test("limiter NULL arguments were not rejected");

  params->minmod_theta = NAN;
  if(ghl_m1_compute_four_point_flux_limiter(
         params, 1.0, 1.0, 1.0, &phi, &sawtooth) !=
         ghl_error_m1_invalid_state)
    fail_test("nonfinite limiter theta was accepted");
  params->minmod_theta = -0.1;
  if(ghl_m1_compute_four_point_flux_limiter(
         params, 1.0, 1.0, 1.0, &phi, &sawtooth) !=
         ghl_error_m1_invalid_state)
    fail_test("negative limiter theta was accepted");
  params->minmod_theta = 2.1;
  if(ghl_m1_compute_four_point_flux_limiter(
         params, 1.0, 1.0, 1.0, &phi, &sawtooth) !=
         ghl_error_m1_invalid_state)
    fail_test("oversized limiter theta was accepted");
  params->minmod_theta = 1.0;
  if(ghl_m1_compute_four_point_flux_limiter(
         params, NAN, 1.0, 1.0, &phi, &sawtooth) !=
         ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_flux_limiter(
         params, 1.0, NAN, 1.0, &phi, &sawtooth) !=
         ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_flux_limiter(
         params, 1.0, 1.0, NAN, &phi, &sawtooth) !=
         ghl_error_m1_invalid_state)
    fail_test("nonfinite limiter increments were accepted");

  if(ghl_m1_compute_four_point_flux_limiter(
         params, -1.0, -1.0, -1.0, &phi, &sawtooth) != ghl_success ||
     !m1_nearly_equal(phi, 1.0, 0.0, 0.0) || sawtooth)
    fail_test("negative monotone limiter branch failed");
  if(ghl_m1_compute_four_point_flux_limiter(
         params, 1.0, -1.0, 1.0, &phi, &sawtooth) != ghl_success ||
     phi != 0.0 || !sawtooth)
    fail_test("positive-to-negative sawtooth branch failed");
  if(ghl_m1_compute_four_point_flux_limiter(
         params, -1.0, 1.0, -1.0, &phi, &sawtooth) != ghl_success ||
     phi != 0.0 || !sawtooth)
    fail_test("negative-to-positive sawtooth branch failed");

  /* Zero and both signs are the sign predicates' complete input classes. */
  const double signs[] = {-1.0, 0.0, 1.0};
  for(size_t i = 0; i < sizeof(signs)/sizeof(signs[0]); ++i)
    for(size_t j = 0; j < sizeof(signs)/sizeof(signs[0]); ++j)
      for(size_t k = 0; k < sizeof(signs)/sizeof(signs[0]); ++k) {
        const bool monotone = signs[i] == signs[j] && signs[j] == signs[k]
                           && signs[j] != 0.0;
        const bool alternating = signs[i] == -signs[j]
                              && signs[k] == -signs[j] && signs[j] != 0.0;
        if(ghl_m1_compute_four_point_flux_limiter(params, signs[i], signs[j],
            signs[k], &phi, &sawtooth) != ghl_success ||
           phi != (monotone ? 1.0 : 0.0) || sawtooth != alternating)
          fail_test("limiter sign-class oracle mismatch");
      }

  double A = 42.0;
  if(ghl_m1_compute_four_point_opacity_suppression(
         NULL, 1.0, 1.0, &A) != ghl_error_m1_null_pointer ||
     ghl_m1_compute_four_point_opacity_suppression(
         params, 1.0, 1.0, NULL) != ghl_error_m1_null_pointer)
    fail_test("opacity NULL arguments were not rejected");
  if(ghl_m1_compute_four_point_opacity_suppression(
         params, NAN, 1.0, &A) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_opacity_suppression(
         params, -1.0, 1.0, &A) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_opacity_suppression(
         params, 1.0, NAN, &A) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_opacity_suppression(
         params, 1.0, 0.0, &A) != ghl_error_m1_invalid_state)
    fail_test("invalid opacity geometry was accepted");
  params->mindiss = NAN;
  if(ghl_m1_compute_four_point_opacity_suppression(
         params, 1.0, 1.0, &A) != ghl_error_m1_invalid_state)
    fail_test("nonfinite minimum dissipation was accepted");
  params->mindiss = -0.1;
  if(ghl_m1_compute_four_point_opacity_suppression(
         params, 1.0, 1.0, &A) != ghl_error_m1_invalid_state)
    fail_test("negative minimum dissipation was accepted");
  params->mindiss = 1.1;
  if(ghl_m1_compute_four_point_opacity_suppression(
         params, 1.0, 1.0, &A) != ghl_error_m1_invalid_state)
    fail_test("oversized minimum dissipation was accepted");
  params->mindiss = 0.0;

  double blended = 43.0;
  if(ghl_m1_compute_four_point_blended_flux(
         NAN, 1.0, 0.5, false, 0.5, &blended) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_blended_flux(
         1.0, NAN, 0.5, false, 0.5, &blended) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_blended_flux(
         1.0, 1.0, NAN, false, 0.5, &blended) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_blended_flux(
         1.0, 1.0, -0.1, false, 0.5, &blended) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_blended_flux(
         1.0, 1.0, 0.5, false, NAN, &blended) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_blended_flux(
         1.0, 1.0, 0.5, false, -0.1, &blended) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_blended_flux(
         1.0, 1.0, 0.5, false, 1.1, &blended) != ghl_error_m1_invalid_state ||
     ghl_m1_compute_four_point_blended_flux(
         1.0, 1.0, 0.5, false, 0.5, NULL) != ghl_error_m1_null_pointer)
    fail_test("invalid blended-flux arguments were accepted");
  if(ghl_m1_compute_four_point_blended_flux(
         DBL_MAX, -DBL_MAX, 0.0, true, 1.0, &blended) !=
         ghl_error_m1_invalid_state || blended != 43.0)
    fail_test("overflowing blended flux was published");

  double state_stencil[4][ghl_m1_neutrino_transport_component_count] = {{0.0}};
  double physical_flux_L[ghl_m1_neutrino_transport_component_count] = {0.0};
  double physical_flux_R[ghl_m1_neutrino_transport_component_count] = {0.0};
  double flux_tilde[ghl_m1_neutrino_transport_component_count] = {
      101.0, 102.0, 103.0, 104.0, 105.0};
  const double flux_before[ghl_m1_neutrino_transport_component_count] = {
      101.0, 102.0, 103.0, 104.0, 105.0};
  const double speeds[4] = {NAN, -1.0, 0.5, 0.5};
  const double right_speeds[4] = {0.5, 0.5, NAN, -1.0};
  /* Finite low/high fluxes can still overflow when their difference is
   * rounded: these binary64 operands make low=high-DBL_MAX finite but
   * high-low infinite. Check propagation from the final blend, independently
   * of earlier limiter/Rusanov rejection. */
  const double high = 0x1.026650c8e9dc3p+1022;
  const double low = high - DBL_MAX;
  if(!isfinite(low) || isfinite(high - low))
    fail_test("blend-overflow witness does not have its required arithmetic");
  state_stencil[2][0] = state_stencil[3][0] = DBL_MAX;
  physical_flux_L[0] = 2.0 * high;
  if(m1_thcm1_call_volume_weighted_transport(params, state_stencil,
      physical_flux_L, physical_flux_R, 2.0, 2.0, 0.0, 1.0, false,
      flux_tilde, NULL) != ghl_error_m1_invalid_state ||
     memcmp(flux_tilde, flux_before, sizeof(flux_tilde)) != 0)
    fail_test("prepared final-blend overflow was not transactional");
  memset(state_stencil, 0, sizeof(state_stencil));
  physical_flux_L[0] = 0.0;
  /* Reach the core's policy and helper-error propagation via the prepared
   * API, whose wrapper does not precheck the pointwise policy. */
  if(m1_thcm1_call_volume_weighted_transport(params, state_stencil,
      physical_flux_L, physical_flux_R, 0.5, 0.5, 0.5, 1.0, true,
      flux_tilde, NULL) != ghl_error_m1_incompatible_transport_policy ||
     memcmp(flux_tilde, flux_before, sizeof(flux_tilde)) != 0)
    fail_test("prepared transport policy failure was not transactional");
  const double saved_mindiss = params->mindiss;
  params->mindiss = NAN;
  if(m1_thcm1_call_volume_weighted_transport(params, state_stencil,
      physical_flux_L, physical_flux_R, 0.5, 0.5, 0.5, 1.0, false,
      flux_tilde, NULL) != ghl_error_m1_invalid_state ||
     memcmp(flux_tilde, flux_before, sizeof(flux_tilde)) != 0)
    fail_test("prepared opacity-helper error was not transactional");
  params->mindiss = saved_mindiss;
  state_stencil[0][0] = -DBL_MAX;
  state_stencil[1][0] = state_stencil[2][0] = state_stencil[3][0] = DBL_MAX;
  if(m1_thcm1_call_volume_weighted_transport(params, state_stencil,
      physical_flux_L, physical_flux_R, 0.5, 0.5, 0.5, 1.0, false,
      flux_tilde, NULL) != ghl_error_m1_invalid_state ||
     memcmp(flux_tilde, flux_before, sizeof(flux_tilde)) != 0)
    fail_test("prepared limiter-difference overflow was not transactional");
  for(int cell = 0; cell < 4; ++cell) state_stencil[cell][0] = 0.0;
  for(int case_index = 0; case_index < 4; ++case_index) {
    memcpy(flux_tilde, flux_before, sizeof(flux_tilde));
    if(m1_thcm1_call_volume_weighted_transport(
           params, state_stencil, physical_flux_L, physical_flux_R,
           speeds[case_index], right_speeds[case_index], 0.5, 1.0, false,
           flux_tilde, NULL) != ghl_error_m1_invalid_state ||
       memcmp(flux_tilde, flux_before, sizeof(flux_tilde)) != 0)
      fail_test("invalid four-point speed was not transactional");
  }
  const double kappas[2] = {NAN, -1.0};
  for(int case_index = 0; case_index < 2; ++case_index) {
    memcpy(flux_tilde, flux_before, sizeof(flux_tilde));
    if(m1_thcm1_call_volume_weighted_transport(
           params, state_stencil, physical_flux_L, physical_flux_R,
           0.5, 0.5, kappas[case_index], 1.0, false, flux_tilde, NULL) !=
           ghl_error_m1_invalid_state ||
       memcmp(flux_tilde, flux_before, sizeof(flux_tilde)) != 0)
      fail_test("invalid four-point opacity was not transactional");
  }
  const double deltas[2] = {NAN, 0.0};
  for(int case_index = 0; case_index < 2; ++case_index) {
    memcpy(flux_tilde, flux_before, sizeof(flux_tilde));
    if(m1_thcm1_call_volume_weighted_transport(
           params, state_stencil, physical_flux_L, physical_flux_R,
           0.5, 0.5, 0.5, deltas[case_index], false, flux_tilde, NULL) !=
           ghl_error_m1_invalid_state ||
       memcmp(flux_tilde, flux_before, sizeof(flux_tilde)) != 0)
      fail_test("invalid four-point spacing was not transactional");
  }

  if(m1_thcm1_call_volume_weighted_transport(
         NULL, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, flux_tilde, NULL) !=
         ghl_error_m1_null_pointer ||
     m1_thcm1_call_volume_weighted_transport(
         params, NULL, physical_flux_L, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, flux_tilde, NULL) !=
         ghl_error_m1_null_pointer ||
     m1_thcm1_call_volume_weighted_transport(
         params, state_stencil, NULL, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, flux_tilde, NULL) !=
         ghl_error_m1_null_pointer ||
     m1_thcm1_call_volume_weighted_transport(
         params, state_stencil, physical_flux_L, NULL,
         0.5, 0.5, 0.5, 1.0, false, flux_tilde, NULL) !=
         ghl_error_m1_null_pointer ||
     m1_thcm1_call_volume_weighted_transport(
         params, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, NULL, NULL) !=
         ghl_error_m1_null_pointer)
    fail_test("prepared four-point NULL arguments were not rejected");

  for(int cell = 0; cell < 4; ++cell) {
    state_stencil[cell][0] = NAN;
    memcpy(flux_tilde, flux_before, sizeof(flux_tilde));
    if(m1_thcm1_call_volume_weighted_transport(
           params, state_stencil, physical_flux_L, physical_flux_R,
           0.5, 0.5, 0.5, 1.0, false, flux_tilde, NULL) !=
           ghl_error_m1_invalid_state ||
       memcmp(flux_tilde, flux_before, sizeof(flux_tilde)) != 0)
      fail_test("invalid prepared stencil state was accepted");
    state_stencil[cell][0] = 0.0;
  }
  physical_flux_L[0] = NAN;
  if(m1_thcm1_call_volume_weighted_transport(
         params, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, flux_tilde, NULL) !=
         ghl_error_m1_invalid_state)
    fail_test("invalid left physical flux was accepted");
  physical_flux_L[0] = 0.0;
  physical_flux_R[0] = NAN;
  if(m1_thcm1_call_volume_weighted_transport(
         params, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, flux_tilde, NULL) !=
         ghl_error_m1_invalid_state)
    fail_test("invalid right physical flux was accepted");
  physical_flux_R[0] = 0.0;

  state_stencil[1][0] = DBL_MAX;
  state_stencil[2][0] = -DBL_MAX;
  memcpy(flux_tilde, flux_before, sizeof(flux_tilde));
  if(m1_thcm1_call_volume_weighted_transport(
         params, state_stencil, physical_flux_L, physical_flux_R,
         1.0, 1.0, 0.5, 1.0, false, flux_tilde, NULL) !=
         ghl_error_m1_invalid_state ||
     memcmp(flux_tilde, flux_before, sizeof(flux_tilde)) != 0)
    fail_test("overflowing low-order flux was published");
  state_stencil[1][0] = 0.0;
  state_stencil[2][0] = 0.0;

  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  double pointwise_flux[ghl_m1_neutrino_transport_component_count] = {
      201.0, 202.0, 203.0, 204.0, 205.0};
  const double pointwise_flux_before[
      ghl_m1_neutrino_transport_component_count] = {
      201.0, 202.0, 203.0, 204.0, 205.0};
  if(ghl_m1_compute_neutrino_four_point_transport_flux(
         NULL, &metric, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, pointwise_flux, NULL) !=
         ghl_error_m1_null_pointer ||
     ghl_m1_compute_neutrino_four_point_transport_flux(
         params, NULL, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, pointwise_flux, NULL) !=
         ghl_error_m1_null_pointer ||
     ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &metric, NULL, physical_flux_L, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, pointwise_flux, NULL) !=
         ghl_error_m1_null_pointer ||
     ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &metric, state_stencil, NULL, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, pointwise_flux, NULL) !=
         ghl_error_m1_null_pointer ||
     ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &metric, state_stencil, physical_flux_L, NULL,
         0.5, 0.5, 0.5, 1.0, false, pointwise_flux, NULL) !=
         ghl_error_m1_null_pointer ||
     ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &metric, state_stencil, physical_flux_L, physical_flux_R,
         0.5, 0.5, 0.5, 1.0, false, NULL, NULL) !=
         ghl_error_m1_null_pointer)
    fail_test("pointwise four-point NULL arguments were not rejected");

  ghl_metric_quantities invalid_metric = metric;
  invalid_metric.gammaDD[0][0] = -1.0;
  if(ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &invalid_metric, state_stencil, physical_flux_L,
         physical_flux_R, 0.5, 0.5, 0.5, 1.0, false, pointwise_flux, NULL) !=
         ghl_error_m1_invalid_metric)
    fail_test("non-SPD four-point metric was accepted");
  invalid_metric = metric;
  invalid_metric.sqrt_detgamma = NAN;
  if(ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &invalid_metric, state_stencil, physical_flux_L,
         physical_flux_R, 0.5, 0.5, 0.5, 1.0, false, pointwise_flux, NULL) !=
         ghl_error_m1_invalid_metric)
    fail_test("nonfinite four-point volume factor was accepted");

  physical_flux_L[0] = 2.0;
  physical_flux_R[0] = 2.0;
  /* Keep the metric coherent while making the finite volume-weighted
   * candidate large enough that densitization, rather than metric validation,
   * reaches the publication guard. */
  ghl_metric_quantities large_volume_metric = metric;
  large_volume_metric.gammaDD[0][0] = 16.0;
  large_volume_metric.gammaUU[0][0] = 1.0 / 16.0;
  large_volume_metric.detgamma = 16.0;
  large_volume_metric.sqrt_detgamma = 4.0;
  physical_flux_L[0] = DBL_MAX;
  physical_flux_R[0] = 0.0;
  if(ghl_m1_compute_neutrino_four_point_transport_flux(
         params, &large_volume_metric, state_stencil, physical_flux_L,
         physical_flux_R, 0.0, 0.0, 0.5, 1.0, false, pointwise_flux, NULL) !=
         ghl_error_m1_invalid_state ||
     memcmp(pointwise_flux, pointwise_flux_before, sizeof(pointwise_flux)) != 0)
    fail_test("overflowing densitized four-point flux was published");
}

int main(int argc, char **argv) {
  const char *fixture_dir = "Unit_Tests/data/m1_thcm1";
  if(argc == 3 && strcmp(argv[1], "--fixture-dir") == 0)
    fixture_dir = argv[2];
  else if(argc != 1) {
    fail_test("usage: [--fixture-dir PATH]");
    return 1;
  }
  ghl_m1_parameters params = {0};
  if(ghl_m1_initialize(
         1.0e-10, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-12,
         20, 1.0e-10, &params) != ghl_success)
    fail_test("M1 initialization failed");

  double phi = NAN;
  bool sawtooth = false;
  if(ghl_m1_compute_four_point_flux_limiter(
         &params, 1.0, 1.0, 1.0, &phi, &sawtooth) != ghl_success ||
     phi != 1.0 || sawtooth)
    fail_test("smooth monotone limiter branch failed");
  if(ghl_m1_compute_four_point_flux_limiter(
         &params, 0.2, 0.5, 0.4, &phi, &sawtooth) != ghl_success ||
     !m1_nearly_equal(phi, 0.4, 2.0e-13, 2.0e-14) || sawtooth)
    fail_test("limited monotone limiter branch failed");
  if(ghl_m1_compute_four_point_flux_limiter(
         &params, -1.0, 1.0, -0.5, &phi, &sawtooth) != ghl_success ||
     phi != 0.0 || !sawtooth)
    fail_test("sawtooth limiter branch failed");
  if(ghl_m1_compute_four_point_flux_limiter(
         &params, 0.0, 1.0, 0.5, &phi, &sawtooth) != ghl_success ||
     phi != 0.0 || sawtooth)
    fail_test("zero-product limiter branch failed");

  double A = NAN;
  if(ghl_m1_compute_four_point_opacity_suppression(
         &params, 0.5, 1.0, &A) != ghl_success || A != 1.0)
    fail_test("thin opacity branch failed");
  if(ghl_m1_compute_four_point_opacity_suppression(
         &params, 4.0, 1.0, &A) != ghl_success || A != 0.25)
    fail_test("thick opacity branch failed");
  params.mindiss = 0.2;
  if(ghl_m1_compute_four_point_opacity_suppression(
         &params, 100.0, 1.0, &A) != ghl_success || A != 0.2)
    fail_test("minimum dissipation branch failed");
  params.mindiss = 0.0;

  double blended = NAN;
  if(ghl_m1_compute_four_point_blended_flux(
         2.0, 1.0, 0.25, false, 0.4, &blended) != ghl_success)
    fail_test("ordinary blend failed");
  check_close(blended, 2.0 - 0.4 * 0.75, "ordinary blend formula failed");
  if(ghl_m1_compute_four_point_blended_flux(
         2.0, 1.0, 0.25, true, 0.4, &blended) != ghl_success)
    fail_test("sawtooth blend failed");
  check_close(blended, 2.0 - 0.75, "sawtooth blend formula failed");
  if(ghl_m1_compute_four_point_blended_flux(
         2.0, 1.0, 1.0, false, 0.0, &blended) != ghl_success ||
     blended != 2.0)
    fail_test("unit limiter blend failed");

  /* Each public stage validates before publication. */
  phi = 41.0;
  sawtooth = true;
  if(ghl_m1_compute_four_point_flux_limiter(
         &params, NAN, 1.0, 1.0, &phi, &sawtooth) !=
         ghl_error_m1_invalid_state || phi != 41.0 || !sawtooth)
    fail_test("invalid limiter input was not transactional");
  A = 42.0;
  if(ghl_m1_compute_four_point_opacity_suppression(
         &params, 1.0, 0.0, &A) != ghl_error_m1_invalid_state || A != 42.0)
    fail_test("invalid opacity input was not transactional");
  blended = 43.0;
  if(ghl_m1_compute_four_point_blended_flux(
         2.0, 1.0, 1.1, false, 0.5, &blended) !=
         ghl_error_m1_invalid_state || blended != 43.0)
    fail_test("invalid blend input was not transactional");

  check_four_point_branch_boundaries(&params);
  check_full_four_point_operator(&params);
  check_prepared_transport_local_contract();
  check_transport_fixtures(fixture_dir);
  check_variable_transport_fixtures(fixture_dir);

  ghl_info("unit_test_m1_thcm1_blended_rusanov: "
           "limiter/opacity/blend/transport branches passed\n");
  return 0;
}
