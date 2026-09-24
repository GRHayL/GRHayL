#include <math.h>
#include <stdio.h>
#include <string.h>

#include "m1_test_utils.h"
#include "m1_thcm1_fixture_utils.h"

/*
 * Direct property and boundary coverage for the optional M1 diffusion
 * helpers. The test checks Jthick, harmonic diffusion coefficients, metric
 * face distances, and transactional rejection of invalid inputs.
 */

static void fail_test(const char *message) {
  ghl_error("unit_test_m1_diffusion_flux: %s\n", message);
}

static void check_close(const double actual, const double expected, const char *label) {
  if(!m1_nearly_equal(actual, expected, 3.0e-12, 2.0e-13)) {
    fail_test(label);
  }
}

static double expected_diffusion_flux(
      const ghl_m1_parameters *restrict params,
      const ghl_metric_quantities *restrict metric,
      const int direction,
      const double hll_flux_tildeE,
      const double E_star,
      const double J_left,
      const double J_right,
      const double gradJ[3],
      const double W,
      const double V[3],
      const double chi_tr,
      const double D,
      const double delta_l,
      double *restrict a_out) {

  const double tau = chi_tr * delta_l;
  const double a = tanh(1.0 / tau);
  const double J_face = 0.5 * (J_left + J_right);
  const double W2 = W * W;
  const double F_adv = (4.0 / 3.0) * W2 * V[direction] * J_face;
  double gamma_grad = 0.0;
  double VdotGrad = 0.0;
  for(int i = 0; i < 3; ++i) {
    gamma_grad += metric->gammaUU[direction][i] * gradJ[i];
    VdotGrad += V[i] * gradJ[i];
  }
  const double F_diff = W * D * (gamma_grad + V[direction] * VdotGrad);
  const double F_asym = F_adv - F_diff;
  const double scalar_shift = fmax(E_star, params->E_floor);
  const double flux_asym
        = metric->sqrt_detgamma
          * (metric->lapse * F_asym - metric->betaU[direction] * scalar_shift);
  if(a_out != NULL) {
    *a_out = a;
  }
  return a * hll_flux_tildeE + (1.0 - a) * flux_asym;
}

static void setup_primitives(
      ghl_primitive_quantities *restrict prims,
      const double vx,
      const double vy,
      const double vz) {
  *prims = (ghl_primitive_quantities){ 0 };
  prims->rho = 1.0;
  prims->press = 0.4;
  prims->eps = 0.2;
  prims->u0 = 1.0;
  prims->vU[0] = vx;
  prims->vU[1] = vy;
  prims->vU[2] = vz;
}

static void check_diffusion_transport_boundaries(
      ghl_m1_parameters *restrict params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict state) {
  double Jthick = 31.0;
  bool valid = true;
  if(ghl_m1_compute_Jthick(NULL, metric, prims, state, &Jthick, &valid)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_Jthick(params, NULL, prims, state, &Jthick, &valid)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_Jthick(params, metric, NULL, state, &Jthick, &valid)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_Jthick(params, metric, prims, NULL, &Jthick, &valid)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_Jthick(params, metric, prims, state, NULL, &valid)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_Jthick(params, metric, prims, state, &Jthick, NULL)
              != ghl_error_m1_null_pointer) {
    fail_test("Jthick NULL arguments were not rejected");
  }

  ghl_primitive_quantities high_energy_prims = *prims;
  high_energy_prims.vU[0] = 0.9;
  ghl_m1_rad_state high_energy = { .E = DBL_MAX, .F = { 0.0, 0.0, 0.0 } };
  if(ghl_m1_compute_Jthick(
           params, metric, &high_energy_prims, &high_energy, &Jthick, &valid)
           != ghl_success
     || isfinite(Jthick) || valid) {
    fail_test("finite high-energy Jthick overflow was not reported");
  }
  ghl_m1_rad_state invalid_state = *state;
  invalid_state.F[0] = 2.0 * invalid_state.E;
  if(ghl_m1_compute_Jthick(params, metric, prims, &invalid_state, &Jthick, &valid)
     != ghl_error_m1_invalid_state) {
    fail_test("unrealizable Jthick state was accepted");
  }
  ghl_primitive_quantities invalid_prims = *prims;
  invalid_prims.vU[0] = NAN;
  if(ghl_m1_compute_Jthick(params, metric, &invalid_prims, state, &Jthick, &valid)
     != ghl_error_m1_invalid_state) {
    fail_test("invalid Jthick primitive velocity was accepted");
  }

  double raw_minus = 17.0, raw_plus = 18.0;
  if(ghl_m1_compute_raw_lightcone_speeds(NULL, ghl_m1_dirn0, &raw_minus, &raw_plus)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_raw_lightcone_speeds(
              metric, (ghl_m1_direction_t)3, &raw_minus, &raw_plus)
              != ghl_error_m1_invalid_state
     || ghl_m1_compute_raw_lightcone_speeds(metric, ghl_m1_dirn0, NULL, &raw_plus)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_raw_lightcone_speeds(metric, ghl_m1_dirn0, &raw_minus, NULL)
              != ghl_error_m1_null_pointer) {
    fail_test("raw light-cone speed boundary was not rejected");
  }
  ghl_metric_quantities invalid_metric = *metric;
  invalid_metric.gammaDD[0][0] = -1.0;
  if(ghl_m1_compute_raw_lightcone_speeds(
           &invalid_metric, ghl_m1_dirn0, &raw_minus, &raw_plus)
     != ghl_error_m1_invalid_metric) {
    fail_test("non-SPD raw-speed metric was accepted");
  }
  invalid_metric = *metric;
  invalid_metric.lapse = NAN;
  if(ghl_m1_compute_raw_lightcone_speeds(
           &invalid_metric, ghl_m1_dirn0, &raw_minus, &raw_plus)
     != ghl_error_m1_invalid_metric) {
    fail_test("nonfinite raw-speed scale was accepted");
  }
  invalid_metric.lapse = 0.0;
  if(ghl_m1_compute_raw_lightcone_speeds(
           &invalid_metric, ghl_m1_dirn0, &raw_minus, &raw_plus)
     != ghl_error_m1_invalid_metric) {
    fail_test("zero raw-speed scale was accepted");
  }
  invalid_metric = *metric;
  invalid_metric.lapse = DBL_MAX;
  invalid_metric.betaU[0] = DBL_MAX;
  if(ghl_m1_compute_raw_lightcone_speeds(
           &invalid_metric, ghl_m1_dirn0, &raw_minus, &raw_plus)
     != ghl_error_m1_invalid_metric) {
    fail_test("overflowing raw light-cone speed was accepted");
  }
  invalid_metric.betaU[0] = -DBL_MAX;
  if(ghl_m1_compute_raw_lightcone_speeds(
           &invalid_metric, ghl_m1_dirn0, &raw_minus, &raw_plus)
     != ghl_error_m1_invalid_metric) {
    fail_test("positive raw light-cone speed overflow was accepted");
  }
  invalid_metric = *metric;
  invalid_metric.lapse = DBL_MAX;
  invalid_metric.gammaDD[0][0] = 1.0e-100;
  invalid_metric.gammaUU[0][0] = 1.0e100;
  invalid_metric.detgamma = 1.0e-100;
  invalid_metric.sqrt_detgamma = 1.0e-50;
  if(ghl_m1_compute_raw_lightcone_speeds(
           &invalid_metric, ghl_m1_dirn0, &raw_minus, &raw_plus)
     != ghl_error_m1_invalid_metric) {
    fail_test("overflowing raw light-cone scale was accepted");
  }

  double clipped_minus = 41.0, clipped_plus = 42.0;
  if(ghl_m1_clip_hll_speeds(NAN, 1.0, &clipped_minus, &clipped_plus)
           != ghl_error_m1_invalid_state
     || ghl_m1_clip_hll_speeds(1.0, NAN, &clipped_minus, &clipped_plus)
              != ghl_error_m1_invalid_state
     || ghl_m1_clip_hll_speeds(2.0, 1.0, &clipped_minus, &clipped_plus)
              != ghl_error_m1_invalid_state
     || ghl_m1_clip_hll_speeds(-2.0, -1.0, &clipped_minus, &clipped_plus) != ghl_success
     || clipped_minus != -2.0 || clipped_plus != 0.0
     || ghl_m1_clip_hll_speeds(1.0, 2.0, &clipped_minus, &clipped_plus) != ghl_success
     || clipped_minus != 0.0 || clipped_plus != 2.0
     || ghl_m1_clip_hll_speeds(1.0, 2.0, NULL, &clipped_plus)
              != ghl_error_m1_null_pointer
     || ghl_m1_clip_hll_speeds(1.0, 2.0, &clipped_minus, NULL)
              != ghl_error_m1_null_pointer) {
    fail_test("HLL speed clipping boundary failed");
  }

  double speed_minus = 51.0, speed_plus = 52.0;
  if(ghl_m1_compute_wavespeeds(NULL, ghl_m1_dirn0, &speed_minus, &speed_plus)
           != ghl_error_m1_null_pointer
     || ghl_m1_compute_wavespeeds(metric, ghl_m1_dirn0, NULL, &speed_plus)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_wavespeeds(metric, ghl_m1_dirn0, &speed_minus, NULL)
              != ghl_error_m1_null_pointer
     || ghl_m1_compute_wavespeeds(
              metric, (ghl_m1_direction_t)3, &speed_minus, &speed_plus)
              != ghl_error_m1_invalid_state) {
    fail_test("combined wave-speed boundary failed");
  }

  const double grad[3] = { 0.2, -0.1, 0.05 };
  const double V[3] = { 0.0, 0.0, 0.0 };
  const double valid_W = 1.0;
  double corrected = 61.0;
  double a_face = 62.0;
  if(ghl_m1_compute_diffusion_flux(
           NULL, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W, V,
           1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_null_pointer
     || corrected != 61.0 || a_face != 62.0) {
    fail_test("diffusion NULL arguments were not rejected");
  }
  corrected = 63.0;
  a_face = 64.0;
  if(ghl_m1_compute_diffusion_flux(
           params, NULL, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W, V,
           1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_null_pointer
     || corrected != 63.0 || a_face != 64.0) {
    fail_test("diffusion NULL metric changed outputs");
  }
  corrected = 65.0;
  a_face = 66.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, NULL, valid_W,
           V, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_null_pointer
     || corrected != 65.0 || a_face != 66.0) {
    fail_test("diffusion NULL gradient changed outputs");
  }
  corrected = 67.0;
  a_face = 68.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W,
           NULL, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_null_pointer
     || corrected != 67.0 || a_face != 68.0) {
    fail_test("diffusion NULL velocity changed outputs");
  }
  a_face = 69.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W,
           V, 1.0, 0.2, 1.0, NULL, &a_face)
           != ghl_error_m1_null_pointer
     || a_face != 69.0) {
    fail_test("diffusion NULL corrected output changed optional output");
  }

  corrected = 61.0;
  a_face = 62.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, (ghl_m1_direction_t)3, 0.4, 0.1, 1.0, true, 1.0, true, grad,
           valid_W, V, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 61.0 || a_face != 62.0) {
    fail_test("invalid diffusion direction changed outputs");
  }
  invalid_metric = *metric;
  invalid_metric.gammaDD[0][0] = -1.0;
  corrected = 63.0;
  a_face = 64.0;
  if(ghl_m1_compute_diffusion_flux(
           params, &invalid_metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad,
           valid_W, V, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_metric
     || corrected != 63.0 || a_face != 64.0) {
    fail_test("invalid diffusion metric changed outputs");
  }
  corrected = 81.0;
  a_face = 82.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, NAN, 0.1, 1.0, true, 1.0, true, grad, valid_W,
           V, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 81.0 || a_face != 82.0) {
    fail_test("nonfinite diffusion HLL/star inputs were accepted");
  }
  corrected = 83.0;
  a_face = 84.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, NAN, 1.0, true, 1.0, true, grad, valid_W,
           V, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 83.0 || a_face != 84.0) {
    fail_test("nonfinite diffusion star energy changed outputs");
  }
  const double invalid_chi[] = { NAN, -1.0 };
  const double invalid_D[] = { NAN, -1.0 };
  const double invalid_delta[] = { NAN, 0.0 };
  for(int case_index = 0; case_index < 2; ++case_index) {
    corrected = 85.0;
    a_face = 86.0;
    if(ghl_m1_compute_diffusion_flux(
             params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W,
             V, invalid_chi[case_index], 0.2, 1.0, &corrected, &a_face)
             != ghl_error_m1_invalid_state
       || corrected != 85.0 || a_face != 86.0) {
      fail_test("invalid diffusion opacity changed outputs");
    }
    corrected = 87.0;
    a_face = 88.0;
    if(ghl_m1_compute_diffusion_flux(
             params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W,
             V, 1.0, invalid_D[case_index], 1.0, &corrected, &a_face)
             != ghl_error_m1_invalid_state
       || corrected != 87.0 || a_face != 88.0) {
      fail_test("invalid diffusion coefficient changed outputs");
    }
    corrected = 89.0;
    a_face = 90.0;
    if(ghl_m1_compute_diffusion_flux(
             params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W,
             V, 1.0, 0.2, invalid_delta[case_index], &corrected, &a_face)
             != ghl_error_m1_invalid_state
       || corrected != 89.0 || a_face != 90.0) {
      fail_test("invalid diffusion cell width changed outputs");
    }
  }

  corrected = 71.0;
  a_face = 72.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W,
           V, DBL_MAX, 0.2, DBL_MAX, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 71.0 || a_face != 72.0) {
    fail_test("overflowing diffusion optical width was accepted");
  }

  const bool validity_cases[3] = { false, true, true };
  const double J_cases[3] = { 1.0, NAN, 0.0 };
  for(int case_index = 0; case_index < 3; ++case_index) {
    if(ghl_m1_compute_diffusion_flux(
             params, metric, ghl_m1_dirn0, 0.4, 0.1, J_cases[case_index],
             validity_cases[case_index], 1.0, true, grad, valid_W, V, 1.0, 0.2, 1.0,
             &corrected, &a_face)
             != ghl_success
       || corrected != 0.4 || a_face != 1.0
       || ghl_m1_compute_diffusion_flux(
                params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, J_cases[case_index],
                validity_cases[case_index], grad, valid_W, V, 1.0, 0.2, 1.0, &corrected,
                &a_face)
                != ghl_success
       || corrected != 0.4 || a_face != 1.0) {
      fail_test("invalid Jthick validity did not take the no-op path");
    }
  }

  const double bad_grad[3] = { NAN, 0.0, 0.0 };
  corrected = 73.0;
  a_face = 74.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, bad_grad,
           valid_W, V, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 73.0 || a_face != 74.0) {
    fail_test("nonfinite diffusion gradient was accepted");
  }
  const double huge_grad[3] = { 2.0, 0.0, 0.0 };
  corrected = 65.0;
  a_face = 66.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, huge_grad,
           valid_W, V, 1.0, DBL_MAX, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 65.0 || a_face != 66.0) {
    fail_test("overflowing diffusion asymptotic flux changed outputs");
  }
  ghl_metric_quantities huge_volume_metric = *metric;
  huge_volume_metric.gammaDD[0][0] = 4.0;
  huge_volume_metric.gammaUU[0][0] = 0.25;
  huge_volume_metric.detgamma = 4.0;
  huge_volume_metric.sqrt_detgamma = 2.0;
  huge_volume_metric.betaU[0] = -2.0;
  corrected = 75.0;
  a_face = 76.0;
  if(ghl_m1_compute_diffusion_flux(
           params, &huge_volume_metric, ghl_m1_dirn0, 0.4, DBL_MAX, 1.0, true, 1.0, true,
           grad, valid_W, V, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 75.0 || a_face != 76.0) {
    fail_test("overflowing diffusion volume flux was accepted");
  }

  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, NAN, V,
           1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_success
     || corrected != 0.4 || a_face != 1.0
     || ghl_m1_compute_diffusion_flux(
              params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, 0.5, V,
              1.0, 0.2, 1.0, &corrected, &a_face)
              != ghl_success
     || corrected != 0.4 || a_face != 1.0) {
    fail_test("invalid face Lorentz factor did not take no-op path");
  }
  const double superluminal_V[3] = { DBL_MAX, 0.0, 0.0 };
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, 1.0,
           superluminal_V, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_success
     || corrected != 0.4 || a_face != 1.0) {
    fail_test("nonfinite face velocity norm did not take no-op path");
  }
  const double bad_velocity_component[3] = { NAN, 0.0, 0.0 };
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, 1.0,
           bad_velocity_component, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_success
     || corrected != 0.4 || a_face != 1.0) {
    fail_test("nonfinite face velocity component did not take no-op path");
  }
  const double moving_velocity[3] = { 0.2, 0.1, 0.0 };
  const double moving_W = 1.0 / sqrt(1.0 - 0.2 * 0.2 - 0.1 * 0.1);
  double corrected_without_a = NAN;
  double corrected_with_a = NAN;
  double observed_a = NAN;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, moving_W,
           moving_velocity, 1.0, 0.2, 1.0, &corrected_without_a, NULL)
           != ghl_success
     || ghl_m1_compute_diffusion_flux(
              params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad,
              moving_W, moving_velocity, 1.0, 0.2, 1.0, &corrected_with_a, &observed_a)
              != ghl_success
     || !isfinite(corrected_without_a) || !isfinite(corrected_with_a)
     || !isfinite(observed_a) || corrected_without_a != corrected_with_a) {
    fail_test("optional diffusion a_face output changed the flux");
  }
  corrected = 77.0;
  a_face = 78.0;
  if(ghl_m1_compute_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, DBL_MAX, true, DBL_MAX, true, grad,
           moving_W, moving_velocity, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 77.0 || a_face != 78.0) {
    fail_test("overflowing positive Jthick average was accepted");
  }

  ghl_m1_neutrino_rates rates = { 0 };
  rates.kappa_tr = 1.0;
  corrected = 79.0;
  a_face = 80.0;
  if(ghl_m1_compute_neutrino_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W,
           V, NULL, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_null_pointer
     || corrected != 79.0 || a_face != 80.0) {
    fail_test("neutrino diffusion NULL rates were accepted");
  }
  if(ghl_m1_compute_neutrino_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true, grad, valid_W,
           V, &rates, 0.2, 1.0, &corrected, &a_face)
     != ghl_success) {
    fail_test("valid neutrino diffusion wrapper failed");
  }
  corrected = 67.0;
  a_face = 68.0;
  if(ghl_m1_compute_neutrino_diffusion_flux(
           params, metric, (ghl_m1_direction_t)3, 0.4, 0.1, 1.0, true, 1.0, true, grad,
           valid_W, V, &rates, 0.2, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 67.0 || a_face != 68.0) {
    fail_test("invalid neutrino diffusion direction changed outputs");
  }
  const double wrapper_huge_grad[3] = { 2.0, 0.0, 0.0 };
  corrected = 69.0;
  a_face = 70.0;
  if(ghl_m1_compute_neutrino_diffusion_flux(
           params, metric, ghl_m1_dirn0, 0.4, 0.1, 1.0, true, 1.0, true,
           wrapper_huge_grad, valid_W, V, &rates, DBL_MAX, 1.0, &corrected, &a_face)
           != ghl_error_m1_invalid_state
     || corrected != 69.0 || a_face != 70.0) {
    fail_test("overflowing neutrino diffusion flux changed outputs");
  }
}

/* E, F_i, lapse, beta^i, gammaDD[9] row-major, and coordinate velocity. */
enum { jthick_fixture_input_count = 20 };

static int setup_jthick_fixture_state(
      const double input[jthick_fixture_input_count],
      ghl_metric_quantities *restrict metric,
      ghl_primitive_quantities *restrict prims,
      ghl_m1_rad_state *restrict rad_state) {
  ghl_initialize_metric(
        input[4], input[5], input[6], input[7], input[8], input[9], input[10], input[12],
        input[13], input[16], metric);
  *prims = (ghl_primitive_quantities){ 0 };
  for(int i = 0; i < 3; ++i) {
    prims->vU[i] = input[17 + i];
  }
  double V_con[3];
  for(int i = 0; i < 3; ++i) {
    V_con[i] = (prims->vU[i] + metric->betaU[i]) / metric->lapse;
  }
  double v2 = 0.0;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      v2 += metric->gammaDD[i][j] * V_con[i] * V_con[j];
    }
  }
  if(!isfinite(v2) || v2 < 0.0 || v2 >= 1.0) {
    return 0;
  }
  prims->u0 = 1.0 / (metric->lapse * sqrt(1.0 - v2));
  *rad_state
        = (ghl_m1_rad_state){ .E = input[0], .F = { input[1], input[2], input[3] } };
  return 1;
}

static void check_jthick_thcm1_fixtures(
      const ghl_m1_parameters *restrict params,
      const char *restrict fixture_path) {
  static const char *const expected_case_ids[] = {
    "rngpkt-v2-rd-radiation-seed-v1-0000-a00", "rngpkt-v2-rd-radiation-seed-v1-0000-a01",
    "rngpkt-v2-rd-radiation-seed-v1-0000-a02", "rngpkt-v2-rd-radiation-seed-v1-0001-a00",
    "rngpkt-v2-rd-radiation-seed-v1-0001-a01", "rngpkt-v2-rd-radiation-seed-v1-0001-a02"
  };
  const size_t expected_case_count
        = sizeof(expected_case_ids) / sizeof(expected_case_ids[0]);
  m1_thcm1_fixture_collection fixtures = { 0 };
  char error[M1_THCM1_FIXTURE_TEXT_MAX] = { 0 };
  if(!m1_thcm1_fixture_load(
           fixture_path, "jthick", jthick_fixture_input_count, 1, &fixtures, error,
           sizeof(error))) {
    ghl_error("%s: %s\n", fixture_path, error);
  }
  if(strcmp(fixtures.policy, "pointwise_a1_a2_v1") != 0) {
    ghl_error("unrecognized Jthick fixture comparison policy\n");
  }
  if(fixtures.record_count != expected_case_count) {
    ghl_error("Jthick fixture does not contain the six seeded radiation pairs\n");
  }

  for(size_t i = 0; i < fixtures.record_count; ++i) {
    const m1_thcm1_fixture_record *record = &fixtures.records[i];
    if(strcmp(record->case_id, expected_case_ids[i]) != 0
       || !m1_thcm1_fixture_pair_id_matches_case_suffix(record, "-pair")
       || strcmp(record->origin, "seeded") != 0
       || strcmp(record->family, "radiation") != 0
       || strcmp(record->perturbation, "radiation.E") != 0
       || record->sensitivity_start != 0 || record->sensitivity_count != 1) {
      ghl_error("Jthick fixture pair metadata differs at %zu\n", i);
    }

    const double normalization
          = fmax(fabs(record->baseline_input[0]), fabs(record->perturbed_input[0]));
    if(normalization != record->normalization[0]) {
      ghl_error("Jthick fixture normalization differs from its input at %zu\n", i);
    }

    const double *inputs[2] = { record->baseline_input, record->perturbed_input };
    const char *const labels[2] = { "baseline", "perturbed" };
    double actual[2] = { NAN, NAN };
    for(int state_index = 0; state_index < 2; ++state_index) {
      ghl_metric_quantities metric;
      ghl_primitive_quantities prims;
      ghl_m1_rad_state rad_state;
      if(!setup_jthick_fixture_state(inputs[state_index], &metric, &prims, &rad_state)) {
        ghl_error(
              "Jthick fixture %s state is not an admissible primitive state: %s\n",
              labels[state_index], record->case_id);
      }
      bool valid = false;
      const ghl_error_codes_t jthick_error = ghl_m1_compute_Jthick(
            params, &metric, &prims, &rad_state, &actual[state_index], &valid);
      if(jthick_error != ghl_success || !valid || !isfinite(actual[state_index])) {
        ghl_error(
              "Jthick fixture %s evaluation failed: %s error=%d actual=%.17g valid=%d\n",
              labels[state_index], record->case_id, jthick_error, actual[state_index],
              valid);
      }
    }
    const double computed_normalization[1] = { normalization };
    m1_thcm1_fixture_comparison_report report = { 0 };
    char compare_error[M1_THCM1_FIXTURE_TEXT_MAX] = { 0 };
    if(!m1_thcm1_fixture_compare_paired(
             record, computed_normalization, actual, actual + 1, &report, compare_error,
             sizeof(compare_error))) {
      ghl_error(
            "Jthick fixture paired comparison failed: %s gate=%s component=%zu\n",
            compare_error, report.gate != NULL ? report.gate : "unknown",
            report.component);
    }
  }
  ghl_info(
        "unit_test_m1_diffusion_flux: %zu baseline/perturbed THC_M1 Jthick pairs "
        "passed baseline, perturbed, and response comparisons\n",
        fixtures.record_count);
  m1_thcm1_fixture_free(&fixtures);
}

int main(int argc, char **argv) {
  const char *jthick_fixture_path = "Unit_Tests/data/jthick_thcm1.fixture";
  if(argc == 3 && strcmp(argv[1], "--fixture") == 0) {
    jthick_fixture_path = argv[2];
  }
  else if(argc != 1) {
    ghl_error("Usage: %s [--fixture PATH]\n", argv[0]);
  }
  const ghl_m1_neutrino_state *volatile no_state = NULL;
  double rejected_J = -1.0;
  bool rejected_valid = true;
  if(ghl_m1_compute_neutrino_Jthick(
           NULL, NULL, NULL, no_state, &rejected_J, &rejected_valid)
           != ghl_error_m1_null_pointer
     || rejected_J != -1.0 || !rejected_valid) {
    fail_test("NULL neutrino state did not preserve Jthick outputs");
  }
  ghl_m1_parameters params = { 0 };
  if(ghl_m1_initialize(1.0e-10, 0.1, 0.25, 1.0e-6, 1.0e-12, 20, 1.0e-10, &params)
     != ghl_success) {
    fail_test("M1 initialization failed");
  }

  ghl_metric_quantities flat_metric;
  m1_setup_flat_metric(&flat_metric);
  ghl_primitive_quantities static_prims;
  setup_primitives(&static_prims, 0.0, 0.0, 0.0);
  const ghl_m1_rad_state static_state = { .E = 1.2, .F = { 0.2, 0.1, 0.0 } };
  double Jthick = NAN;
  bool valid = false;
  if(ghl_m1_compute_Jthick(
           &params, &flat_metric, &static_prims, &static_state, &Jthick, &valid)
           != ghl_success
     || !valid) {
    fail_test("static Jthick computation failed");
  }
  check_close(Jthick, static_state.E, "static Jthick should equal E");

  ghl_primitive_quantities moving_prims;
  setup_primitives(&moving_prims, 0.2, 0.1, 0.0);
  double moving_J = NAN;
  bool moving_valid = false;
  if(ghl_m1_compute_neutrino_Jthick(
           &params, &flat_metric, &moving_prims,
           (const ghl_m1_neutrino_state *)&(const ghl_m1_neutrino_state){
                 .N = 0.4,
                 .E = static_state.E,
                 .F = { static_state.F[0], static_state.F[1], static_state.F[2] } },
           &moving_J, &moving_valid)
           != ghl_success
     || !moving_valid) {
    fail_test("moving Jthick computation failed");
  }
  const double W2 = 1.0 / (1.0 - 0.2 * 0.2 - 0.1 * 0.1);
  const double FdotV = static_state.F[0] * 0.2 + static_state.F[1] * 0.1;
  const double expected_moving_J
        = 3.0 / (2.0 * W2 + 1.0)
          * ((2.0 * W2 - 1.0) * static_state.E - 2.0 * W2 * FdotV);
  check_close(moving_J, expected_moving_J, "moving Jthick formula failed");

  double D_face = NAN;
  if(ghl_m1_compute_harmonic_diffusion_coefficient(1.0, 2.0, &D_face) != ghl_success) {
    fail_test("harmonic diffusion coefficient failed");
  }
  check_close(D_face, 2.0 / 9.0, "harmonic diffusion coefficient mismatch");
  if(ghl_m1_compute_harmonic_diffusion_coefficient(DBL_MAX, DBL_MAX, &D_face)
           != ghl_success
     || D_face != (1.0 / 3.0) / DBL_MAX) {
    fail_test("representable harmonic coefficient was rejected at DBL_MAX opacity");
  }
  if(ghl_m1_compute_harmonic_diffusion_coefficient(1.0e-200, 1.0e-200, &D_face)
           != ghl_success
     || !m1_nearly_equal(D_face, 1.0 / (3.0e-200), 3.0e-12, 0.0)) {
    fail_test("harmonic coefficient overflowed before finite division");
  }
  D_face = 17.0;
  if(ghl_m1_compute_harmonic_diffusion_coefficient(0.0, 2.0, &D_face)
           != ghl_error_m1_invalid_state
     || D_face != 17.0) {
    fail_test("invalid harmonic coefficient was not transactional");
  }

  const double grad[3] = { 0.3, 0.0, 0.0 };
  const double V_static[3] = { 0.0, 0.0, 0.0 };
  double corrected = NAN;
  double a_face = NAN;
  if(ghl_m1_compute_diffusion_flux(
           &params, &flat_metric, ghl_m1_dirn0, 0.7, 0.06, 0.8, true, 0.92, true, grad,
           1.0, V_static, 0.4, 1.0 / 6.0, 0.5, &corrected, &a_face)
           != ghl_success
     || corrected != 0.7 || a_face != 1.0) {
    fail_test("below-threshold diffusion gate was not an exact no-op");
  }

  /* At the threshold the active branch is selected, including E_floor in the
   * shift term. This is a pinned independent oracle for the entire blend. */
  if(ghl_m1_compute_diffusion_flux(
           &params, &flat_metric, ghl_m1_dirn0, 0.7, 0.06, 0.8, true, 0.92, true, grad,
           1.0, V_static, 0.5, 1.0 / 6.0, 0.5, &corrected, &a_face)
     != ghl_success) {
    fail_test("threshold diffusion branch failed");
  }
  double expected_a = NAN;
  const double expected = expected_diffusion_flux(
        &params, &flat_metric, 0, 0.7, 0.06, 0.8, 0.92, grad, 1.0, V_static, 0.5,
        1.0 / 6.0, 0.5, &expected_a);
  check_close(a_face, expected_a, "threshold blend factor mismatch");
  check_close(corrected, expected, "threshold corrected flux mismatch");

  ghl_metric_quantities curved_metric;
  m1_setup_metric_anchor_B(&curved_metric);
  for(int direction = ghl_m1_dirn0; direction <= ghl_m1_dirn2; ++direction) {
    double delta_l = NAN;
    if(ghl_m1_compute_face_normal_delta_l(
             &curved_metric, (ghl_m1_direction_t)direction, 0.5, &delta_l)
       != ghl_success) {
      fail_test("face-normal proper thickness failed");
    }
    check_close(
          delta_l, 0.5 / sqrt(curved_metric.gammaUU[direction][direction]),
          "face-normal proper thickness mismatch");
  }
  double unchanged_delta_l = 19.0;
  if(ghl_m1_compute_face_normal_delta_l(
           &curved_metric, (ghl_m1_direction_t)3, 0.5, &unchanged_delta_l)
           != ghl_error_m1_invalid_state
     || unchanged_delta_l != 19.0) {
    fail_test("invalid face direction was not transactional");
  }

  const double V_moving[3] = { 0.2, 0.1, 0.0 };
  const double W_moving
        = 1.0 / sqrt(1.0 - ghl_compute_vec2_from_vec3D(curved_metric.gammaDD, V_moving));
  const double moving_grad[3] = { 0.3, -0.12, 0.07 };
  if(ghl_m1_compute_diffusion_flux(
           &params, &curved_metric, ghl_m1_dirn0, 0.22, 0.8, 0.8, true, 0.92, true,
           moving_grad, W_moving, V_moving, 1.4, 0.20, 0.5, &corrected, &a_face)
     != ghl_success) {
    fail_test("moving active diffusion branch failed");
  }
  const double moving_expected = expected_diffusion_flux(
        &params, &curved_metric, 0, 0.22, 0.8, 0.8, 0.92, moving_grad, W_moving,
        V_moving, 1.4, 0.20, 0.5, &expected_a);
  check_close(a_face, expected_a, "moving blend factor mismatch");
  check_close(corrected, moving_expected, "moving corrected flux mismatch");

  /* An invalid face velocity is a deliberate no-op, not a published flux. */
  const double invalid_V[3] = { 1.0, 0.0, 0.0 };
  if(ghl_m1_compute_diffusion_flux(
           &params, &flat_metric, ghl_m1_dirn0, 0.42, 1.0, 0.8, true, 0.9, true, grad,
           1.0, invalid_V, 1.0, 0.2, 1.0, &corrected, &a_face)
           != ghl_success
     || corrected != 0.42 || a_face != 1.0) {
    fail_test("invalid face velocity did not take the documented no-op path");
  }

  ghl_m1_neutrino_rates rates = { 0 };
  rates.kappa_tr = 0.5;
  if(ghl_m1_compute_neutrino_diffusion_flux(
           &params, &flat_metric, ghl_m1_dirn0, 0.7, 0.06, 0.8, true, 0.92, true, grad,
           1.0, V_static, &rates, 1.0 / 6.0, 0.5, &corrected, &a_face)
     != ghl_success) {
    fail_test("neutrino diffusion wrapper failed");
  }
  check_close(corrected, expected, "neutrino diffusion wrapper mismatch");

  double rejected = 37.0;
  if(ghl_m1_compute_diffusion_flux(
           &params, &flat_metric, ghl_m1_dirn0, NAN, 0.1, 0.8, true, 0.9, true, grad,
           1.0, V_static, 0.5, 0.2, 0.5, &rejected, NULL)
           != ghl_error_m1_invalid_state
     || rejected != 37.0) {
    fail_test("invalid diffusion input was not transactional");
  }

  check_diffusion_transport_boundaries(
        &params, &flat_metric, &static_prims, &static_state);
  ghl_m1_parameters fixture_params = { 0 };
  if(ghl_m1_initialize(
           1.0e-10, 1.0e-12, 1.0e-8, 1.0e-6, 1.0e-12, 20, 1.0e-10, &fixture_params)
     != ghl_success) {
    fail_test("Jthick fixture parameters failed to initialize");
  }
  check_jthick_thcm1_fixtures(&fixture_params, jthick_fixture_path);

  ghl_info("unit_test_m1_diffusion_flux: Jthick/harmonic/diffusion cases passed\n");
  return 0;
}
