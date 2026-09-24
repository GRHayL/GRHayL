#include "ghl_m1.h"
#include "ghl_m1_utils.h"

static bool ghl_m1_same_nonzero_sign(const double left, const double right) {
  return (left > 0.0 && right > 0.0) || (left < 0.0 && right < 0.0);
}

static bool ghl_m1_opposite_nonzero_sign(const double left, const double right) {
  return (left > 0.0 && right < 0.0) || (left < 0.0 && right > 0.0);
}

/* Both operands are finite. Keep the original rounding, including subnormal
 * sums, unless the addition itself would overflow. */
static double ghl_m1_average_finite(const double left, const double right) {
  const double sum = left + right;
  return isfinite(sum) ? 0.5 * sum : 0.5 * left + 0.5 * right;
}

static double ghl_m1_half_scaled_rusanov_candidate(
      const double physical_flux_L, const double physical_flux_R,
      const double speed, const double jump, const double jump_low,
      const bool jump_is_half) {
  const double quarter_L = 0.25 * physical_flux_L;
  const double quarter_R = 0.25 * physical_flux_R;
  const double half_average = quarter_L + quarter_R;
  const double average_part = half_average - quarter_L;
  const double average_low = (quarter_L - (half_average - average_part))
                             + (quarter_R - average_part);

  const double coefficient = jump_is_half ? -0.5 * speed : -0.25 * speed;
  const double half_product = coefficient * jump;
  if(!isfinite(half_product)) {
    return half_product;
  }
  const double product_low = fma(coefficient, jump, -half_product);
  const double half_sum = half_product + half_average;
  const double sum_part = half_sum - half_product;
  const double sum_low = (half_product - (half_sum - sum_part))
                         + (half_average - sum_part);
  const double jump_low_product = coefficient * jump_low;
  return 2.0 * (half_sum + (((sum_low + average_low) + product_low)
                            + jump_low_product));
}

static double ghl_m1_rusanov_candidate_finite(
      const double state_L, const double state_R,
      const double physical_flux_L, const double physical_flux_R,
      const double speed) {
  const double average = ghl_m1_average_finite(physical_flux_L, physical_flux_R);
  const double jump = state_R - state_L;
  if(isfinite(jump)) {
    const double candidate = average - 0.5 * speed * jump;
    if(isfinite(candidate)) {
      return candidate;
    }
    const double fused = fma(-0.5 * speed, jump, average);
    return isfinite(fused) ? fused
                           : ghl_m1_half_scaled_rusanov_candidate(
                                   physical_flux_L, physical_flux_R, speed,
                                   jump, 0.0, false);
  }
  const double half_R = 0.5 * state_R;
  const double minus_half_L = -0.5 * state_L;
  const double half_jump = half_R + minus_half_L;
  const double rounded_part = half_jump - half_R;
  const double low_part = (half_R - (half_jump - rounded_part))
                          + (minus_half_L - rounded_part);
  const double candidate = fma(-speed, half_jump, average);
  const double corrected = fma(-speed, low_part, candidate);
  if(isfinite(corrected)) {
    return corrected;
  }
  return ghl_m1_half_scaled_rusanov_candidate(
        physical_flux_L, physical_flux_R, speed, half_jump, low_part, true);
}

typedef struct {
  double mantissa;
  int exponent;
} ghl_m1_half_slope;

/* Represent every stencil slope at the same half scale. frexp retains the
 * sign of a small slope even when a different slope exceeds double range. */
static ghl_m1_half_slope ghl_m1_scaled_half_slope(
      const double left, const double right) {
  ghl_m1_half_slope slope;
  const double difference = right - left;
  if(isfinite(difference)) {
    slope.mantissa = frexp(difference, &slope.exponent);
    --slope.exponent;
  }
  else {
    const double half_difference = 0.5 * right - 0.5 * left;
    slope.mantissa = frexp(half_difference, &slope.exponent);
  }
  return slope;
}

static double ghl_m1_half_slope_ratio(
      const ghl_m1_half_slope numerator, const ghl_m1_half_slope denominator) {
  return scalbn(numerator.mantissa / denominator.mantissa,
                numerator.exponent - denominator.exponent);
}

static ghl_error_codes_t ghl_m1_compute_overflowing_stencil_limiter(
      const ghl_m1_parameters *restrict m1_params,
      const double state_0, const double state_1,
      const double state_2, const double state_3,
      double *restrict phi, bool *restrict sawtooth) {
  if(!isfinite(m1_params->minmod_theta) || m1_params->minmod_theta < 0.0
     || m1_params->minmod_theta > 2.0) {
    return ghl_error_m1_invalid_state;
  }

  const ghl_m1_half_slope dum = ghl_m1_scaled_half_slope(state_0, state_1);
  const ghl_m1_half_slope duc = ghl_m1_scaled_half_slope(state_1, state_2);
  const ghl_m1_half_slope dup = ghl_m1_scaled_half_slope(state_2, state_3);
  double candidate_phi = 0.0;
  bool candidate_sawtooth = false;
  if(ghl_m1_same_nonzero_sign(dup.mantissa, duc.mantissa)
     && ghl_m1_same_nonzero_sign(dum.mantissa, duc.mantissa)) {
    if(m1_params->minmod_theta != 0.0) {
      const double ratio_left = ghl_m1_half_slope_ratio(dum, duc);
      const double ratio_right = ghl_m1_half_slope_ratio(dup, duc);
      candidate_phi = fmin(1.0, fmin(m1_params->minmod_theta * ratio_left,
                                     m1_params->minmod_theta * ratio_right));
    }
  }
  else if(ghl_m1_opposite_nonzero_sign(dup.mantissa, duc.mantissa)
          && ghl_m1_opposite_nonzero_sign(dum.mantissa, duc.mantissa)) {
    candidate_sawtooth = true;
  }
  *phi = candidate_phi;
  *sawtooth = candidate_sawtooth;
  return ghl_success;
}

static ghl_error_codes_t
ghl_m1_validate_four_point_transport_policy(const bool diffusion_correction_enabled) {
  /* The separate diffusion correction is not part of this transport
   * operator. Rejecting it here prevents a caller from accidentally
   * obtaining a different face flux. */
  if(diffusion_correction_enabled) {
    return ghl_error_m1_incompatible_transport_policy;
  }
  return ghl_success;
}

/* The public Flux_Source Rusanov helper documents undensitized operands.
 * Volume-weighted transport uses the same componentwise arithmetic, but its
 * state and physical-flux arrays are intentionally measured in the
 * volume-weighted units used by
 * ghl_m1_compute_neutrino_four_point_volume_weighted_transport_flux. Keep
 * that use local so the public helper contract is not broadened implicitly. */
static ghl_error_codes_t ghl_m1_compute_neutrino_four_point_rusanov_flux(
      const double state_L[ghl_m1_neutrino_transport_component_count],
      const double state_R[ghl_m1_neutrino_transport_component_count],
      const double physical_flux_L[ghl_m1_neutrino_transport_component_count],
      const double physical_flux_R[ghl_m1_neutrino_transport_component_count],
      const double speed,
      double flux[ghl_m1_neutrino_transport_component_count]) {
  if(state_L == NULL || state_R == NULL || physical_flux_L == NULL
     || physical_flux_R == NULL || flux == NULL || !isfinite(speed) || speed < 0.0) {
    return ghl_error_m1_invalid_state;
  }

  double candidate_flux[ghl_m1_neutrino_transport_component_count];
  for(int component = 0; component < ghl_m1_neutrino_transport_component_count;
      ++component) {
    if(!isfinite(state_L[component]) || !isfinite(state_R[component])
       || !isfinite(physical_flux_L[component])
       || !isfinite(physical_flux_R[component])) {
      return ghl_error_m1_invalid_state;
    }

    candidate_flux[component] = ghl_m1_rusanov_candidate_finite(
          state_L[component], state_R[component], physical_flux_L[component],
          physical_flux_R[component], speed);
    if(!isfinite(candidate_flux[component])) {
      return ghl_error_m1_invalid_state;
    }
  }

  for(int component = 0; component < ghl_m1_neutrino_transport_component_count;
      ++component) {
    flux[component] = candidate_flux[component];
  }
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_four_point_flux_limiter(
      const ghl_m1_parameters *restrict m1_params,
      const double dum,
      const double duc,
      const double dup,
      double *restrict phi,
      bool *restrict sawtooth) {
  if(m1_params == NULL || phi == NULL || sawtooth == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(dum) || !isfinite(duc) || !isfinite(dup)
     || !isfinite(m1_params->minmod_theta) || m1_params->minmod_theta < 0.0
     || m1_params->minmod_theta > 2.0) {
    return ghl_error_m1_invalid_state;
  }

  double candidate_phi = 0.0;
  bool candidate_sawtooth = false;
  if(ghl_m1_same_nonzero_sign(dup, duc) && ghl_m1_same_nonzero_sign(dum, duc)) {
    /* theta == 0 selects the fully dissipative low-order flux.  Forming the
     * ratios first lets an overflowing ratio turn 0 * inf into NaN, and fmin
     * then returns its other argument -- silently publishing the undissipated
     * phi = 1.  Decide on theta before any division. */
    if(m1_params->minmod_theta == 0.0) {
      candidate_phi = 0.0;
    }
    else {
      const double ratio_left = dum / duc;
      const double ratio_right = dup / duc;
      candidate_phi = fmin(
            1.0, fmin(m1_params->minmod_theta * ratio_left,
                      m1_params->minmod_theta * ratio_right));
    }
  }
  else if(
        ghl_m1_opposite_nonzero_sign(dup, duc)
        && ghl_m1_opposite_nonzero_sign(dum, duc)) {
    candidate_sawtooth = true;
  }

  if(!isfinite(candidate_phi)) {
    return ghl_error_m1_invalid_state;
  }
  *phi = candidate_phi;
  *sawtooth = candidate_sawtooth;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_four_point_opacity_suppression(
      const ghl_m1_parameters *restrict m1_params,
      const double kappa_face,
      const double delta,
      double *restrict A) {
  if(m1_params == NULL || A == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(kappa_face) || kappa_face < 0.0 || !isfinite(delta) || delta <= 0.0
     || !isfinite(m1_params->mindiss) || m1_params->mindiss < 0.0
     || m1_params->mindiss > 1.0) {
    return ghl_error_m1_invalid_state;
  }

  double candidate_A = 1.0;
  const double optical_width = delta * kappa_face;
  if(optical_width > 1.0) {
    candidate_A = fmin(1.0, 1.0 / optical_width);
    candidate_A = fmax(candidate_A, m1_params->mindiss);
  }
  if(!isfinite(candidate_A)) {
    return ghl_error_m1_invalid_state;
  }
  *A = candidate_A;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_four_point_blended_flux(
      const double flux_high,
      const double flux_low,
      const double phi,
      const bool sawtooth,
      const double A,
      double *restrict flux_num) {
  if(flux_num == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(flux_high) || !isfinite(flux_low) || !isfinite(phi) || phi < 0.0
     || phi > 1.0 || !isfinite(A) || A < 0.0 || A > 1.0) {
    return ghl_error_m1_invalid_state;
  }

  const double dissipation = sawtooth ? 1.0 : A;
  const double weight = dissipation * (1.0 - phi);
  double candidate_flux;
  if(weight == 0.0) {
    candidate_flux = flux_high;
  }
  else if(weight == 1.0) {
    candidate_flux = flux_low;
  }
  else {
    const double difference = flux_high - flux_low;
    /* The difference may overflow while the convex blend remains finite. */
    candidate_flux = isfinite(difference)
                           ? flux_high - weight * difference
                           : (1.0 - weight) * flux_high + weight * flux_low;
  }
  if(!isfinite(candidate_flux)) {
    return ghl_error_m1_invalid_state;
  }
  *flux_num = candidate_flux;
  return ghl_success;
}

static ghl_error_codes_t ghl_m1_compute_neutrino_four_point_transport_core(
      const ghl_m1_parameters *restrict m1_params,
      const double state_stencil[4][ghl_m1_neutrino_transport_component_count],
      const double physical_flux_L[ghl_m1_neutrino_transport_component_count],
      const double physical_flux_R[ghl_m1_neutrino_transport_component_count],
      const double speed_L,
      const double speed_R,
      const double kappa_face,
      const double delta_x,
      const bool diffusion_correction_enabled,
      double flux_tilde[ghl_m1_neutrino_transport_component_count],
      ghl_m1_four_point_transport_diagnostics *restrict candidate_diagnostics) {
  /* Keep all candidates local until every component has passed validation so
   * that a failed volume-weighted operation cannot publish a partial result. */
  if(m1_params == NULL || state_stencil == NULL || physical_flux_L == NULL
     || physical_flux_R == NULL || flux_tilde == NULL || candidate_diagnostics == NULL) {
    return ghl_error_m1_null_pointer;
  }

  const ghl_error_codes_t policy_error
        = ghl_m1_validate_four_point_transport_policy(diffusion_correction_enabled);
  if(policy_error != ghl_success) {
    return policy_error;
  }

  if(!isfinite(speed_L) || speed_L < 0.0 || !isfinite(speed_R) || speed_R < 0.0
     || !isfinite(kappa_face) || kappa_face < 0.0 || !isfinite(delta_x)
     || delta_x <= 0.0) {
    return ghl_error_m1_invalid_state;
  }

  const double face_speed = fmax(speed_L, speed_R);
  if(!isfinite(face_speed) || face_speed < 0.0) {
    return ghl_error_m1_invalid_state;
  }

  /* Validate all volume-weighted operands before invoking the private
   * componentwise Rusanov arithmetic. */
  double state_L[ghl_m1_neutrino_transport_component_count];
  double state_R[ghl_m1_neutrino_transport_component_count];
  double flux_L[ghl_m1_neutrino_transport_component_count];
  double flux_R[ghl_m1_neutrino_transport_component_count];
  for(int component = 0; component < ghl_m1_neutrino_transport_component_count;
      ++component) {
    state_L[component] = state_stencil[1][component];
    state_R[component] = state_stencil[2][component];
    flux_L[component] = physical_flux_L[component];
    flux_R[component] = physical_flux_R[component];
    if(!isfinite(state_stencil[0][component]) || !isfinite(state_stencil[1][component])
       || !isfinite(state_stencil[2][component])
       || !isfinite(state_stencil[3][component]) || !isfinite(flux_L[component])
       || !isfinite(flux_R[component])) {
      return ghl_error_m1_invalid_state;
    }
  }

  double flux_low[ghl_m1_neutrino_transport_component_count];
  const ghl_error_codes_t rusanov_error
        = ghl_m1_compute_neutrino_four_point_rusanov_flux(
              state_L, state_R, flux_L, flux_R, face_speed, flux_low);
  if(rusanov_error != ghl_success) {
    return ghl_error_m1_invalid_state;
  }

  double candidate_flux[ghl_m1_neutrino_transport_component_count];
  candidate_diagnostics->opacity_suppression = 0.0;
  candidate_diagnostics->face_speed = face_speed;

  double opacity_suppression = 0.0;
  ghl_error_codes_t error = ghl_m1_compute_four_point_opacity_suppression(
        m1_params, kappa_face, delta_x, &opacity_suppression);
  if(error != ghl_success) {
    return error;
  }

  for(int component = 0; component < ghl_m1_neutrino_transport_component_count;
      ++component) {
    double phi = 0.0;
    bool sawtooth = false;
    const double dum = state_stencil[1][component] - state_stencil[0][component];
    const double duc = state_stencil[2][component] - state_stencil[1][component];
    const double dup = state_stencil[3][component] - state_stencil[2][component];
    if(isfinite(dum) && isfinite(duc) && isfinite(dup)) {
      error = ghl_m1_compute_four_point_flux_limiter(
            m1_params, dum, duc, dup, &phi, &sawtooth);
    }
    else {
      error = ghl_m1_compute_overflowing_stencil_limiter(
            m1_params, state_stencil[0][component], state_stencil[1][component],
            state_stencil[2][component], state_stencil[3][component], &phi,
            &sawtooth);
    }
    if(error != ghl_success) {
      return error;
    }

    const double flux_high
          = ghl_m1_average_finite(flux_L[component], flux_R[component]);
    if(!isfinite(flux_high)) {
      return ghl_error_m1_invalid_state;
    }

    double blended_flux = 0.0;
    error = ghl_m1_compute_four_point_blended_flux(
          flux_high, flux_low[component], phi, sawtooth, opacity_suppression,
          &blended_flux);
    if(error != ghl_success) {
      /* All arguments were validated above; a finite convex blend cannot
       * fail unless a future change breaks that invariant. */
      return error; // GCOVR_EXCL_LINE
    }

    candidate_flux[component] = blended_flux;
    candidate_diagnostics->phi[component] = phi;
    candidate_diagnostics->sawtooth[component] = sawtooth;
  }
  candidate_diagnostics->opacity_suppression = opacity_suppression;

  for(int component = 0; component < ghl_m1_neutrino_transport_component_count;
      ++component) {
    flux_tilde[component] = candidate_flux[component];
  }
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_neutrino_four_point_volume_weighted_transport_flux(
      const ghl_m1_parameters *restrict m1_params,
      const double state_stencil[4][ghl_m1_neutrino_transport_component_count],
      const double physical_flux_L[ghl_m1_neutrino_transport_component_count],
      const double physical_flux_R[ghl_m1_neutrino_transport_component_count],
      const double speed_L,
      const double speed_R,
      const double kappa_face,
      const double delta_x,
      const bool diffusion_correction_enabled,
      double flux_tilde[ghl_m1_neutrino_transport_component_count],
      ghl_m1_four_point_transport_diagnostics *restrict diagnostics) {
  ghl_m1_four_point_transport_diagnostics candidate_diagnostics;
  const ghl_error_codes_t error = ghl_m1_compute_neutrino_four_point_transport_core(
        m1_params, state_stencil, physical_flux_L, physical_flux_R, speed_L, speed_R,
        kappa_face, delta_x, diffusion_correction_enabled, flux_tilde,
        &candidate_diagnostics);
  if(error != ghl_success) {
    return error;
  }
  if(diagnostics != NULL) {
    *diagnostics = candidate_diagnostics;
  }
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_neutrino_four_point_transport_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric_face,
      const double state_stencil[4][ghl_m1_neutrino_transport_component_count],
      const double physical_flux_L[ghl_m1_neutrino_transport_component_count],
      const double physical_flux_R[ghl_m1_neutrino_transport_component_count],
      const double speed_L,
      const double speed_R,
      const double kappa_face,
      const double delta_x,
      const bool diffusion_correction_enabled,
      double flux_tilde[ghl_m1_neutrino_transport_component_count],
      ghl_m1_four_point_transport_diagnostics *restrict diagnostics) {
  /* Preserve the established pointwise validation order before invoking the
   * shared prepared-operand core. */
  if(m1_params == NULL || metric_face == NULL || state_stencil == NULL
     || physical_flux_L == NULL || physical_flux_R == NULL || flux_tilde == NULL) {
    return ghl_error_m1_null_pointer;
  }

  const ghl_error_codes_t policy_error
        = ghl_m1_validate_four_point_transport_policy(diffusion_correction_enabled);
  if(policy_error != ghl_success) {
    return policy_error;
  }

  if(!ghl_m1_metric_is_symmetric_spd(metric_face)
     || !isfinite(metric_face->sqrt_detgamma) || metric_face->sqrt_detgamma <= 0.0) {
    return ghl_error_m1_invalid_metric;
  }

  double blended_flux[ghl_m1_neutrino_transport_component_count];
  ghl_m1_four_point_transport_diagnostics candidate_diagnostics;
  const ghl_error_codes_t error = ghl_m1_compute_neutrino_four_point_transport_core(
        m1_params, state_stencil, physical_flux_L, physical_flux_R, speed_L, speed_R,
        kappa_face, delta_x, diffusion_correction_enabled, blended_flux,
        &candidate_diagnostics);
  if(error != ghl_success) {
    return error;
  }

  double candidate_flux[ghl_m1_neutrino_transport_component_count];
  for(int component = 0; component < ghl_m1_neutrino_transport_component_count;
      ++component) {
    const double densitized_flux = metric_face->sqrt_detgamma * blended_flux[component];
    if(!isfinite(densitized_flux)) {
      return ghl_error_m1_invalid_state;
    }
    candidate_flux[component] = densitized_flux;
  }

  for(int component = 0; component < ghl_m1_neutrino_transport_component_count;
      ++component) {
    flux_tilde[component] = candidate_flux[component];
  }
  if(diagnostics != NULL) {
    *diagnostics = candidate_diagnostics;
  }
  return ghl_success;
}
