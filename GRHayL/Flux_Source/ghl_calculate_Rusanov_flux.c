#include "ghl_flux_source.h"

#include <math.h>

static double rusanov_average_finite(const double left, const double right) {
  const double sum = left + right;
  /* Preserve the usual rounding unless the sum overflows before halving. */
  return isfinite(sum) ? 0.5 * sum : 0.5 * left + 0.5 * right;
}

/* Retry a nonfinite finite-jump candidate at half scale. The two low parts
 * retain cancellation that a rounded physical-flux average can otherwise
 * lose at the edge of double range. */
static double rusanov_half_scaled_candidate(
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

static double rusanov_candidate_finite(const double state_L, const double state_R,
                                       const double physical_flux_L,
                                       const double physical_flux_R,
                                       const double speed) {
  const double average = rusanov_average_finite(physical_flux_L, physical_flux_R);
  const double jump = state_R - state_L;
  if(isfinite(jump)) {
    /* Halve the jump when halving a subnormal speed would round to zero. */
    if(speed > 0.0 && 0.5 * speed == 0.0) {
      return fma(-speed, 0.5 * jump, average);
    }
    const double candidate = average - 0.5 * speed * jump;
    /* A fused product and sum can remain finite when the product overflows. */
    if(isfinite(candidate)) {
      return candidate;
    }
    const double fused = fma(-0.5 * speed, jump, average);
    return isfinite(fused) ? fused
                           : rusanov_half_scaled_candidate(
                                   physical_flux_L, physical_flux_R, speed,
                                   jump, 0.0, false);
  }
  /* The jump alone overflows; its half and the rounding error of that half
   * are both representable. Keep the low part for cancellation against the
   * physical-flux average. */
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
  /* The first fused result or the rounded average may overflow despite a
   * finite exact result. Retry the complete expression at half scale. */
  return rusanov_half_scaled_candidate(
        physical_flux_L, physical_flux_R, speed, half_jump, low_part, true);
}

ghl_error_codes_t ghl_calculate_Rusanov_flux(
      const double *restrict state_L,
      const double *restrict state_R,
      const double *restrict physical_flux_L,
      const double *restrict physical_flux_R,
      const int component_count,
      const double speed,
      double *restrict flux) {
  if(state_L == NULL || state_R == NULL || physical_flux_L == NULL
     || physical_flux_R == NULL || flux == NULL || component_count <= 0
     || !isfinite(speed) || speed < 0.0) {
    return ghl_error_flux_source_invalid_input;
  }

  /* Validate and evaluate every candidate before publishing any component. */
  for(int component = 0; component < component_count; ++component) {
    if(!isfinite(state_L[component]) || !isfinite(state_R[component])
       || !isfinite(physical_flux_L[component])
       || !isfinite(physical_flux_R[component])) {
      return ghl_error_flux_source_invalid_input;
    }

    const double candidate = rusanov_candidate_finite(
          state_L[component], state_R[component], physical_flux_L[component],
          physical_flux_R[component], speed);
    if(!isfinite(candidate)) {
      return ghl_error_flux_source_invalid_input;
    }
  }

  for(int component = 0; component < component_count; ++component) {
    flux[component] = rusanov_candidate_finite(
          state_L[component], state_R[component], physical_flux_L[component],
          physical_flux_R[component], speed);
  }
  return ghl_success;
}
