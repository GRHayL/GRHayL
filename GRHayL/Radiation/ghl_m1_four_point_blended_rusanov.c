#include "ghl_m1.h"
#include "ghl_m1_utils.h"
#include "ghl_flux_source.h"

static bool ghl_m1_same_nonzero_sign(const double left, const double right) {
  return (left > 0.0 && right > 0.0) ||
         (left < 0.0 && right < 0.0);
}

static bool ghl_m1_opposite_nonzero_sign(const double left, const double right) {
  return (left > 0.0 && right < 0.0) ||
         (left < 0.0 && right > 0.0);
}

static ghl_error_codes_t ghl_m1_validate_four_point_transport_policy(
      const bool diffusion_correction_enabled) {
  /* The separate diffusion correction is not part of this transport
   * operator. Rejecting it here prevents a caller from accidentally
   * obtaining a different face flux. */
  if(diffusion_correction_enabled)
    return ghl_error_m1_incompatible_transport_policy;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_four_point_flux_limiter(
      const ghl_m1_parameters *restrict m1_params,
      const double dum,
      const double duc,
      const double dup,
      double *restrict phi,
      bool *restrict sawtooth) {
  if(m1_params == NULL || phi == NULL || sawtooth == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(dum) || !isfinite(duc) || !isfinite(dup) ||
     !isfinite(m1_params->minmod_theta) ||
     m1_params->minmod_theta < 0.0 || m1_params->minmod_theta > 2.0)
    return ghl_error_m1_invalid_state;

  double candidate_phi = 0.0;
  bool candidate_sawtooth = false;
  if(ghl_m1_same_nonzero_sign(dup, duc) &&
     ghl_m1_same_nonzero_sign(dum, duc)) {
    const double ratio_left = dum / duc;
    const double ratio_right = dup / duc;
    candidate_phi = fmin(
        1.0, fmin(m1_params->minmod_theta * ratio_left,
                  m1_params->minmod_theta * ratio_right));
  } else if(ghl_m1_opposite_nonzero_sign(dup, duc) &&
            ghl_m1_opposite_nonzero_sign(dum, duc)) {
    candidate_sawtooth = true;
  }

  if(!isfinite(candidate_phi))
    return ghl_error_m1_invalid_state;
  *phi = candidate_phi;
  *sawtooth = candidate_sawtooth;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_four_point_opacity_suppression(
      const ghl_m1_parameters *restrict m1_params,
      const double kappa_face,
      const double delta,
      double *restrict A) {
  if(m1_params == NULL || A == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(kappa_face) || kappa_face < 0.0 ||
     !isfinite(delta) || delta <= 0.0 ||
     !isfinite(m1_params->mindiss) || m1_params->mindiss < 0.0 ||
     m1_params->mindiss > 1.0)
    return ghl_error_m1_invalid_state;

  double candidate_A = 1.0;
  const double optical_width = delta * kappa_face;
  if(optical_width > 1.0) {
    candidate_A = fmin(1.0, 1.0 / optical_width);
    candidate_A = fmax(candidate_A, m1_params->mindiss);
  }
  if(!isfinite(candidate_A))
    return ghl_error_m1_invalid_state;
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
  if(flux_num == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(flux_high) || !isfinite(flux_low) ||
     !isfinite(phi) || phi < 0.0 || phi > 1.0 ||
     !isfinite(A) || A < 0.0 || A > 1.0)
    return ghl_error_m1_invalid_state;

  const double dissipation = sawtooth ? 1.0 : A;
  const double candidate_flux = flux_high - dissipation * (1.0 - phi) *
      (flux_high - flux_low);
  if(!isfinite(candidate_flux))
    return ghl_error_m1_invalid_state;
  *flux_num = candidate_flux;
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
  /*
   * This is the canonical prepared, pointwise operation.  The host owns the
   * stencil assembly and supplies the physical fluxes and uncapped adjacent
   * light-cone speeds.  Keep all candidates local until every component has
   * passed validation so that a failed face operation cannot publish a
   * partially updated flux or diagnostic packet.
   */
  if(m1_params == NULL || metric_face == NULL || state_stencil == NULL ||
     physical_flux_L == NULL || physical_flux_R == NULL || flux_tilde == NULL)
    return ghl_error_m1_null_pointer;

  const ghl_error_codes_t policy_error =
      ghl_m1_validate_four_point_transport_policy(
          diffusion_correction_enabled);
  if(policy_error != ghl_success)
    return policy_error;

  if(!ghl_m1_metric_is_symmetric_spd(metric_face)
     || !isfinite(metric_face->sqrt_detgamma)
     || metric_face->sqrt_detgamma <= 0.0)
    return ghl_error_m1_invalid_metric;
  if(!isfinite(speed_L) || speed_L < 0.0 || !isfinite(speed_R) ||
     speed_R < 0.0 || !isfinite(kappa_face) || kappa_face < 0.0 ||
     !isfinite(delta_x) || delta_x <= 0.0)
    return ghl_error_m1_invalid_state;

  const double face_speed = fmax(speed_L, speed_R);
  if(!isfinite(face_speed) || face_speed < 0.0)
    return ghl_error_m1_invalid_state;

  /* Validate all prepared operands before invoking the generic Rusanov
   * arithmetic.  The latter owns the component-wise low-flux formula used by
   * the ordinary GRHayL transport wrappers. */
  double state_L[ghl_m1_neutrino_transport_component_count];
  double state_R[ghl_m1_neutrino_transport_component_count];
  double flux_L[ghl_m1_neutrino_transport_component_count];
  double flux_R[ghl_m1_neutrino_transport_component_count];
  for(int component = 0;
      component < ghl_m1_neutrino_transport_component_count; ++component) {
    state_L[component] = state_stencil[1][component];
    state_R[component] = state_stencil[2][component];
    flux_L[component] = physical_flux_L[component];
    flux_R[component] = physical_flux_R[component];
    if(!isfinite(state_stencil[0][component]) ||
       !isfinite(state_stencil[1][component]) ||
       !isfinite(state_stencil[2][component]) ||
       !isfinite(state_stencil[3][component]) ||
       !isfinite(flux_L[component]) || !isfinite(flux_R[component]))
      return ghl_error_m1_invalid_state;
  }

  double flux_low[ghl_m1_neutrino_transport_component_count];
  const ghl_error_codes_t rusanov_error = ghl_calculate_Rusanov_flux(
      state_L, state_R, flux_L, flux_R,
      ghl_m1_neutrino_transport_component_count, face_speed, flux_low);
  if(rusanov_error != ghl_success)
    return ghl_error_m1_invalid_state;

  double candidate_flux[ghl_m1_neutrino_transport_component_count];
  ghl_m1_four_point_transport_diagnostics candidate_diagnostics;
  candidate_diagnostics.opacity_suppression = 0.0;
  candidate_diagnostics.face_speed = face_speed;

  double opacity_suppression = 0.0;
  ghl_error_codes_t error = ghl_m1_compute_four_point_opacity_suppression(
      m1_params, kappa_face, delta_x, &opacity_suppression);
  if(error != ghl_success)
    return error;

  for(int component = 0;
      component < ghl_m1_neutrino_transport_component_count; ++component) {
    double phi = 0.0;
    bool sawtooth = false;
    error = ghl_m1_compute_four_point_flux_limiter(
        m1_params,
        state_stencil[1][component] - state_stencil[0][component],
        state_stencil[2][component] - state_stencil[1][component],
        state_stencil[3][component] - state_stencil[2][component],
        &phi, &sawtooth);
    if(error != ghl_success)
      return error;

    const double flux_high = 0.5 * (flux_L[component] + flux_R[component]);
    if(!isfinite(flux_high))
      return ghl_error_m1_invalid_state;

    double blended_flux = 0.0;
    error = ghl_m1_compute_four_point_blended_flux(
        flux_high, flux_low[component], phi, sawtooth,
        opacity_suppression, &blended_flux);
    if(error != ghl_success)
      return error;

    const double densitized_flux = metric_face->sqrt_detgamma * blended_flux;
    if(!isfinite(densitized_flux))
      return ghl_error_m1_invalid_state;
    candidate_flux[component] = densitized_flux;
    candidate_diagnostics.phi[component] = phi;
    candidate_diagnostics.sawtooth[component] = sawtooth;
  }
  candidate_diagnostics.opacity_suppression = opacity_suppression;

  for(int component = 0;
      component < ghl_m1_neutrino_transport_component_count; ++component)
    flux_tilde[component] = candidate_flux[component];
  if(diagnostics != NULL)
    *diagnostics = candidate_diagnostics;
  return ghl_success;
}
