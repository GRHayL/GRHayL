#include "ghl_m1.h"
#include "../ghl_m1_utils.h"
#include "ghl_flux_source.h"
#include <float.h>

/*
 * Combined scalar N and E/F_i symmetric Rusanov flux across a single face.
 * The caller supplies one nonnegative, undensitized speed that is applied to
 * all five components in the order {N, E, Fx, Fy, Fz}; the returned fluxes are
 * densitized with metric_face->sqrt_detgamma. No HLL envelope, star-state,
 * speed cap, or diffusion intermediate is constructed.
 */

ghl_error_codes_t ghl_m1_compute_neutrino_rusanov_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const ghl_m1_neutrino_state *restrict state_L,
      const ghl_m1_neutrino_state *restrict state_R,
      const ghl_m1_closure *restrict closure_L,
      const ghl_m1_closure *restrict closure_R,
      const double number_flux_L[3],
      const double number_flux_R[3],
      const double number_transport_velocity_L[3],
      const double number_transport_velocity_R[3],
      const double speed,
      double *restrict flux_tildeN,
      double *restrict flux_tildeE,
      double flux_tildeF[3]) {
  if(m1_params == NULL || nu_params == NULL || metric_face == NULL ||
     state_L == NULL || state_R == NULL || closure_L == NULL ||
     closure_R == NULL || number_flux_L == NULL || number_flux_R == NULL ||
     number_transport_velocity_L == NULL ||
     number_transport_velocity_R == NULL || flux_tildeN == NULL ||
     flux_tildeE == NULL || flux_tildeF == NULL)
    return ghl_error_m1_null_pointer;

  if(direction < ghl_m1_dirn0 || direction > ghl_m1_dirn2 ||
     !isfinite(speed) || speed < 0.0)
    return ghl_error_m1_invalid_state;
  if(!ghl_m1_metric_is_symmetric_spd(metric_face))
    return ghl_error_m1_invalid_metric;
  if(!isfinite(nu_params->N_floor) || nu_params->N_floor < 0.0)
    return ghl_error_m1_invalid_state;

  const ghl_m1_rad_state rad_state_L =
      ghl_m1_neutrino_project_rad_state(state_L);
  const ghl_m1_rad_state rad_state_R =
      ghl_m1_neutrino_project_rad_state(state_R);
  ghl_error_codes_t error = ghl_m1_validate_realizability(
      m1_params, metric_face, &rad_state_L, 64.0, NULL);
  if(error != ghl_success)
    return error;
  error = ghl_m1_validate_realizability(
      m1_params, metric_face, &rad_state_R, 64.0, NULL);
  if(error != ghl_success)
    return error;
  error = ghl_m1_validate_closure_tensor(metric_face, &rad_state_L, closure_L);
  if(error != ghl_success)
    return error;
  error = ghl_m1_validate_closure_tensor(metric_face, &rad_state_R, closure_R);
  if(error != ghl_success)
    return error;
  if(!isfinite(state_L->N) || !isfinite(state_R->N) ||
     state_L->N < nu_params->N_floor || state_R->N < nu_params->N_floor)
    return ghl_error_m1_invalid_state;

  for(int i = 0; i < 3; ++i) {
    if(!isfinite(number_flux_L[i]) || !isfinite(number_flux_R[i]) ||
       !isfinite(number_transport_velocity_L[i]) ||
       !isfinite(number_transport_velocity_R[i]))
      return ghl_error_m1_invalid_state;
  }
  error = ghl_m1_validate_transport_velocity(
      metric_face, number_transport_velocity_L);
  if(error != ghl_success)
    return error;
  error = ghl_m1_validate_transport_velocity(
      metric_face, number_transport_velocity_R);
  if(error != ghl_success)
    return error;
  for(int i = 0; i < 3; ++i) {
    const double expected_L = state_L->N * number_transport_velocity_L[i];
    const double expected_R = state_R->N * number_transport_velocity_R[i];
    const double scale_L = ghl_m1_max(fabs(number_flux_L[i]), fabs(expected_L));
    const double scale_R = ghl_m1_max(fabs(number_flux_R[i]), fabs(expected_R));
    if(!isfinite(expected_L) || !isfinite(expected_R) ||
       (state_L->N == 0.0 ? number_flux_L[i] != 0.0 :
        fabs(number_flux_L[i] - expected_L) >
            128.0 * DBL_EPSILON * scale_L) ||
       (state_R->N == 0.0 ? number_flux_R[i] != 0.0 :
        fabs(number_flux_R[i] - expected_R) >
            128.0 * DBL_EPSILON * scale_R))
      return ghl_error_m1_invalid_state;
  }

  const int d = (int)direction;
  const double physical_number_flux_L =
      metric_face->lapse * number_flux_L[d]
      - metric_face->betaU[d] * state_L->N;
  const double physical_number_flux_R =
      metric_face->lapse * number_flux_R[d]
      - metric_face->betaU[d] * state_R->N;
  if(!isfinite(physical_number_flux_L) ||
     !isfinite(physical_number_flux_R))
    return ghl_error_m1_invalid_state;

  double physical_E_L = 0.0;
  double physical_E_R = 0.0;
  double physical_F_L[3];
  double physical_F_R[3];
  error = ghl_m1_compute_physical_flux(
      metric_face, direction, &rad_state_L, closure_L,
      &physical_E_L, physical_F_L);
  if(error != ghl_success)
    return error;
  error = ghl_m1_compute_physical_flux(
      metric_face, direction, &rad_state_R, closure_R,
      &physical_E_R, physical_F_R);
  if(error != ghl_success)
    return error;

  const double state_components_L[5] = {
      state_L->N, state_L->E, state_L->F[0], state_L->F[1], state_L->F[2] };
  const double state_components_R[5] = {
      state_R->N, state_R->E, state_R->F[0], state_R->F[1], state_R->F[2] };
  const double physical_flux_components_L[5] = {
      physical_number_flux_L, physical_E_L, physical_F_L[0],
      physical_F_L[1], physical_F_L[2] };
  const double physical_flux_components_R[5] = {
      physical_number_flux_R, physical_E_R, physical_F_R[0],
      physical_F_R[1], physical_F_R[2] };
  double flux_local[5];
  error = ghl_calculate_Rusanov_flux(
      state_components_L, state_components_R, physical_flux_components_L,
      physical_flux_components_R, 5, speed, flux_local);
  if(error != ghl_success)
    return ghl_error_m1_invalid_state;

  for(int component = 0; component < 5; ++component) {
    flux_local[component] *= metric_face->sqrt_detgamma;
    if(!isfinite(flux_local[component]))
      return ghl_error_m1_invalid_state;
  }

  *flux_tildeN = flux_local[0];
  *flux_tildeE = flux_local[1];
  for(int i = 0; i < 3; ++i)
    flux_tildeF[i] = flux_local[i + 2];
  return ghl_success;
}
