#include "ghl_m1.h"
#include "ghl_m1_utils.h"
#include "ghl_flux_source.h"
#include <float.h>

ghl_error_codes_t ghl_m1_compute_physical_flux(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      double *restrict flux_E,
      double flux_F[3]) {
  if(metric_face == NULL || rad_state == NULL || closure == NULL ||
     flux_E == NULL || flux_F == NULL)
    return ghl_error_m1_null_pointer;
  ghl_error_codes_t error = ghl_m1_validate_direction(direction);
  if(error != ghl_success)
    return error;
  if(!ghl_m1_metric_is_symmetric_spd(metric_face))
    return ghl_error_m1_invalid_metric;
  error = ghl_m1_validate_closure_tensor(metric_face, rad_state, closure);
  if(error != ghl_success)
    return error;

  double F_faceU[3];
  ghl_raise_lower_vector_3D(metric_face->gammaUU, rad_state->F, F_faceU);

  const double candidate_E = metric_face->lapse * F_faceU[direction]
                           - metric_face->betaU[direction] * rad_state->E;
  if(!isfinite(candidate_E))
    return ghl_error_m1_invalid_state;

  double candidate_F[3];
  for(int i = 0; i < 3; i++) {
    double P_mixed = 0.0;
    for(int k = 0; k < 3; k++)
      P_mixed += closure->P[direction][k] * metric_face->gammaDD[k][i];

    candidate_F[i] = metric_face->lapse * P_mixed
                   - metric_face->betaU[direction] * rad_state->F[i];
    if(!isfinite(candidate_F[i]))
      return ghl_error_m1_invalid_state;
  }

  *flux_E = candidate_E;
  for(int i = 0; i < 3; ++i)
    flux_F[i] = candidate_F[i];
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_rusanov_flux(
      const ghl_m1_rad_state *restrict state_L,
      const ghl_m1_rad_state *restrict state_R,
      const double physical_flux_E_L,
      const double physical_flux_F_L[3],
      const double physical_flux_E_R,
      const double physical_flux_F_R[3],
      const double speed,
      double *restrict flux_E,
      double flux_F[3]) {
  if(state_L == NULL || state_R == NULL || physical_flux_F_L == NULL ||
     physical_flux_F_R == NULL || flux_E == NULL || flux_F == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(speed) || speed < 0.0 || !isfinite(physical_flux_E_L) ||
     !isfinite(physical_flux_E_R) || !isfinite(state_L->E) ||
     !isfinite(state_R->E))
    return ghl_error_m1_invalid_state;

  const double state_components_L[4] = {
      state_L->E, state_L->F[0], state_L->F[1], state_L->F[2] };
  const double state_components_R[4] = {
      state_R->E, state_R->F[0], state_R->F[1], state_R->F[2] };
  const double physical_flux_components_L[4] = {
      physical_flux_E_L, physical_flux_F_L[0], physical_flux_F_L[1],
      physical_flux_F_L[2] };
  const double physical_flux_components_R[4] = {
      physical_flux_E_R, physical_flux_F_R[0], physical_flux_F_R[1],
      physical_flux_F_R[2] };
  double candidate[4];
  const ghl_error_codes_t error = ghl_calculate_Rusanov_flux(
      state_components_L, state_components_R, physical_flux_components_L,
      physical_flux_components_R, 4, speed, candidate);
  if(error != ghl_success)
    return ghl_error_m1_invalid_state;

  *flux_E = candidate[0];
  for(int i = 0; i < 3; ++i)
    flux_F[i] = candidate[i + 1];
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_number_rusanov_flux(
      const double N_L,
      const double N_R,
      const double physical_flux_L,
      const double physical_flux_R,
      const double speed,
      double *restrict number_flux) {
  if(number_flux == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(N_L) || !isfinite(N_R) || !isfinite(physical_flux_L) ||
     !isfinite(physical_flux_R) || !isfinite(speed) || speed < 0.0)
    return ghl_error_m1_invalid_state;

  const double state_components_L[1] = { N_L };
  const double state_components_R[1] = { N_R };
  const double physical_flux_components_L[1] = { physical_flux_L };
  const double physical_flux_components_R[1] = { physical_flux_R };
  double candidate[1];
  const ghl_error_codes_t error = ghl_calculate_Rusanov_flux(
      state_components_L, state_components_R, physical_flux_components_L,
      physical_flux_components_R, 1, speed, candidate);
  if(error != ghl_success)
    return ghl_error_m1_invalid_state;
  *number_flux = candidate[0];
  return ghl_success;
}
