#include "ghl_m1.h"
#include "ghl_m1_utils.h"

static ghl_error_codes_t ghl_m1_validate_inputs(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      const bool configuration_validated) {

  if(!configuration_validated && !ghl_m1_metric_is_symmetric_spd(metric))
    return ghl_error_m1_invalid_metric;

  ghl_error_codes_t error = configuration_validated
      ? ghl_m1_validate_realizability_state(m1_params, metric, rad_state, 64.0, NULL)
      : ghl_m1_validate_realizability(m1_params, metric, rad_state, 64.0, NULL);
  if(error != ghl_success)
    return error;

  return ghl_m1_validate_closure_tensor(metric, rad_state, closure);
}

static ghl_error_codes_t ghl_m1_compute_comoving_moments_internal(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving,
      double V_con[3],
      double V_cov[3],
      double *restrict W_out,
      const bool configuration_validated) {
  if(m1_params == NULL || metric == NULL || prims == NULL ||
     rad_state == NULL || closure == NULL || comoving == NULL ||
     V_con == NULL || V_cov == NULL || W_out == NULL)
    return ghl_error_m1_null_pointer;
  ghl_error_codes_t error = ghl_m1_validate_inputs(
      m1_params, metric, rad_state, closure, configuration_validated);
  if(error != ghl_success)
    return error;

  double W = 0.0;
  error = ghl_m1_compute_eulerian_velocity(metric, prims, V_con, NULL, &W);
  if(error != ghl_success)
    return error;

  ghl_raise_lower_vector_3D(metric->gammaDD, V_con, V_cov);

  double F_con[3];
  ghl_raise_lower_vector_3D(metric->gammaUU, rad_state->F, F_con);
  const double FdotV = rad_state->F[0] * V_con[0]
                     + rad_state->F[1] * V_con[1]
                     + rad_state->F[2] * V_con[2];

  double P_DD[3][3];
  ghl_m1_lower_spatial_tensor(metric, closure->P, P_DD);
  const double PVV = P_DD[0][0] * V_con[0] * V_con[0]
                   + P_DD[1][1] * V_con[1] * V_con[1]
                   + P_DD[2][2] * V_con[2] * V_con[2]
                   + 2.0 * (P_DD[0][1] * V_con[0] * V_con[1]
                          + P_DD[0][2] * V_con[0] * V_con[2]
                          + P_DD[1][2] * V_con[1] * V_con[2]);

  ghl_m1_comoving candidate = {0};
  const double J = SQR(W) * (rad_state->E - 2.0 * FdotV + PVV);
  // Keep this admissibility guard in all builds: negative/non-finite comoving
  // energy density is unphysical and must be rejected deterministically.
  if(!isfinite(J) || J < 0.0)
    return ghl_error_m1_invalid_state;
  candidate.J = J;

  for(int i = 0; i < 3; i++) {
    const double PijVj = closure->P[i][0] * V_cov[0]
                       + closure->P[i][1] * V_cov[1]
                       + closure->P[i][2] * V_cov[2];
    candidate.HU[i] = W * (F_con[i] - PijVj - J * V_con[i]);
    if(!isfinite(candidate.HU[i]))
      return ghl_error_m1_invalid_state;
  }

  ghl_raise_lower_vector_3D(metric->gammaDD, candidate.HU, candidate.HD);
  for(int i = 0; i < 3; i++) {
    if(!isfinite(candidate.HD[i]))
      return ghl_error_m1_invalid_state;
  }

  candidate.Hn = -(V_cov[0] * candidate.HU[0]
                 + V_cov[1] * candidate.HU[1]
                 + V_cov[2] * candidate.HU[2]);
  if(!isfinite(candidate.Hn))
    return ghl_error_m1_invalid_state;

  *comoving = candidate;
  *W_out = W;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_comoving_moments_with_velocity(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving,
      double V_con[3],
      double V_cov[3],
      double *restrict W_out) {
  return ghl_m1_compute_comoving_moments_internal(
      m1_params, metric, prims, rad_state, closure, comoving,
      V_con, V_cov, W_out, false);
}

ghl_error_codes_t ghl_m1_compute_comoving_moments_validated(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving,
      double V_con[3],
      double V_cov[3],
      double *restrict W_out) {
  return ghl_m1_compute_comoving_moments_internal(
      m1_params, metric, prims, rad_state, closure, comoving,
      V_con, V_cov, W_out, true);
}

ghl_error_codes_t ghl_m1_compute_comoving_moments(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving) {
  double V_con[3], V_cov[3], W;
  return ghl_m1_compute_comoving_moments_with_velocity(
      m1_params, metric, prims, rad_state, closure, comoving,
      V_con, V_cov, &W);
}
