#include "ghl_m1.h"
#include "ghl_m1_utils.h"

static bool ghl_m1_validate_metric_derivatives(
      const ghl_metric_quantities *restrict metric_derivs) {

  if(!isfinite(metric_derivs->lapse))
    return false;
  for(int j = 0; j < 3; j++) {
    if(!isfinite(metric_derivs->betaU[j]))
      return false;
    for(int k = 0; k < 3; k++) {
      if(!isfinite(metric_derivs->gammaDD[j][k]))
        return false;
    }
  }
  return true;
}

static ghl_error_codes_t ghl_m1_validate_inputs(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_derivs_x,
      const ghl_metric_quantities *restrict metric_derivs_y,
      const ghl_metric_quantities *restrict metric_derivs_z,
      const ghl_extrinsic_curvature *restrict curv,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure) {

  if(!ghl_m1_metric_is_symmetric_spd(metric))
    return ghl_error_m1_invalid_metric;
  if(!ghl_m1_validate_metric_derivatives(metric_derivs_x) ||
     !ghl_m1_validate_metric_derivatives(metric_derivs_y) ||
     !ghl_m1_validate_metric_derivatives(metric_derivs_z))
    return ghl_error_m1_invalid_metric;

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(!isfinite(curv->K[i][j]))
        return ghl_error_m1_invalid_metric;
    }
  }

  ghl_error_codes_t error = ghl_m1_validate_realizability(
      m1_params, metric, rad_state, 64.0, NULL);
  if(error != ghl_success)
    return error;

  return ghl_m1_validate_closure_tensor(metric, rad_state, closure);
}

ghl_error_codes_t ghl_m1_compute_geometry_sources(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_derivs_x,
      const ghl_metric_quantities *restrict metric_derivs_y,
      const ghl_metric_quantities *restrict metric_derivs_z,
      const ghl_extrinsic_curvature *restrict curv,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_sources *restrict geometry_sources) {
  if(m1_params == NULL || metric == NULL ||
     metric_derivs_x == NULL || metric_derivs_y == NULL || metric_derivs_z == NULL ||
     curv == NULL || rad_state == NULL || closure == NULL || geometry_sources == NULL)
    return ghl_error_m1_null_pointer;

  ghl_error_codes_t error = ghl_m1_validate_inputs(
      m1_params, metric, metric_derivs_x, metric_derivs_y, metric_derivs_z, curv, rad_state, closure);
  if(error != ghl_success)
    return error;

  const ghl_metric_quantities *metric_derivs[3] = {
      metric_derivs_x, metric_derivs_y, metric_derivs_z};

  double F_con[3];
  ghl_raise_lower_vector_3D(metric->gammaUU, rad_state->F, F_con);

  const double alpha = metric->lapse;
  const double sqrt_detgamma = metric->sqrt_detgamma;

  double P_contracted_K = 0.0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++)
      P_contracted_K += closure->P[i][j] * curv->K[i][j];
  }

  double F_dot_dalpha = 0.0;
  for(int j = 0; j < 3; j++)
    F_dot_dalpha += F_con[j] * metric_derivs[j]->lapse;

  /* Cancel alpha * d_i(log(alpha)) analytically. Dividing by a floored
   * lapse changes the continuum source in strong-field, low-lapse cells. */
  const double S_E = sqrt_detgamma
      * (alpha * P_contracted_K - F_dot_dalpha);
  if(!isfinite(S_E))
    return ghl_error_m1_invalid_state;

  double S[3];
  for(int i = 0; i < 3; i++) {
    const ghl_metric_quantities *restrict d_i = metric_derivs[i];
    double Fj_dibeta_j = 0.0;
    for(int j = 0; j < 3; j++)
      Fj_dibeta_j += rad_state->F[j] * d_i->betaU[j];

    double Pjk_digamma_jk = 0.0;
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++)
        Pjk_digamma_jk += closure->P[j][k] * d_i->gammaDD[j][k];
    }

    S[i] = sqrt_detgamma * (-rad_state->E * d_i->lapse
                          + Fj_dibeta_j
                          + 0.5 * alpha * Pjk_digamma_jk);
    if(!isfinite(S[i]))
      return ghl_error_m1_invalid_state;
  }

  geometry_sources->S_E = S_E;
  for(int i = 0; i < 3; i++)
    geometry_sources->S[i] = S[i];

  return ghl_success;
}
