#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"
#include "../ghl_m1_utils.h"
#include <float.h>

/*
 * Internal residual helper for the grey three-species neutrino M1 implicit
 * E/F_i source solve. The coupling uses frozen rates and frozen primitives;
 * it does not call Con2Prim or opacity callbacks.
 *
 * The implicit unknowns are U[4] = (tildeE, tildeF_i), the densitized radiation
 * energy and covariant flux. For a supplied substep base U_base, the residual is
 *
 *   residual[0]    = tildeE    - U_base[0]   - dt_alpha * sqrt(gamma) * S_E
 *   residual[1..3] = tildeF_i  - U_base[1+i] - dt_alpha * sqrt(gamma) * S_i
 *
 * where S_E and S_i are the undensitized neutrino interaction sources
 * (ghl_m1_compute_neutrino_interaction_sources) evaluated at the trial
 * undensitized rad_state with FROZEN matter primitives and FROZEN rates.
 *
 * Frozen rates make eta_N, eta_E, kappa_a_N, kappa_a_E, kappa_s, kappa_tr,
 * n_eq, and J_eq constant throughout the Newton solve.
 */

static ghl_error_codes_t ghl_m1_neutrino_build_trial_state_core(
      const ghl_metric_quantities *restrict metric,
      const double U[4],
      ghl_m1_rad_state *restrict rad_state) {

  const double inv_sqrt_detgamma = 1.0 / metric->sqrt_detgamma;
  rad_state->E = U[0] * inv_sqrt_detgamma;
  for(int i = 0; i < 3; i++)
    rad_state->F[i] = U[i + 1] * inv_sqrt_detgamma;

  if(!isfinite(rad_state->E))
    return ghl_error_m1_implicit_admissibility;
  for(int i = 0; i < 3; i++) {
    if(!isfinite(rad_state->F[i]))
      return ghl_error_m1_implicit_admissibility;
  }

  return ghl_success;
}

ghl_error_codes_t ghl_m1_neutrino_build_trial_state(
      const ghl_metric_quantities *restrict metric,
      const double U[4],
      ghl_m1_rad_state *restrict rad_state) {

  if(metric == NULL || U == NULL || rad_state == NULL)
    return ghl_error_m1_null_pointer;
  if(!ghl_m1_metric_is_symmetric_spd(metric))
    return ghl_error_m1_invalid_metric;
  return ghl_m1_neutrino_build_trial_state_core(metric, U, rad_state);
}

static ghl_error_codes_t ghl_m1_neutrino_check_trial_admissibility_validated(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state) {

  const ghl_error_codes_t error = ghl_m1_validate_realizability_state(
      m1_params, metric, rad_state, 128.0, NULL);
  if(error == ghl_error_m1_invalid_state)
    return ghl_error_m1_implicit_admissibility;
  return error;
}

ghl_error_codes_t ghl_m1_neutrino_check_trial_admissibility(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state) {

  if(m1_params == NULL || metric == NULL || rad_state == NULL)
    return ghl_error_m1_null_pointer;

  const ghl_error_codes_t error = ghl_m1_validate_realizability(
      m1_params, metric, rad_state, 128.0, NULL);
  if(error == ghl_error_m1_invalid_state)
    return ghl_error_m1_implicit_admissibility;
  return error;
}

static ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_core(
      const ghl_m1_neutrino_implicit_context *restrict context,
      const double dt,
      const double U_base[4],
      const double U[4],
      const bool validate_configuration,
      const bool validate_rates,
      bool *restrict closure_fallback_observed,
      double residual[4]) {

  if(context == NULL || context->m1_params == NULL || context->metric == NULL ||
     context->prims_frozen == NULL || context->rates == NULL || U_base == NULL ||
     U == NULL || residual == NULL)
    return ghl_error_m1_null_pointer;
  const ghl_m1_parameters *restrict m1_params = context->m1_params;
  const ghl_metric_quantities *restrict metric = context->metric;
  const ghl_primitive_quantities *restrict prims_frozen = context->prims_frozen;
  const ghl_m1_neutrino_rates *restrict rates = context->rates;

  if(!isfinite(dt) || dt < 0.0)
    return ghl_error_m1_invalid_state;
  for(int i = 0; i < 4; i++) {
    if(!isfinite(U_base[i]))
      return ghl_error_m1_invalid_state;
  }

  ghl_m1_rad_state rad_state = {0};
  ghl_error_codes_t error = ghl_m1_neutrino_build_trial_state_core(
      metric, U, &rad_state);
  if(error != ghl_success)
    return error;

  if(validate_configuration) {
    /* The checked wrapper has already validated the metric. Keep the same
     * parameter error boundary without performing that metric walk again. */
    error = ghl_m1_validate_parameters(m1_params);
    if(error != ghl_success)
      return error;
  }
  error = ghl_m1_neutrino_check_trial_admissibility_validated(
      m1_params, metric, &rad_state);
  if(error != ghl_success)
    return error;

  if(validate_rates) {
    error = ghl_m1_neutrino_validate_single_species_rates(rates, NULL);
    if(error != ghl_success)
      return error;
  }

  /* The Newton projection is intentionally E/F-only: no number-current floor,
   * Gamma_N, or N_source is inspected here. */
  ghl_m1_sources EF_sources = {0};
  error = ghl_m1_neutrino_compute_EF_interaction_sources_validated(
      m1_params, metric, prims_frozen, &rad_state, rates,
      closure_fallback_observed, &EF_sources);
  if(error != ghl_success)
    return error;

  const double sqrt_detgamma = metric->sqrt_detgamma;
  const double dt_alpha = metric->lapse * dt;
  if(!isfinite(dt_alpha) || dt_alpha < 0.0)
    return ghl_error_m1_invalid_state;
  const double dt_alpha_sqrt_detgamma = dt_alpha * sqrt_detgamma;
  if(!isfinite(dt_alpha_sqrt_detgamma))
    return ghl_error_m1_invalid_state;

  double residual_local[4];
  residual_local[0] = U[0] - U_base[0]
                    - dt_alpha_sqrt_detgamma * EF_sources.S_E;
  for(int i = 0; i < 3; i++) {
    residual_local[i + 1] = U[i + 1] - U_base[i + 1]
                          - dt_alpha_sqrt_detgamma * EF_sources.S[i];
  }

  for(int i = 0; i < 4; i++) {
    if(!isfinite(residual_local[i]))
      return ghl_error_m1_invalid_state;
  }
  for(int i = 0; i < 4; i++)
    residual[i] = residual_local[i];

  return ghl_success;
}

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_validated(
      const ghl_m1_neutrino_implicit_context *restrict context,
      const double dt,
      const double U_base[4],
      const double U[4],
      bool *restrict closure_fallback_observed,
      double residual[4]) {
  return ghl_m1_neutrino_compute_implicit_residual_core(
      context, dt, U_base, U, false, false,
      closure_fallback_observed, residual);
}

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_with_base(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double U_base[4],
      const double U[4],
      double residual[4]) {

  return ghl_m1_neutrino_compute_implicit_residual_with_base_diagnostics(
      m1_params, metric, prims_frozen, rates, dt, U_base, U, NULL, residual);
}

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_with_base_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double U_base[4],
      const double U[4],
      bool *restrict closure_fallback_observed,
      double residual[4]) {

  if(m1_params == NULL || metric == NULL ||
     prims_frozen == NULL || rates == NULL ||
     U_base == NULL || U == NULL || residual == NULL)
    return ghl_error_m1_null_pointer;

  if(!isfinite(dt) || dt < 0.0)
    return ghl_error_m1_invalid_state;

  if(!ghl_m1_metric_is_symmetric_spd(metric))
    return ghl_error_m1_invalid_metric;

  for(int i = 0; i < 4; i++) {
    if(!isfinite(U_base[i]))
      return ghl_error_m1_invalid_state;
  }
  const ghl_m1_neutrino_implicit_context context = {
        .m1_params = m1_params,
        .metric = metric,
        .prims_frozen = prims_frozen,
        .rates = rates };
  return ghl_m1_neutrino_compute_implicit_residual_core(
      &context, dt, U_base, U, true, true,
      closure_fallback_observed, residual);
}

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in,
      const double dt,
      const double U[4],
      double residual[4]) {

  if(metric == NULL || state_in == NULL)
    return ghl_error_m1_null_pointer;

  if(!isfinite(metric->sqrt_detgamma) || metric->sqrt_detgamma <= 0.0)
    return ghl_error_m1_invalid_metric;

  const double sqrt_detgamma = metric->sqrt_detgamma;
  const double U_base[4] = {
    state_in->E * sqrt_detgamma,
    state_in->F[0] * sqrt_detgamma,
    state_in->F[1] * sqrt_detgamma,
    state_in->F[2] * sqrt_detgamma
  };

  if(nu_params == NULL)
    return ghl_error_m1_null_pointer;

  return ghl_m1_neutrino_compute_implicit_residual_with_base(
      m1_params, metric, prims_frozen, rates, dt, U_base, U, residual);
}
