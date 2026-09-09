#include "ghl_m1.h"
#include "ghl_m1_utils.h"

static void ghl_m1_reset_optional_diagnostics(
      ghl_m1_diagnostics *restrict diagnostics) {

  diagnostics->realizability_repaired = false;
  diagnostics->Jthick = NAN;
  diagnostics->Jthick_is_valid = false;
  for(int side = 0; side < 2; side++) {
    for(int dir = 0; dir < 3; dir++)
      diagnostics->diffusion_blend_factor[side][dir] = NAN;
  }
}

static ghl_error_codes_t ghl_m1_validate_inputs(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure) {

  if(!ghl_m1_metric_is_symmetric_spd(metric))
    return ghl_error_m1_invalid_metric;

  ghl_error_codes_t error = ghl_m1_validate_realizability(
      m1_params, metric, rad_state, 64.0, NULL);
  if(error != ghl_success)
    return error;

  error = ghl_m1_validate_closure_tensor(metric, rad_state, closure);
  if(error != ghl_success)
    return error;

  if(!isfinite(closure->xi) || closure->xi < 0.0 || closure->xi > 1.0 ||
     !isfinite(closure->chi) || closure->chi < 1.0 / 3.0 ||
     closure->chi > 1.0 || !isfinite(closure->root_residual) ||
     closure->root_residual < 0.0 || closure->root_iterations < 0 ||
     (closure->solve_status != ghl_m1_closure_solve_converged &&
      closure->solve_status != ghl_m1_closure_solve_endpoint_fallback &&
      closure->solve_status != ghl_m1_closure_solve_iteration_exhausted))
    return ghl_error_m1_invalid_state;

  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_diagnostics *restrict diagnostics) {
  if(m1_params == NULL || metric == NULL || rad_state == NULL ||
     closure == NULL || diagnostics == NULL)
    return ghl_error_m1_null_pointer;

  ghl_error_codes_t error = ghl_m1_validate_inputs(m1_params, metric, rad_state, closure);
  if(error != ghl_success)
    return error;

  double flux_factor;
  error = ghl_m1_scaled_covector_norm_ratio(
      metric->gammaUU, rad_state->F, rad_state->E, &flux_factor);
  if(error != ghl_success)
    return error;
  const double r_from_state = flux_factor * flux_factor;
  /* r is the squared reduced-flux magnitude, so its admissible limit is
   * 1 - epsilon_c. */
  const double r_limit = 1.0 - m1_params->epsilon_c;
  const double r_diag = ghl_m1_min(ghl_m1_max(r_from_state, 0.0), r_limit);
  diagnostics->closure_xi = closure->xi;
  diagnostics->closure_root_residual = closure->root_residual;
  diagnostics->closure_root_iterations = closure->root_iterations;
  diagnostics->closure_solve_status = closure->solve_status;
  diagnostics->r = r_diag;
  diagnostics->chi_eddington = closure->chi;

  // These fields are workflow-owned and may be filled by callers after closure diagnostics.
  ghl_m1_reset_optional_diagnostics(diagnostics);

  return ghl_success;
}
