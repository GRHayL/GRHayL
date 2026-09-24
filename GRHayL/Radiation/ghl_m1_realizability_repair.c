#include "ghl_m1.h"
#include "ghl_m1_utils.h"

ghl_error_codes_t ghl_m1_apply_energy_floor(
      const ghl_m1_parameters *restrict m1_params,
      const double E_in,
      double *restrict E_out,
      bool *restrict floor_applied) {
  if(m1_params == NULL || E_out == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(m1_params->E_floor) || m1_params->E_floor <= 0.0) {
    return ghl_error_m1_invalid_E_floor;
  }
  if(!isfinite(E_in)) {
    return ghl_error_m1_invalid_state;
  }
  const bool applied = E_in < m1_params->E_floor;
  *E_out = applied ? m1_params->E_floor : E_in;
  if(floor_applied != NULL) {
    *floor_applied = applied;
  }
  return ghl_success;
}

ghl_error_codes_t ghl_m1_realizability_repair(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      ghl_m1_rad_state *restrict rad_state) {

  if(m1_params == NULL || metric == NULL || rad_state == NULL) {
    return ghl_error_m1_null_pointer;
  }

  ghl_error_codes_t error = ghl_m1_validate_configuration(m1_params, metric);
  if(error != ghl_success) {
    return error;
  }

  if(!isfinite(rad_state->E) || !isfinite(rad_state->F[0]) || !isfinite(rad_state->F[1])
     || !isfinite(rad_state->F[2])) {
    return ghl_error_m1_invalid_state;
  }

  double E_local;
  /* Configuration, finite energy and the local output satisfy every check
   * in the public floor helper. */
  (void)ghl_m1_apply_energy_floor(m1_params, rad_state->E, &E_local, NULL);
  double F_local[3] = { rad_state->F[0], rad_state->F[1], rad_state->F[2] };
  const double cone_factor = 1.0 - m1_params->epsilon_c;
  /* The canonical repair uses the squared-ratio flux rescale. */
  /* The factor below is a norm ratio.  The squared-ratio rescale places the
   * repaired state at or below sqrt(cone_factor), not cone_factor itself. */
  const double permitted_flux_factor = sqrt(cone_factor);

  double flux_factor;
  error = ghl_m1_scaled_covector_norm_ratio(
        metric->gammaUU, F_local, E_local, &flux_factor);
  if(error != ghl_success) {
    return error;
  }
  const bool flux_rescaled = flux_factor > permitted_flux_factor;
  if(flux_rescaled) {
    const double applied_scale = (cone_factor / flux_factor) / flux_factor;
    /* cone_factor is in (0,1], and flux_factor is finite and greater
     * than sqrt(cone_factor). Both divisions yield a finite scale in [0,1];
     * underflow to zero is the permitted complete flux repair. */
    for(int i = 0; i < 3; ++i) {
      F_local[i] *= applied_scale;
    }

    error = ghl_m1_scaled_covector_norm_ratio(
          metric->gammaUU, F_local, E_local, &flux_factor);
    if(error != ghl_success) {
      return error;
    }
  }

  const ghl_m1_rad_state repaired
        = { .E = E_local, .F = { F_local[0], F_local[1], F_local[2] } };
  /* sqrt(1-epsilon_c) is at most one for a validated epsilon_c. */
  const double scale = 1.0;
  if(flux_factor > permitted_flux_factor + 64.0 * DBL_EPSILON * scale) {
    return ghl_error_m1_invalid_state;
  }

  if(repaired.E != rad_state->E || repaired.F[0] != rad_state->F[0]
     || repaired.F[1] != rad_state->F[1] || repaired.F[2] != rad_state->F[2]) {
    ghl_m1_record_closure_downstream_repair();
  }
  *rad_state = repaired;
  return ghl_success;
}
