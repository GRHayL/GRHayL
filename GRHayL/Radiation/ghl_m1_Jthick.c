#include "ghl_m1.h"
#include "ghl_m1_utils.h"

ghl_error_codes_t ghl_m1_compute_Jthick(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      double *restrict Jthick,
      bool *restrict Jthick_is_valid) {
  if(m1_params == NULL || metric == NULL || prims == NULL || rad_state == NULL ||
     Jthick == NULL || Jthick_is_valid == NULL)
    return ghl_error_m1_null_pointer;

  ghl_error_codes_t error = ghl_m1_validate_realizability(
      m1_params, metric, rad_state, 64.0, NULL);
  if(error != ghl_success)
    return error;

  double V_con[3];
  double W = 0.0;
  error = ghl_m1_compute_eulerian_velocity(
      metric, prims, V_con, NULL, &W);
  if(error != ghl_success)
    return error;

  const double FdotV = rad_state->F[0] * V_con[0]
                     + rad_state->F[1] * V_con[1]
                     + rad_state->F[2] * V_con[2];
  const double W2 = SQR(W);
  const double denom = 2.0 * W2 + 1.0;
  if(!isfinite(denom) || denom <= 0.0)
    return ghl_error_m1_invalid_state;

  const double prefactor = 3.0 / denom;
  const double bracket = (2.0 * W2 - 1.0) * rad_state->E - 2.0 * W2 * FdotV;
  const double J_local = prefactor * bracket;
  if(!isfinite(J_local)) {
    *Jthick = J_local;
    *Jthick_is_valid = false;
    return ghl_success;
  }

  *Jthick = J_local;
  *Jthick_is_valid = (J_local > 0.0);
  return ghl_success;
}
