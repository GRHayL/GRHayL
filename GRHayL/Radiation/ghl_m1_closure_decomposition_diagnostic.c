#include "ghl_m1.h"
#include "ghl_m1_utils.h"

ghl_error_codes_t ghl_m1_compute_closure_decomposition_diagnostic(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure_decomposition_diagnostic *restrict diagnostic) {
  if(m1_params == NULL || metric == NULL || prims == NULL ||
     rad_state == NULL || diagnostic == NULL)
    return ghl_error_m1_null_pointer;

  ghl_m1_closure closure = {0};
  ghl_error_codes_t error = ghl_m1_compute_closure_with_primitives(
      m1_params, metric, prims, rad_state, &closure);
  if(error != ghl_success)
    return error;

  /* The public diagnostic field is explicitly H_mu H^mu / J^2. Minerbo's
   * solved xi is the reduced-flux magnitude, so its square belongs in that
   * field even though chi is evaluated from the unsquared magnitude. */
  const double xi_squared = SQR(closure.xi);
  const double chi = closure.chi;
  if(!isfinite(chi))
    return ghl_error_m1_invalid_state;

  diagnostic->chi_minerbo_xi = chi;
  diagnostic->xi_HaHa_over_J2 = xi_squared;
  diagnostic->dthin_scalar = 0.5 * (3.0 * chi - 1.0);
  diagnostic->dthick_scalar = 1.5 * (1.0 - chi);

  double Pthin_UU[3][3], Pthick_UU[3][3];
  error = ghl_m1_compute_minerbo_decomposition(
      m1_params, metric, prims, rad_state, Pthin_UU, Pthick_UU);
  if(error != ghl_success)
    return error;
  double Pthin_DD[3][3], Pthick_DD[3][3];
  ghl_m1_lower_spatial_tensor(metric, Pthin_UU, Pthin_DD);
  ghl_m1_lower_spatial_tensor(metric, Pthick_UU, Pthick_DD);

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      diagnostic->Pthin_dd[i][j] = Pthin_DD[i][j];
      diagnostic->Pthick_dd[i][j] = Pthick_DD[i][j];
    }
    diagnostic->Pthick_minus_Pthin_dd[i] = diagnostic->Pthick_dd[i][i]
                                         - diagnostic->Pthin_dd[i][i];
  }

  double Pthick_UU_trace = 0.0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++)
      Pthick_UU_trace += metric->gammaUU[i][j] * diagnostic->Pthick_dd[i][j];
  }
  diagnostic->Pth_dd_3_3_UU = Pthick_UU_trace / 3.0;
  diagnostic->Pth_dd_0_0_DD = diagnostic->Pthick_dd[0][0];

  return ghl_success;
}
