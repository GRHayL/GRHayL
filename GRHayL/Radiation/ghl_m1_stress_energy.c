#include "ghl_m1.h"
#include "ghl_m1_utils.h"

static ghl_error_codes_t ghl_m1_validate_inputs(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure) {

  ghl_error_codes_t error = ghl_m1_validate_realizability(
      m1_params, metric, rad_state, 64.0, NULL);
  if(error != ghl_success)
    return error;

  return ghl_m1_validate_closure_tensor(metric, rad_state, closure);
}

ghl_error_codes_t ghl_m1_compute_stress_energy(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_stress_energy *restrict Rmunu) {
  if(m1_params == NULL || metric == NULL || rad_state == NULL || closure == NULL || Rmunu == NULL)
    return ghl_error_m1_null_pointer;

  ghl_error_codes_t error = ghl_m1_validate_inputs(m1_params, metric, rad_state, closure);
  if(error != ghl_success)
    return error;

  if(metric->lapse <= 0.0)
    return ghl_error_m1_invalid_metric;

  const double inv_alpha = 1.0 / metric->lapse;
  const double inv_alpha_sq = SQR(inv_alpha);
  if(!isfinite(inv_alpha_sq))
    return ghl_error_m1_invalid_metric;

  const double E = rad_state->E;

  double F_con[3];
  ghl_raise_lower_vector_3D(metric->gammaUU, rad_state->F, F_con);

  double beta_over_alpha[3];
  for(int i = 0; i < 3; i++)
    beta_over_alpha[i] = metric->betaU[i] * inv_alpha;

  ghl_stress_energy Rmunu_local = {0};
  for(int mu = 0; mu < 4; mu++) {
    for(int nu = 0; nu < 4; nu++)
      Rmunu_local.T4[mu][nu] = 0.0;
  }

  Rmunu_local.T4[0][0] = E * inv_alpha_sq;

  for(int i = 0; i < 3; i++) {
    const double R0i = inv_alpha * (F_con[i] - beta_over_alpha[i] * E);
    Rmunu_local.T4[0][i + 1] = R0i;
    Rmunu_local.T4[i + 1][0] = R0i;
  }

  for(int i = 0; i < 3; i++) {
    for(int j = i; j < 3; j++) {
      const double Rij = closure->P[i][j]
                       - beta_over_alpha[i] * F_con[j]
                       - beta_over_alpha[j] * F_con[i]
                       + beta_over_alpha[i] * beta_over_alpha[j] * E;
      Rmunu_local.T4[i + 1][j + 1] = Rij;
      Rmunu_local.T4[j + 1][i + 1] = Rij;
    }
  }

  for(int mu = 0; mu < 4; mu++) {
    for(int nu = 0; nu < 4; nu++) {
      if(!isfinite(Rmunu_local.T4[mu][nu]))
        return ghl_error_m1_invalid_state;
    }
  }

  *Rmunu = Rmunu_local;

  return ghl_success;
}
