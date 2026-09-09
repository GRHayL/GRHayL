#include "ghl_m1.h"
#include "ghl_m1_utils.h"

#include <float.h>

static bool ghl_m1_face_velocity_is_valid(
      const ghl_metric_quantities *restrict metric_face,
      const double W_face,
      const double V_face[3]) {
  if(metric_face == NULL || !isfinite(W_face) || W_face < 1.0 || V_face == NULL)
    return false;
  for(int i = 0; i < 3; ++i) {
    if(!isfinite(V_face[i]))
      return false;
  }
  const double v2 = ghl_compute_vec2_from_vec3D(metric_face->gammaDD, V_face);
  if(!isfinite(v2) || v2 < 0.0 || v2 >= 1.0)
    return false;
  const double W_expected = 1.0 / sqrt(fmax(1.0 - v2, DBL_MIN));
  if(!isfinite(W_expected))
    return false;
  const double scale = fmax(1.0, fmax(fabs(W_face), fabs(W_expected)));
  const double tol = fmax(1024.0 * DBL_EPSILON * scale, 1e-8 * scale);
  return fabs(W_face - W_expected) <= tol;
}

ghl_error_codes_t ghl_m1_compute_diffusion_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const double hll_flux_tildeE,
      const double E_star,
      const double Jthick_L,
      const bool Jthick_L_is_valid,
      const double Jthick_R,
      const bool Jthick_R_is_valid,
      const double gradJthick_face[3],
      const double W_face,
      const double V_face[3],
      const double chi_tr_face,
      const double D_face,
      const double delta_l,
      double *restrict corrected_flux_tildeE,
      double *restrict a_face) {
  if(m1_params == NULL || metric_face == NULL || gradJthick_face == NULL ||
     V_face == NULL || corrected_flux_tildeE == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(hll_flux_tildeE) || !isfinite(E_star) ||
     !isfinite(chi_tr_face) || chi_tr_face < 0.0 ||
     !isfinite(D_face) || D_face < 0.0 ||
     !isfinite(delta_l) || delta_l <= 0.0)
    return ghl_error_m1_invalid_state;

  *corrected_flux_tildeE = hll_flux_tildeE;
  if(a_face != NULL)
    *a_face = 1.0;

  ghl_error_codes_t error = ghl_m1_validate_direction(direction);
  if(error != ghl_success)
    return error;
  if(!ghl_m1_metric_is_symmetric_spd(metric_face))
    return ghl_error_m1_invalid_metric;

  const double tau_face = chi_tr_face * delta_l;
  if(!isfinite(tau_face) || tau_face < 0.0)
    return ghl_error_m1_invalid_state;
  if(tau_face < m1_params->zeta_min)
    return ghl_success;
  if(!(Jthick_L_is_valid && isfinite(Jthick_L) && Jthick_L > 0.0) ||
     !(Jthick_R_is_valid && isfinite(Jthick_R) && Jthick_R > 0.0))
    return ghl_success;
  if(!ghl_m1_face_velocity_is_valid(metric_face, W_face, V_face))
    return ghl_success;
  for(int i = 0; i < 3; ++i) {
    if(!isfinite(gradJthick_face[i]))
      return ghl_error_m1_invalid_state;
  }

  const double a = tanh(1.0 / tau_face);
  if(!isfinite(a) || a <= 0.0 || a > 1.0)
    return ghl_error_m1_invalid_state;
  const int d = (int)direction;
  const double Jthick_face = 0.5 * (Jthick_L + Jthick_R);
  if(!isfinite(Jthick_face) || Jthick_face <= 0.0)
    return ghl_error_m1_invalid_state;

  const double W2 = SQR(W_face);
  const double F_adv = (4.0 / 3.0) * W2 * V_face[d] * Jthick_face;
  double gamma_grad = 0.0;
  double VdotGrad = 0.0;
  for(int i = 0; i < 3; ++i) {
    gamma_grad += metric_face->gammaUU[d][i] * gradJthick_face[i];
    VdotGrad += V_face[i] * gradJthick_face[i];
  }
  const double F_diff = W_face * D_face *
      (gamma_grad + V_face[d] * VdotGrad);
  const double F_asym = F_adv - F_diff;
  if(!isfinite(F_asym))
    return ghl_error_m1_invalid_state;

  const double scalar_shift = fmax(E_star, m1_params->E_floor);
  const double flux_asym = metric_face->sqrt_detgamma *
      (metric_face->lapse * F_asym - metric_face->betaU[d] * scalar_shift);
  if(!isfinite(flux_asym))
    return ghl_error_m1_invalid_state;
  *corrected_flux_tildeE = a * hll_flux_tildeE + (1.0 - a) * flux_asym;
  if(!isfinite(*corrected_flux_tildeE))
    return ghl_error_m1_invalid_state;
  if(a_face != NULL)
    *a_face = a;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_neutrino_diffusion_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const double hll_flux_tildeE,
      const double E_star,
      const double Jthick_L,
      const bool Jthick_L_is_valid,
      const double Jthick_R,
      const bool Jthick_R_is_valid,
      const double gradJthick_face[3],
      const double W_face,
      const double V_face[3],
      const ghl_m1_neutrino_rates *restrict rates_face,
      const double D_face,
      const double delta_l,
      double *restrict corrected_flux_tildeE,
      double *restrict a_face) {
  if(rates_face == NULL)
    return ghl_error_m1_null_pointer;
  return ghl_m1_compute_diffusion_flux(
      m1_params, metric_face, direction, hll_flux_tildeE, E_star,
      Jthick_L, Jthick_L_is_valid, Jthick_R, Jthick_R_is_valid,
      gradJthick_face, W_face, V_face, rates_face->kappa_tr, D_face,
      delta_l, corrected_flux_tildeE, a_face);
}
