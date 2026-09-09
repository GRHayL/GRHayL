#include "ghl_m1_utils.h"
#include <float.h>

#ifdef GRHAYL_M1_DEBUG
ghl_error_codes_t ghl_m1_validate_runtime_params(
      const ghl_m1_parameters *restrict m1_params) {

  if(m1_params == NULL)
    return ghl_error_m1_null_pointer;

  if(!isfinite(m1_params->epsilon_c) || m1_params->epsilon_c <= 0.0 || m1_params->epsilon_c >= 1.0)
    return ghl_error_m1_invalid_epsilon_c;
  if(!isfinite(m1_params->E_floor) || m1_params->E_floor <= 0.0)
    return ghl_error_m1_invalid_E_floor;
  if(m1_params->repair_policy != ghl_m1_repair_linear_factor_compatibility)
    return ghl_error_m1_invalid_repair_policy;
  if(!isfinite(m1_params->closure_root_tolerance) ||
     m1_params->closure_root_tolerance <= 0.0 ||
     m1_params->closure_root_tolerance > 1.0)
    return ghl_error_m1_invalid_closure_tolerance;
  if(m1_params->closure_root_max_iterations <= 0)
    return ghl_error_m1_invalid_closure_max_iterations;
  if(!isfinite(m1_params->closure_root_residual_tolerance) ||
     m1_params->closure_root_residual_tolerance <= 0.0)
    return ghl_error_m1_invalid_closure_tolerance;
  if(!isfinite(m1_params->one_minus_epsilon_c_sq) || m1_params->one_minus_epsilon_c_sq < 0.0)
    return ghl_error_m1_invalid_epsilon_c;
  if(!isfinite(m1_params->zeta_min) || m1_params->zeta_min <= 0.0)
    return ghl_error_m1_invalid_zeta_min;
  if(!isfinite(m1_params->fd_epsilon_rel) || m1_params->fd_epsilon_rel <= 0.0)
    return ghl_error_m1_invalid_fd_epsilon_rel;
  if(!isfinite(m1_params->fd_epsilon_abs) || m1_params->fd_epsilon_abs <= 0.0)
    return ghl_error_m1_invalid_fd_epsilon_abs;
  if(m1_params->newton_max_iterations <= 0)
    return ghl_error_m1_invalid_newton_max_iterations;
  if(!isfinite(m1_params->newton_tolerance) || m1_params->newton_tolerance <= 0.0)
    return ghl_error_m1_invalid_newton_tolerance;
  if(!isfinite(m1_params->newton_absolute_tolerance) ||
     m1_params->newton_absolute_tolerance <= 0.0)
    return ghl_error_m1_invalid_newton_absolute_tolerance;
  const double expected = 1.0 - m1_params->epsilon_c;
  const double scale = ghl_m1_max(1.0, fabs(expected));
  const double tol = 64.0 * DBL_EPSILON * scale;
  if(fabs(m1_params->one_minus_epsilon_c_sq - expected) > tol)
    return ghl_error_m1_invalid_epsilon_c;

  return ghl_success;
}
#endif

typedef struct {
  double mantissa;
  int exponent;
} ghl_m1_positive_scaled_value;

static bool ghl_m1_is_symmetric_spd_3x3(
      const double A[3][3],
      double L[3][3],
      double *restrict matrix_scale) {

  *matrix_scale = 0.0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(!isfinite(A[i][j]))
        return false;
      *matrix_scale = ghl_m1_max(*matrix_scale, fabs(A[i][j]));
    }
  }
  if(*matrix_scale == 0.0)
    return false;

  for(int i = 0; i < 3; ++i) {
    for(int j = i + 1; j < 3; ++j) {
      const double pair_scale = ghl_m1_max(fabs(A[i][j]), fabs(A[j][i]));
      if(fabs(A[i][j] - A[j][i]) >
         64.0 * DBL_EPSILON * pair_scale)
        return false;
    }
  }

  /* Cholesky on a uniformly scaled matrix avoids overflow in minors and the
   * determinant while retaining the SPD decision. */
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      L[i][j] = 0.0;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j <= i; ++j) {
      double value = A[i][j] / *matrix_scale;
      for(int k = 0; k < j; ++k)
        value -= L[i][k] * L[j][k];
      if(i == j) {
        if(!isfinite(value) || value <= 0.0)
          return false;
        L[i][j] = sqrt(value);
      } else {
        L[i][j] = value / L[j][j];
        if(!isfinite(L[i][j]))
          return false;
      }
    }
  }
  return true;
}

static ghl_m1_positive_scaled_value ghl_m1_scaled_determinant_from_cholesky(
      const double matrix_scale,
      const double L[3][3]) {

  int scale_exponent;
  const double scale_mantissa = frexp(matrix_scale, &scale_exponent);
  double diagonal_mantissa = 1.0;
  int diagonal_exponent = 0;
  for(int i = 0; i < 3; ++i) {
    int factor_exponent;
    diagonal_mantissa *= frexp(L[i][i], &factor_exponent);
    diagonal_exponent += factor_exponent;
  }
  int normalization_exponent;
  const double mantissa = frexp(
      scale_mantissa * scale_mantissa * scale_mantissa
    * diagonal_mantissa * diagonal_mantissa,
      &normalization_exponent);
  const ghl_m1_positive_scaled_value determinant = {
      .mantissa = mantissa,
      .exponent = 3 * scale_exponent + 2 * diagonal_exponent
                + normalization_exponent};
  return determinant;
}

static ghl_m1_positive_scaled_value ghl_m1_scaled_square(
      const double value) {

  int value_exponent, normalization_exponent;
  const double value_mantissa = frexp(value, &value_exponent);
  const double mantissa = frexp(
      value_mantissa * value_mantissa, &normalization_exponent);
  const ghl_m1_positive_scaled_value square = {
      .mantissa = mantissa,
      .exponent = 2 * value_exponent + normalization_exponent};
  return square;
}

static bool ghl_m1_scaled_value_agrees_with_double(
      const ghl_m1_positive_scaled_value computed,
      const double stored) {

  int stored_exponent;
  const double stored_mantissa = frexp(stored, &stored_exponent);
  const int exponent_difference = computed.exponent - stored_exponent;
  if(exponent_difference < -1 || exponent_difference > 1)
    return false;
  const double computed_at_stored_exponent =
      scalbn(computed.mantissa, exponent_difference);
  const double scale = ghl_m1_max(
      fabs(computed_at_stored_exponent), fabs(stored_mantissa));
  return fabs(computed_at_stored_exponent - stored_mantissa)
      <= 512.0 * DBL_EPSILON * scale;
}

bool ghl_m1_metric_is_symmetric_spd(
      const ghl_metric_quantities *restrict metric) {

  if(!isfinite(metric->lapse) || metric->lapse <= 0.0)
    return false;
  if(!isfinite(metric->detgamma) || metric->detgamma <= 0.0)
    return false;
  if(!isfinite(metric->sqrt_detgamma) || metric->sqrt_detgamma <= 0.0)
    return false;

  for(int i = 0; i < 3; i++) {
    if(!isfinite(metric->betaU[i]))
      return false;
  }

  double gammaDD_cholesky[3][3], gammaUU_cholesky[3][3];
  double gammaDD_scale, gammaUU_scale;
  if(!ghl_m1_is_symmetric_spd_3x3(
         metric->gammaDD, gammaDD_cholesky, &gammaDD_scale))
    return false;
  if(!ghl_m1_is_symmetric_spd_3x3(
         metric->gammaUU, gammaUU_cholesky, &gammaUU_scale))
    return false;

  const ghl_m1_positive_scaled_value computed_detgamma =
      ghl_m1_scaled_determinant_from_cholesky(
          gammaDD_scale, gammaDD_cholesky);
  if(!ghl_m1_scaled_value_agrees_with_double(
         computed_detgamma, metric->detgamma))
    return false;
  const ghl_m1_positive_scaled_value stored_sqrt_square =
      ghl_m1_scaled_square(metric->sqrt_detgamma);
  if(!ghl_m1_scaled_value_agrees_with_double(
         stored_sqrt_square, metric->detgamma))
    return false;

  long double max_abs_product_entry = 1.0L;
  long double max_abs_error = 0.0L;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      long double product_ij = 0.0L;
      for(int k = 0; k < 3; k++)
        product_ij += (long double)metric->gammaDD[i][k] *
                      metric->gammaUU[k][j];

      max_abs_product_entry = fmaxl(max_abs_product_entry, fabsl(product_ij));
      const long double target = (i == j) ? 1.0L : 0.0L;
      max_abs_error = fmaxl(max_abs_error, fabsl(product_ij - target));
    }
  }

  const long double inverse_tol =
      256.0L * DBL_EPSILON * max_abs_product_entry;
  return max_abs_error <= inverse_tol;
}

ghl_error_codes_t ghl_m1_compute_face_normal_delta_l(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const double delta_x_d,
      double *restrict delta_l) {

  if(metric_face == NULL || delta_l == NULL)
    return ghl_error_m1_null_pointer;

  ghl_error_codes_t error = ghl_m1_validate_direction(direction);
  if(error != ghl_success)
    return error;

  if(!isfinite(delta_x_d) || delta_x_d <= 0.0)
    return ghl_error_m1_invalid_state;

  if(!ghl_m1_metric_is_symmetric_spd(metric_face))
    return ghl_error_m1_invalid_metric;

  const double gammaUU_dd = metric_face->gammaUU[(int)direction][(int)direction];
  if(!isfinite(gammaUU_dd) || gammaUU_dd <= 0.0)
    return ghl_error_m1_invalid_metric;

  *delta_l = delta_x_d / sqrt(gammaUU_dd);
  if(!isfinite(*delta_l) || *delta_l <= 0.0)
    return ghl_error_m1_invalid_state;

  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_harmonic_diffusion_coefficient(
      const double chi_tr_L,
      const double chi_tr_R,
      double *restrict D_face) {
  if(D_face == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(chi_tr_L) || !isfinite(chi_tr_R) ||
     chi_tr_L <= 0.0 || chi_tr_R <= 0.0)
    return ghl_error_m1_invalid_state;

  const double D_L = 1.0 / (3.0 * chi_tr_L);
  const double D_R = 1.0 / (3.0 * chi_tr_R);
  const double denom = D_L + D_R;
  if(!isfinite(D_L) || !isfinite(D_R) || !isfinite(denom) || denom <= 0.0)
    return ghl_error_m1_invalid_state;

  *D_face = 2.0 * D_L * D_R / denom;
  if(!isfinite(*D_face) || *D_face <= 0.0)
    return ghl_error_m1_invalid_state;
  return ghl_success;
}
