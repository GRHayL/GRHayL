#include "ghl_m1.h"
#include "ghl_m1_utils.h"

ghl_error_codes_t ghl_m1_compute_raw_lightcone_speeds(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      double *restrict s_minus_raw,
      double *restrict s_plus_raw) {
  if(metric_face == NULL || s_minus_raw == NULL || s_plus_raw == NULL)
    return ghl_error_m1_null_pointer;
  ghl_error_codes_t error = ghl_m1_validate_direction(direction);
  if(error != ghl_success)
    return error;
  if(!ghl_m1_metric_is_symmetric_spd(metric_face))
    return ghl_error_m1_invalid_metric;
  const double scale = metric_face->lapse
                     * sqrt(metric_face->gammaUU[direction][direction]);
  if(!isfinite(scale) || scale <= 0.0)
    return ghl_error_m1_invalid_metric;
  const double minus = -metric_face->betaU[direction] - scale;
  const double plus = -metric_face->betaU[direction] + scale;
  if(!isfinite(minus) || !isfinite(plus))
    return ghl_error_m1_invalid_metric;
  *s_minus_raw = minus;
  *s_plus_raw = plus;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_clip_hll_speeds(
      const double s_minus_raw,
      const double s_plus_raw,
      double *restrict s_minus,
      double *restrict s_plus) {
  if(s_minus == NULL || s_plus == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(s_minus_raw) || !isfinite(s_plus_raw) ||
     s_minus_raw > s_plus_raw)
    return ghl_error_m1_invalid_state;
  *s_minus = ghl_m1_min(0.0, s_minus_raw);
  *s_plus = ghl_m1_max(0.0, s_plus_raw);
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_wavespeeds(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      double *restrict s_minus,
      double *restrict s_plus) {
  if(metric_face == NULL || s_minus == NULL || s_plus == NULL)
    return ghl_error_m1_null_pointer;
  double minus_raw, plus_raw;
  ghl_error_codes_t error = ghl_m1_compute_raw_lightcone_speeds(
      metric_face, direction, &minus_raw, &plus_raw);
  if(error != ghl_success)
    return error;
  return ghl_m1_clip_hll_speeds(minus_raw, plus_raw, s_minus, s_plus);
}
