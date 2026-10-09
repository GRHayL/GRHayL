#include "ghl_m1_utils.h"
#include <float.h>

ghl_error_codes_t ghl_m1_compute_face_normal_delta_l(
      const ghl_metric_quantities *restrict metric_face,
      const ghl_m1_direction_t direction,
      const double delta_x_d,
      double *restrict delta_l) {

  if(metric_face == NULL || delta_l == NULL) {
    return ghl_error_m1_null_pointer;
  }

  ghl_error_codes_t error = ghl_m1_validate_direction(direction);
  if(error != ghl_success) {
    return error;
  }

  if(!isfinite(delta_x_d) || delta_x_d <= 0.0) {
    return ghl_error_m1_invalid_state;
  }

  const double gammaUU_dd = metric_face->gammaUU[(int)direction][(int)direction];
  *delta_l = delta_x_d / sqrt(gammaUU_dd);
  if(!isfinite(*delta_l) || *delta_l <= 0.0) {
    return ghl_error_m1_invalid_state;
  }

  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_harmonic_diffusion_coefficient(
      const double chi_tr_L,
      const double chi_tr_R,
      double *restrict D_face) {
  if(D_face == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(chi_tr_L) || !isfinite(chi_tr_R) || chi_tr_L <= 0.0 || chi_tr_R <= 0.0) {
    return ghl_error_m1_invalid_state;
  }

  const double D_L = 1.0 / (3.0 * chi_tr_L);
  const double D_R = 1.0 / (3.0 * chi_tr_R);
  const double denom = D_L + D_R;
  double candidate = NAN;
  if(isfinite(D_L) && isfinite(D_R) && isfinite(denom) && denom > 0.0) {
    candidate = 2.0 * D_L * D_R / denom;
  }
  if(!isfinite(candidate) || candidate <= 0.0) {
    /* The harmonic mean simplifies to 2/[3*(chi_L+chi_R)]. Scale the
     * opacity sum so its intermediate cannot overflow. */
    const double scale = fmax(chi_tr_L, chi_tr_R);
    const double ratio = fmin(chi_tr_L, chi_tr_R) / scale;
    candidate = (2.0 / (3.0 * (1.0 + ratio))) / scale;
  }
  /* On the fallback path the numerator is in [1/3, 2/3] and scale is
   * finite and positive. Even division by DBL_MAX exceeds DBL_TRUE_MIN,
   * so the result cannot round to zero. The fast path already required a
   * positive result. Overflow for subnormal opacities is still rejected. */
  if(!isfinite(candidate)) {
    return ghl_error_m1_invalid_state;
  }
  *D_face = candidate;
  return ghl_success;
}

/* The PSD caller supplies a finite symmetric matrix normalized to unit
 * maximum entry. Keeping one sweep separate makes both the rotation and the
 * bounded iteration contract independently testable. */
static bool ghl_m1_jacobi_sweep(void *workspace) {
  double (*B)[3] = workspace;
  const double eigen_tolerance = 64.0 * DBL_EPSILON;
  double offdiag = 0.0;
  for(int p = 0; p < 3; ++p) {
    for(int q = p + 1; q < 3; ++q) {
      offdiag = fmax(offdiag, fabs(B[p][q]));
    }
  }
  if(offdiag <= eigen_tolerance) {
    return true;
  }
  for(int p = 0; p < 3; ++p) {
    for(int q = p + 1; q < 3; ++q) {
      if(fabs(B[p][q]) <= eigen_tolerance) {
        continue;
      }
      const double tau = (B[q][q] - B[p][p]) / (2.0 * B[p][q]);
      const double t = copysign(1.0 / (fabs(tau) + sqrt(1.0 + tau * tau)), tau);
      const double c = 1.0 / sqrt(1.0 + t * t);
      const double s = t * c;
      const double Bpp = B[p][p];
      const double Bqq = B[q][q];
      B[p][p] = Bpp - t * B[p][q];
      B[q][q] = Bqq + t * B[p][q];
      B[p][q] = B[q][p] = 0.0;
      for(int k = 0; k < 3; ++k) {
        if(k == p || k == q) {
          continue;
        }
        const double Bkp = B[k][p];
        const double Bkq = B[k][q];
        B[k][p] = B[p][k] = c * Bkp - s * Bkq;
        B[k][q] = B[q][k] = s * Bkp + c * Bkq;
      }
    }
  }
  return false;
}

void ghl_m1_jacobi_iteration_driver(void *workspace, bool (*sweep)(void *)) {
  for(int iteration = 0; iteration < 32; ++iteration) {
    if(sweep(workspace)) {
      break;
    }
  }
}

void ghl_m1_jacobi_eigenvalues(double matrix[3][3]) {
  ghl_m1_jacobi_iteration_driver(matrix, ghl_m1_jacobi_sweep);
}

ghl_error_codes_t ghl_m1_finish_scaled_norm_ratio(
      const double x_scale,
      const double A_scale,
      const double denom,
      const double scaled_norm,
      double *restrict ratio) {
  if(!isfinite(scaled_norm) || scaled_norm <= 0.0) {
    return ghl_error_m1_invalid_state;
  }

  const double value = (x_scale / denom) * sqrt(A_scale) * scaled_norm;
  if(!isfinite(value)) {
    return ghl_error_m1_invalid_state;
  }
  *ratio = value;
  return ghl_success;
}
