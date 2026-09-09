#ifndef GHL_M1_UTILS_H
#define GHL_M1_UTILS_H

#include "ghl_m1.h"
#include <float.h>

/* Radiation-private Newton entry point with a separate initial iterate. */
ghl_error_codes_t ghl_m1_newton_solve_4d_with_initial_guess(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_newton_callbacks *restrict callbacks,
      const void *restrict context,
      const double U_base[4],
      const double U_initial[4],
      double U_out[4],
      ghl_m1_newton_diagnostics *restrict diagnostics);

/*
 * Keep M1's historical ternary semantics private to the radiation
 * implementation.  In particular, the false branch preserves the second
 * operand for equal values, signed zero, and NaN inputs.
 */
static inline double ghl_m1_min(const double A, const double B) {
  return A < B ? A : B;
}

static inline double ghl_m1_max(const double A, const double B) {
  return A > B ? A : B;
}

bool ghl_m1_metric_is_symmetric_spd(
      const ghl_metric_quantities *restrict metric);

/* Private cross-translation-unit observability hook. */
void ghl_m1_record_closure_downstream_repair(void);
void ghl_m1_record_closure_validation_failure(const int reason);

/* Radiation-private variant that returns the Eulerian velocity quantities
 * already computed while transforming the radiation moments. */
ghl_error_codes_t ghl_m1_compute_comoving_moments_with_velocity(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving,
      double V_con[3],
      double V_cov[3],
      double *restrict W);

/* Same transformation for a caller that has already validated the immutable
 * M1 configuration and metric for the current operation. */
ghl_error_codes_t ghl_m1_compute_comoving_moments_validated(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving,
      double V_con[3],
      double V_cov[3],
      double *restrict W);

static inline void ghl_m1_lower_spatial_tensor(
      const ghl_metric_quantities *restrict metric,
      const double tensor_UU[3][3],
      double tensor_DD[3][3]) {

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      double sum = 0.0;
      for(int k = 0; k < 3; k++) {
        for(int l = 0; l < 3; l++)
          sum += metric->gammaDD[i][k] * metric->gammaDD[j][l]
               * tensor_UU[k][l];
      }
      tensor_DD[i][j] = sum;
    }
  }
}

enum {
  GHL_M1_CLOSURE_VALIDATION_NONFINITE = 1,
  GHL_M1_CLOSURE_VALIDATION_TRACE = 2,
  GHL_M1_CLOSURE_VALIDATION_SYMMETRY = 3,
  GHL_M1_CLOSURE_VALIDATION_PSD = 4
};

static inline double ghl_m1_compute_fd_delta(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const double U_base,
      const double U) {

  const double E_floor_scale = metric->sqrt_detgamma * m1_params->E_floor;
  const double component_scale = ghl_m1_max(
      ghl_m1_max(fabs(U), fabs(U_base)), E_floor_scale);
  return m1_params->fd_epsilon_rel * component_scale
         + m1_params->fd_epsilon_abs * E_floor_scale;
}

static inline bool
ghl_m1_fd_error_allows_one_sided_fallback(const ghl_error_codes_t error) {

  return error == ghl_error_m1_implicit_admissibility;
}

static inline bool ghl_m1_schedule_error_allows_retry(
      const ghl_error_codes_t error) {

  return error == ghl_error_m1_implicit_admissibility
         || error == ghl_error_m1_invalid_implicit_jacobian
         || error == ghl_error_m1_implicit_solve_failure;
}

static inline void ghl_m1_initialize_implicit_solve_diagnostics(
      ghl_m1_implicit_solve_diagnostics *restrict diagnostics) {

  if(diagnostics == NULL) {
    return;
  }

  *diagnostics = (ghl_m1_implicit_solve_diagnostics){
        .fallback_substeps = 1,
        .residual_max_norm = INFINITY,
        .residual_scaled_norm = INFINITY };
}

static inline ghl_error_codes_t ghl_m1_validate_parameters(
      const ghl_m1_parameters *restrict m1_params);

static inline ghl_error_codes_t ghl_m1_validate_configuration(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {

  if(!ghl_m1_metric_is_symmetric_spd(metric))
    return ghl_error_m1_invalid_metric;
  return ghl_m1_validate_parameters(m1_params);
}

/* The implicit neutrino solve validates its immutable metric and M1
 * parameters once before Newton iterations. This helper retains the
 * state-dependent half of the public realizability check for each trial. */
static inline ghl_error_codes_t ghl_m1_validate_parameters(
      const ghl_m1_parameters *restrict m1_params) {

  if(!isfinite(m1_params->E_floor) || m1_params->E_floor <= 0.0)
    return ghl_error_m1_invalid_E_floor;
  if(m1_params->repair_policy != ghl_m1_repair_linear_factor_compatibility)
    return ghl_error_m1_invalid_repair_policy;
  if(!isfinite(m1_params->epsilon_c) || m1_params->epsilon_c <= 0.0 ||
     m1_params->epsilon_c >= 1.0 ||
     !isfinite(m1_params->one_minus_epsilon_c_sq))
    return ghl_error_m1_invalid_epsilon_c;
  if(!isfinite(m1_params->closure_root_tolerance) ||
     m1_params->closure_root_tolerance <= 0.0 ||
     m1_params->closure_root_tolerance > 1.0)
    return ghl_error_m1_invalid_closure_tolerance;
  if(m1_params->closure_root_max_iterations <= 0)
    return ghl_error_m1_invalid_closure_max_iterations;
  if(!isfinite(m1_params->closure_root_residual_tolerance) ||
     m1_params->closure_root_residual_tolerance <= 0.0)
    return ghl_error_m1_invalid_closure_tolerance;

  const double expected = 1.0 - m1_params->epsilon_c;
  const double scale = ghl_m1_max(fabs(expected),
                                  fabs(m1_params->one_minus_epsilon_c_sq));
  if(scale == 0.0 ||
     fabs(m1_params->one_minus_epsilon_c_sq - expected) >
       64.0 * DBL_EPSILON * scale)
    return ghl_error_m1_invalid_epsilon_c;

  return ghl_success;
}

/* Evaluate sqrt(x^T A x)/denom without materializing an overflow- or
 * underflow-prone quadratic form. A is assumed finite SPD and denom positive. */
static inline ghl_error_codes_t ghl_m1_scaled_norm_ratio(
      const double A[3][3],
      const double x[3],
      const double denom,
      double *restrict ratio) {

  if(!isfinite(denom) || denom <= 0.0 || ratio == NULL)
    return ghl_error_m1_invalid_state;

  double x_scale = 0.0;
  double A_scale = 0.0;
  for(int i = 0; i < 3; ++i) {
    if(!isfinite(x[i]))
      return ghl_error_m1_invalid_state;
    x_scale = ghl_m1_max(x_scale, fabs(x[i]));
    for(int j = 0; j < 3; ++j) {
      if(!isfinite(A[i][j]))
        return ghl_error_m1_invalid_metric;
      A_scale = ghl_m1_max(A_scale, fabs(A[i][j]));
    }
  }
  if(x_scale == 0.0) {
    *ratio = 0.0;
    return ghl_success;
  }
  if(A_scale == 0.0)
    return ghl_error_m1_invalid_metric;

  /* Factor the uniformly scaled SPD matrix and evaluate ||L^T x||_2.
   * Summing x^T A x directly is range-safe after scaling, but can still lose
   * its small positive result through cancellation for an ill-conditioned
   * metric.  Cholesky turns the contraction into a sum of squares. */
  double L[3][3] = {{0.0}};
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j <= i; ++j) {
      double value = A[i][j] / A_scale;
      for(int k = 0; k < j; ++k)
        value -= L[i][k] * L[j][k];
      if(i == j) {
        if(!isfinite(value) || value <= 0.0)
          return ghl_error_m1_invalid_metric;
        L[i][j] = sqrt(value);
      } else {
        L[i][j] = value / L[j][j];
        if(!isfinite(L[i][j]))
          return ghl_error_m1_invalid_metric;
      }
    }
  }

  double LT_x[3] = {0.0, 0.0, 0.0};
  for(int i = 0; i < 3; ++i)
    for(int j = i; j < 3; ++j)
      LT_x[i] += L[j][i] * (x[j] / x_scale);
  const double scaled_norm = hypot(hypot(LT_x[0], LT_x[1]), LT_x[2]);
  if(!isfinite(scaled_norm) || scaled_norm <= 0.0)
    return ghl_error_m1_invalid_state;

  int x_exp, A_exp, denom_exp;
  const double x_mant = frexp(x_scale, &x_exp);
  double A_mant = frexp(A_scale, &A_exp);
  const double denom_mant = frexp(denom, &denom_exp);
  if((A_exp & 1) != 0) {
    A_mant *= 2.0;
    --A_exp;
  }
  const double mant =
      (x_mant / denom_mant) * sqrt(A_mant) * scaled_norm;
  const int exponent = x_exp - denom_exp + A_exp / 2;
  const double value = scalbn(mant, exponent);
  *ratio = isfinite(value) ? value : DBL_MAX;
  return ghl_success;
}

static inline ghl_error_codes_t ghl_m1_scaled_covector_norm_ratio(
      const double gammaUU[3][3],
      const double F_cov[3],
      const double E,
      double *restrict flux_factor) {
  return ghl_m1_scaled_norm_ratio(gammaUU, F_cov, E, flux_factor);
}

static inline ghl_error_codes_t ghl_m1_scaled_vector_norm(
      const double gammaDD[3][3],
      const double V_con[3],
      double *restrict magnitude) {
  return ghl_m1_scaled_norm_ratio(gammaDD, V_con, 1.0, magnitude);
}

static inline ghl_error_codes_t ghl_m1_validate_direction(
      const ghl_m1_direction_t direction) {

  if(direction < ghl_m1_dirn0 || direction > ghl_m1_dirn2)
    return ghl_error_m1_invalid_state;

  return ghl_success;
}

static inline ghl_error_codes_t ghl_m1_compute_eulerian_velocity(
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      double V_con[3],
      double V_cov[3],
      double *restrict W_out) {

  const double alpha = metric->lapse;
  const double inv_alpha = 1.0 / alpha;

  for(int i = 0; i < 3; i++) {
    if(!isfinite(prims->vU[i]))
      return ghl_error_m1_invalid_state;
  }

  for(int i = 0; i < 3; i++)
    V_con[i] = (prims->vU[i] + metric->betaU[i]) * inv_alpha;

  if(V_cov != NULL)
    ghl_raise_lower_vector_3D(metric->gammaDD, V_con, V_cov);

  /* The Eulerian velocity is already available above.  Computing W from it
   * avoids cancellation in g00 + 2 g0i v^i + gij v^i v^j at Eulerian rest. */
  double V_mag;
  const ghl_error_codes_t norm_error =
      ghl_m1_scaled_vector_norm(metric->gammaDD, V_con, &V_mag);
  if(norm_error != ghl_success || V_mag >= 1.0)
    return ghl_error_u0_singular;
  const double W = 1.0 / sqrt((1.0 - V_mag) * (1.0 + V_mag));
  if(!isfinite(W))
    return ghl_error_u0_singular;

  *W_out = W;
  return ghl_success;
}

static inline ghl_error_codes_t ghl_m1_validate_realizability_state(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const double tol_factor,
      double *restrict flux_factor_sq_out);

static inline ghl_error_codes_t ghl_m1_validate_realizability(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const double tol_factor,
      double *restrict flux_factor_sq_out) {

  ghl_error_codes_t error = ghl_m1_validate_configuration(m1_params, metric);
  if(error != ghl_success)
    return error;
  return ghl_m1_validate_realizability_state(
      m1_params, metric, rad_state, tol_factor, flux_factor_sq_out);
}

static inline ghl_error_codes_t ghl_m1_validate_realizability_state(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const double tol_factor,
      double *restrict flux_factor_sq_out) {

  if(!isfinite(rad_state->E) || rad_state->E < m1_params->E_floor)
    return ghl_error_m1_invalid_state;
  for(int i = 0; i < 3; i++) {
    if(!isfinite(rad_state->F[i]))
      return ghl_error_m1_invalid_state;
  }

  double flux_factor;
  const ghl_error_codes_t error = ghl_m1_scaled_covector_norm_ratio(
      metric->gammaUU, rad_state->F, rad_state->E, &flux_factor);
  if(error != ghl_success)
    return error;
  const double cone_factor = 1.0 - m1_params->epsilon_c;
  const double permitted_flux_factor = sqrt(cone_factor);
  const double scale = ghl_m1_max(flux_factor, permitted_flux_factor);
  if(flux_factor > permitted_flux_factor &&
     flux_factor - permitted_flux_factor >
       tol_factor * DBL_EPSILON * scale)
    return ghl_error_m1_invalid_state;

  if(flux_factor_sq_out != NULL)
    *flux_factor_sq_out = flux_factor * flux_factor;

  return ghl_success;
}

/* Check the spatial pressure tensor in a metric-independent basis. The
 * Jacobi eigensolver is deliberately local and small: unlike a determinant
 * test, it remains reliable for rank-deficient realizable tensors. */
static inline ghl_error_codes_t ghl_m1_validate_closure_tensor_psd(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_closure *restrict closure) {

  double A[3][3] = {{0.0}};
  ghl_m1_lower_spatial_tensor(metric, closure->P, A);
  double scale = 0.0;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) {
      if(!isfinite(A[i][j])) {
        ghl_m1_record_closure_validation_failure(
            GHL_M1_CLOSURE_VALIDATION_NONFINITE);
        return ghl_error_m1_invalid_state;
      }
      scale = ghl_m1_max(scale, fabs(A[i][j]));
    }
  if(!isfinite(scale) || scale <= 0.0) {
    ghl_m1_record_closure_validation_failure(
        GHL_M1_CLOSURE_VALIDATION_PSD);
    return ghl_error_m1_invalid_state;
  }

  /* A congruence preserves positive semidefiniteness and the metric Cholesky
   * factor provides a numerically well-scaled representative. */
  double L[3][3] = {{0.0}};
  double gamma_scale = 0.0;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      gamma_scale = ghl_m1_max(gamma_scale, fabs(metric->gammaDD[i][j]));
  if(!isfinite(gamma_scale) || gamma_scale <= 0.0) {
    ghl_m1_record_closure_validation_failure(
        GHL_M1_CLOSURE_VALIDATION_PSD);
    return ghl_error_m1_invalid_state;
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j <= i; ++j) {
      double value = metric->gammaDD[i][j] / gamma_scale;
      for(int k = 0; k < j; ++k)
        value -= L[i][k] * L[j][k];
      if(i == j) {
        if(!isfinite(value) || value <= 0.0) {
          ghl_m1_record_closure_validation_failure(
              GHL_M1_CLOSURE_VALIDATION_PSD);
          return ghl_error_m1_invalid_state;
        }
        L[i][j] = sqrt(value);
      } else {
        L[i][j] = value / L[j][j];
        if(!isfinite(L[i][j])) {
          ghl_m1_record_closure_validation_failure(
              GHL_M1_CLOSURE_VALIDATION_PSD);
          return ghl_error_m1_invalid_state;
        }
      }
    }
  }
  double Linv[3][3] = {{0.0}};
  for(int col = 0; col < 3; ++col) {
    Linv[col][col] = 1.0 / L[col][col];
    for(int row = col + 1; row < 3; ++row) {
      double value = 0.0;
      for(int k = col; k < row; ++k)
        value -= L[row][k] * Linv[k][col];
      Linv[row][col] = value / L[row][row];
    }
  }
  double B[3][3] = {{0.0}};
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      for(int k = 0; k < 3; ++k)
        for(int l = 0; l < 3; ++l)
          B[i][j] += Linv[i][k] * A[k][l] * Linv[j][l];

  double eigen_scale = 0.0;
  for(int i = 0; i < 3; ++i) {
    for(int j = i + 1; j < 3; ++j) {
      const double symmetric = 0.5 * (B[i][j] + B[j][i]);
      B[i][j] = symmetric;
      B[j][i] = symmetric;
    }
    eigen_scale = ghl_m1_max(eigen_scale, fabs(B[i][i]));
  }
  for(int sweep = 0; sweep < 32; ++sweep) {
    double offdiag = 0.0;
    for(int p = 0; p < 3; ++p)
      for(int q = p + 1; q < 3; ++q)
        offdiag = ghl_m1_max(offdiag, fabs(B[p][q]));
    if(offdiag <= 64.0 * DBL_EPSILON * ghl_m1_max(eigen_scale, 1.0))
      break;
    for(int p = 0; p < 3; ++p) {
      for(int q = p + 1; q < 3; ++q) {
        if(fabs(B[p][q]) <= 64.0 * DBL_EPSILON *
                              ghl_m1_max(eigen_scale, 1.0))
          continue;
        const double tau = (B[q][q] - B[p][p]) / (2.0 * B[p][q]);
        const double t = copysign(1.0 / (fabs(tau) + sqrt(1.0 + tau * tau)),
                                  tau);
        const double c = 1.0 / sqrt(1.0 + t * t);
        const double s = t * c;
        const double Bpp = B[p][p];
        const double Bqq = B[q][q];
        B[p][p] = Bpp - t * B[p][q];
        B[q][q] = Bqq + t * B[p][q];
        B[p][q] = B[q][p] = 0.0;
        for(int k = 0; k < 3; ++k) {
          if(k == p || k == q)
            continue;
          const double Bkp = B[k][p];
          const double Bkq = B[k][q];
          B[k][p] = B[p][k] = c * Bkp - s * Bkq;
          B[k][q] = B[q][k] = s * Bkp + c * Bkq;
        }
      }
    }
  }
  for(int i = 0; i < 3; ++i) {
    if(!isfinite(B[i][i]) || B[i][i] <
       -1024.0 * DBL_EPSILON * ghl_m1_max(eigen_scale, 1.0)) {
      ghl_m1_record_closure_validation_failure(
          GHL_M1_CLOSURE_VALIDATION_PSD);
      return ghl_error_m1_invalid_state;
    }
  }
  return ghl_success;
}

static inline ghl_error_codes_t ghl_m1_validate_closure_tensor(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure) {

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(!isfinite(closure->P[i][j])) {
        ghl_m1_record_closure_validation_failure(
            GHL_M1_CLOSURE_VALIDATION_NONFINITE);
        return ghl_error_m1_invalid_state;
      }
    }
  }

  long double traceP = 0.0L;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) {
      traceP += (long double)metric->gammaDD[i][j] * closure->P[i][j];
    }
  const long double trace_scale = fmaxl(fabsl(traceP), fabsl(rad_state->E));
  if(fabsl(traceP - rad_state->E) >
     128.0L * DBL_EPSILON * trace_scale) {
    ghl_m1_record_closure_validation_failure(
        GHL_M1_CLOSURE_VALIDATION_TRACE);
    return ghl_error_m1_invalid_state;
  }

  for(int i = 0; i < 3; ++i) {
    for(int j = i + 1; j < 3; ++j) {
      const double pair_scale = ghl_m1_max(fabs(closure->P[i][j]),
                                           fabs(closure->P[j][i]));
      if(fabs(closure->P[i][j] - closure->P[j][i]) >
         64.0 * DBL_EPSILON * pair_scale) {
        ghl_m1_record_closure_validation_failure(
            GHL_M1_CLOSURE_VALIDATION_SYMMETRY);
        return ghl_error_m1_invalid_state;
      }
    }
  }

  return ghl_m1_validate_closure_tensor_psd(metric, closure);
}

static inline ghl_error_codes_t ghl_m1_validate_transport_velocity(
      const ghl_metric_quantities *restrict metric,
      const double velocity[3]) {

  double magnitude;
  const ghl_error_codes_t error = ghl_m1_scaled_vector_norm(
      metric->gammaDD, velocity, &magnitude);
  if(error != ghl_success)
    return error;
  if(magnitude > 1.0 + 128.0 * DBL_EPSILON)
    return ghl_error_m1_invalid_state;
  return ghl_success;
}

/* Private production decomposition used by closure diagnostics and unit
 * oracles. Outputs are contravariant spatial tensors and are transactional. */
ghl_error_codes_t ghl_m1_compute_minerbo_decomposition(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      double Pthin[3][3],
      double Pthick[3][3]);

/* Compute the configured closure after the caller has validated the
 * immutable M1 configuration and metric for the current operation. */
ghl_error_codes_t ghl_m1_compute_closure_minerbo_validated(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure);

static inline ghl_error_codes_t ghl_m1_validate_face_transport_scalars(
      const ghl_m1_direction_t direction,
      const double chi_tr_face,
      const double delta_l) {

#ifdef GRHAYL_M1_DEBUG
  ghl_error_codes_t error = ghl_m1_validate_direction(direction);
  if(error != ghl_success) {
    return error;
  }

  if(!isfinite(chi_tr_face) || chi_tr_face < 0.0) {
    return ghl_error_m1_invalid_state;
  }
  if(!isfinite(delta_l) || delta_l <= 0.0) {
    return ghl_error_m1_invalid_state;
  }
#else
  (void)direction;
  (void)chi_tr_face;
  (void)delta_l;
#endif

  return ghl_success;
}

/**
 * Validate M1 runtime parameters after initialization or debug-time mutation.
 *
 * This is a debug-only internal helper used by initialization and debug tests.
 */
#ifdef GRHAYL_M1_DEBUG
ghl_error_codes_t ghl_m1_validate_runtime_params(
      const ghl_m1_parameters *restrict m1_params);
#endif

#endif // GHL_M1_UTILS_H
