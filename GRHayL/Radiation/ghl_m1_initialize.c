#include "ghl_m1.h"
#include "ghl_m1_utils.h"
#include <float.h>

ghl_error_codes_t ghl_m1_initialize(
      const double epsilon_c,
      const double E_floor,
      const double zeta_min,
      const double fd_epsilon_rel,
      const double fd_epsilon_abs,
      const int newton_max_iterations,
      const double newton_tolerance,
      ghl_m1_parameters *restrict m1_params) {

  const double legacy_absolute_tolerance =
      ghl_m1_max(E_floor * newton_tolerance, DBL_MIN);
  return ghl_m1_initialize_with_newton_tolerances(
      epsilon_c, E_floor, zeta_min, fd_epsilon_rel, fd_epsilon_abs,
      newton_max_iterations, newton_tolerance, legacy_absolute_tolerance,
      m1_params);
}

ghl_error_codes_t ghl_m1_set_newton_tolerances(
      const double newton_relative_tolerance,
      const double newton_absolute_tolerance,
      ghl_m1_parameters *restrict m1_params) {

  if(m1_params == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(newton_relative_tolerance) ||
     newton_relative_tolerance <= 0.0)
    return ghl_error_m1_invalid_newton_tolerance;
  if(!isfinite(newton_absolute_tolerance) ||
     newton_absolute_tolerance <= 0.0)
    return ghl_error_m1_invalid_newton_absolute_tolerance;

  m1_params->newton_tolerance = newton_relative_tolerance;
  m1_params->newton_absolute_tolerance = newton_absolute_tolerance;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_set_closure_solver_controls(
      const double root_interval_tolerance,
      const int root_max_iterations,
      ghl_m1_parameters *restrict m1_params) {
  if(m1_params == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(root_interval_tolerance) ||
     root_interval_tolerance <= 0.0 || root_interval_tolerance > 1.0)
    return ghl_error_m1_invalid_closure_tolerance;
  if(root_max_iterations <= 0)
    return ghl_error_m1_invalid_closure_max_iterations;
  m1_params->closure_root_tolerance = root_interval_tolerance;
  m1_params->closure_root_max_iterations = root_max_iterations;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_set_closure_residual_tolerance(
      const double max_normalized_residual,
      ghl_m1_parameters *restrict m1_params) {
  if(m1_params == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(max_normalized_residual) || max_normalized_residual <= 0.0)
    return ghl_error_m1_invalid_closure_tolerance;
  m1_params->closure_root_residual_tolerance = max_normalized_residual;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_initialize_with_newton_tolerances(
      const double epsilon_c,
      const double E_floor,
      const double zeta_min,
      const double fd_epsilon_rel,
      const double fd_epsilon_abs,
      const int newton_max_iterations,
      const double newton_relative_tolerance,
      const double newton_absolute_tolerance,
      ghl_m1_parameters *restrict m1_params) {

  if(m1_params == NULL)
    return ghl_error_m1_null_pointer;

  if(!isfinite(epsilon_c) || epsilon_c <= 0.0 || epsilon_c >= 1.0)
    return ghl_error_m1_invalid_epsilon_c;
  if(!isfinite(E_floor) || E_floor <= 0.0)
    return ghl_error_m1_invalid_E_floor;
  if(!isfinite(zeta_min) || zeta_min <= 0.0)
    return ghl_error_m1_invalid_zeta_min;
  if(!isfinite(fd_epsilon_rel) || fd_epsilon_rel <= 0.0)
    return ghl_error_m1_invalid_fd_epsilon_rel;
  if(!isfinite(fd_epsilon_abs) || fd_epsilon_abs <= 0.0)
    return ghl_error_m1_invalid_fd_epsilon_abs;
  if(newton_max_iterations <= 0)
    return ghl_error_m1_invalid_newton_max_iterations;
  const ghl_error_codes_t tolerance_error = ghl_m1_set_newton_tolerances(
      newton_relative_tolerance, newton_absolute_tolerance, m1_params);
  if(tolerance_error != ghl_success)
    return tolerance_error;

  m1_params->epsilon_c = epsilon_c;
  /* Historical field name: this stores the admissible squared reduced-flux
   * limit r = (F/E)^2 = 1 - epsilon_c. */
  m1_params->one_minus_epsilon_c_sq = 1.0 - epsilon_c;
  m1_params->E_floor = E_floor;
  m1_params->repair_policy = ghl_m1_repair_linear_factor_compatibility;
  m1_params->closure_root_tolerance = 1.0e-12;
  m1_params->closure_root_max_iterations = 100;
  m1_params->closure_root_residual_tolerance = 1.0e-10;
  m1_params->zeta_min = zeta_min;
  m1_params->fd_epsilon_rel = fd_epsilon_rel;
  m1_params->fd_epsilon_abs = fd_epsilon_abs;
  m1_params->newton_max_iterations = newton_max_iterations;
  /* Canonical four-point transport controls. */
  m1_params->minmod_theta = 1.0;
  m1_params->mindiss = 0.0;

#ifdef GRHAYL_M1_DEBUG
  ghl_error_codes_t debug_error = ghl_m1_validate_runtime_params(m1_params);
  if(debug_error != ghl_success)
    return debug_error;
#endif

  return ghl_success;
}
