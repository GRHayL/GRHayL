#include "ghl_m1.h"
#include "ghl_m1_utils.h"
#include <float.h>

static double ghl_m1_newton_residual_max_norm(const double residual[4]) {

  double norm = 0.0;
  for(int i = 0; i < 4; i++) {
    norm = ghl_m1_max(norm, fabs(residual[i]));
  }
  return norm;
}

double ghl_m1_newton_weighted_merit(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const double U[4],
      const double U_base[4],
      const double residual[4]) {

  const double energy_floor_tilde = metric->sqrt_detgamma * m1_params->E_floor;
  const double absolute_tolerance_tilde
        = metric->sqrt_detgamma * m1_params->newton_absolute_tolerance;
  const double energy_scale
        = ghl_m1_max(ghl_m1_max(fabs(U[0]), fabs(U_base[0])), energy_floor_tilde);

  double merit = 0.0;
  for(int i = 0; i < 4; i++) {
    double component_scale = ghl_m1_max(fabs(U[i]), fabs(U_base[i]));
    if(i > 0) {
      component_scale = ghl_m1_max(component_scale, energy_scale);
    }
    else {
      component_scale = ghl_m1_max(component_scale, energy_floor_tilde);
    }

    const double denominator
          = absolute_tolerance_tilde + m1_params->newton_tolerance * component_scale;
    if(!isfinite(denominator) || denominator <= 0.0 || !isfinite(residual[i])) {
      return INFINITY;
    }
    merit = ghl_m1_max(merit, fabs(residual[i]) / denominator);
  }
  return merit;
}

ghl_error_codes_t ghl_m1_newton_project_admissible(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      double U[4]) {

  if(!isfinite(metric->sqrt_detgamma) || metric->sqrt_detgamma <= 0.0) {
    return ghl_error_m1_invalid_metric;
  }

  const double inv_sqrt_detgamma = 1.0 / metric->sqrt_detgamma;
  ghl_m1_rad_state state = { .E = U[0] * inv_sqrt_detgamma };
  for(int i = 0; i < 3; i++) {
    state.F[i] = U[i + 1] * inv_sqrt_detgamma;
  }

  if(!isfinite(state.E) || !isfinite(state.F[0]) || !isfinite(state.F[1])
     || !isfinite(state.F[2])) {
    return ghl_error_m1_implicit_admissibility;
  }

  const ghl_error_codes_t error = ghl_m1_realizability_repair(m1_params, metric, &state);
  if(error != ghl_success) {
    return error;
  }

  U[0] = state.E * metric->sqrt_detgamma;
  for(int i = 0; i < 3; i++) {
    U[i + 1] = state.F[i] * metric->sqrt_detgamma;
  }

  if(!isfinite(U[0]) || !isfinite(U[1]) || !isfinite(U[2]) || !isfinite(U[3])) {
    return ghl_error_m1_implicit_admissibility;
  }
  return ghl_success;
}

static bool ghl_m1_newton_solve_linear_4x4(
      const double matrix_in[4][4],
      const double rhs_in[4],
      double solution[4]) {

  double A[4][4];
  double rhs[4];
  double row_scale[4];
  for(int i = 0; i < 4; i++) {
    row_scale[i] = 0.0;
    rhs[i] = rhs_in[i];
    if(!isfinite(rhs[i])) {
      return false;
    }
    for(int j = 0; j < 4; j++) {
      A[i][j] = matrix_in[i][j];
      if(!isfinite(A[i][j])) {
        return false;
      }
      row_scale[i] = ghl_m1_max(row_scale[i], fabs(A[i][j]));
    }
    if(row_scale[i] == 0.0) {
      return false;
    }
  }

  for(int col = 0; col < 4; col++) {
    int pivot_row = col;
    double best_scaled_pivot = -1.0;
    for(int row = col; row < 4; row++) {
      const double scaled_pivot = fabs(A[row][col]) / row_scale[row];
      if(scaled_pivot > best_scaled_pivot) {
        best_scaled_pivot = scaled_pivot;
        pivot_row = row;
      }
    }

    const double pivot_abs = fabs(A[pivot_row][col]);
    if(!isfinite(best_scaled_pivot) || best_scaled_pivot <= 64.0 * DBL_EPSILON
       || !isfinite(pivot_abs) || pivot_abs < DBL_MIN) {
      return false;
    }

    if(pivot_row != col) {
      for(int j = 0; j < 4; j++) {
        const double tmp = A[col][j];
        A[col][j] = A[pivot_row][j];
        A[pivot_row][j] = tmp;
      }
      const double rhs_tmp = rhs[col];
      rhs[col] = rhs[pivot_row];
      rhs[pivot_row] = rhs_tmp;
      const double scale_tmp = row_scale[col];
      row_scale[col] = row_scale[pivot_row];
      row_scale[pivot_row] = scale_tmp;
    }

    for(int row = col + 1; row < 4; row++) {
      const double factor = A[row][col] / A[col][col];
      if(!isfinite(factor)) {
        return false;
      }
      A[row][col] = 0.0;
      for(int j = col + 1; j < 4; j++) {
        A[row][j] -= factor * A[col][j];
      }
      rhs[row] -= factor * rhs[col];
    }
  }

  for(int i = 3; i >= 0; i--) {
    double value = rhs[i];
    for(int j = i + 1; j < 4; j++) {
      value -= A[i][j] * solution[j];
    }
    if(!isfinite(value) || !isfinite(A[i][i]) || fabs(A[i][i]) < DBL_MIN) {
      return false;
    }
    solution[i] = value / A[i][i];
    if(!isfinite(solution[i])) {
      return false;
    }
  }
  return true;
}

static void ghl_m1_newton_set_diagnostics(
      const int iterations,
      const int backtracks,
      const bool used_projection,
      const double residual_norm,
      const double merit,
      ghl_m1_newton_diagnostics *restrict diagnostics) {

  if(diagnostics == NULL) {
    return;
  }
  diagnostics->iterations = iterations;
  diagnostics->backtracks = backtracks;
  diagnostics->used_projection = used_projection;
  diagnostics->residual_max_norm = residual_norm;
  diagnostics->residual_weighted_merit = merit;
}

static void ghl_m1_newton_notify(
      const ghl_m1_newton_callbacks *restrict callbacks,
      const ghl_m1_solver_stage_t stage,
      const ghl_error_codes_t status,
      const double U[4],
      const double residual[4],
      const ghl_m1_newton_diagnostics *restrict diagnostics) {

  if(callbacks != NULL && callbacks->observer != NULL) {
    callbacks->observer(
          callbacks->observer_context, stage, status, U, residual, diagnostics);
  }
}

ghl_error_codes_t ghl_m1_newton_solve_4d_with_initial_guess(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_newton_callbacks *restrict callbacks,
      const void *restrict context,
      const double U_base[4],
      const double U_initial[4],
      double U_out[4],
      ghl_m1_newton_diagnostics *restrict diagnostics) {

  if(m1_params == NULL || metric == NULL || callbacks == NULL
     || callbacks->residual == NULL || callbacks->jacobian == NULL || U_base == NULL
     || U_initial == NULL || U_out == NULL) {
    return ghl_error_m1_null_pointer;
  }

  ghl_m1_newton_set_diagnostics(0, 0, false, INFINITY, INFINITY, diagnostics);
  double U[4] = { U_initial[0], U_initial[1], U_initial[2], U_initial[3] };
  double residual[4] = { 0.0 };
  bool residual_is_current = false;
  int backtracks = 0;
  bool used_projection = false;

  for(int iter = 0; iter < m1_params->newton_max_iterations; iter++) {
    ghl_error_codes_t error = ghl_success;
    if(!residual_is_current) {
      error = callbacks->residual(context, U, residual);
    }
    ghl_m1_newton_notify(
          callbacks, ghl_m1_solver_stage_residual, error, U, residual, diagnostics);
    if(error != ghl_success) {
      ghl_m1_newton_notify(
            callbacks, ghl_m1_solver_stage_completed_solve, error, U, residual,
            diagnostics);
      return error;
    }

    const double residual_norm = ghl_m1_newton_residual_max_norm(residual);
    const double merit
          = ghl_m1_newton_weighted_merit(m1_params, metric, U, U_base, residual);
    ghl_m1_newton_set_diagnostics(
          iter, backtracks, used_projection, residual_norm, merit, diagnostics);

    if(merit <= 1.0) {
      for(int i = 0; i < 4; i++) {
        U_out[i] = U[i];
      }
      ghl_m1_newton_notify(
            callbacks, ghl_m1_solver_stage_completed_solve, ghl_success, U_out, residual,
            diagnostics);
      return ghl_success;
    }

    double jacobian[4][4] = { { 0.0 } };
    error = callbacks->jacobian(context, U, residual, jacobian);
    ghl_m1_newton_notify(
          callbacks, ghl_m1_solver_stage_jacobian, error, U, residual, diagnostics);
    if(error != ghl_success) {
      ghl_m1_newton_notify(
            callbacks, ghl_m1_solver_stage_completed_solve, error, U, residual,
            diagnostics);
      return error;
    }

    double rhs[4];
    for(int i = 0; i < 4; i++) {
      rhs[i] = -residual[i];
    }
    double delta[4] = { 0.0 };
    if(!ghl_m1_newton_solve_linear_4x4(jacobian, rhs, delta)) {
      ghl_m1_newton_notify(
            callbacks, ghl_m1_solver_stage_newton_step,
            ghl_error_m1_implicit_solve_failure, U, residual, diagnostics);
      ghl_m1_newton_notify(
            callbacks, ghl_m1_solver_stage_completed_solve,
            ghl_error_m1_implicit_solve_failure, U, residual, diagnostics);
      return ghl_error_m1_implicit_solve_failure;
    }

    double step_scale = 1.0;
    bool accepted = false;
    while(step_scale >= 1.0 / 1024.0) {
      double trial_U[4];
      for(int i = 0; i < 4; i++) {
        trial_U[i] = U[i] + step_scale * delta[i];
      }

      double trial_residual[4] = { 0.0 };
      error = callbacks->residual(context, trial_U, trial_residual);
      ghl_m1_newton_notify(
            callbacks, ghl_m1_solver_stage_residual, error, trial_U, trial_residual,
            diagnostics);
      if(error == ghl_error_m1_implicit_admissibility) {
        const ghl_error_codes_t projection_error
              = ghl_m1_newton_project_admissible(m1_params, metric, trial_U);
        if(projection_error == ghl_success) {
          used_projection = true;
          error = callbacks->residual(context, trial_U, trial_residual);
          ghl_m1_newton_notify(
                callbacks, ghl_m1_solver_stage_residual, error, trial_U, trial_residual,
                diagnostics);
          if(error != ghl_success && error != ghl_error_m1_implicit_admissibility) {
            ghl_m1_newton_notify(
                  callbacks, ghl_m1_solver_stage_completed_solve, error, U, residual,
                  diagnostics);
            return error;
          }
        }
        else if(projection_error != ghl_error_m1_implicit_admissibility) {
          ghl_m1_newton_notify(
                callbacks, ghl_m1_solver_stage_completed_solve, projection_error, U,
                residual, diagnostics);
          return projection_error;
        }
      }
      else if(error != ghl_success) {
        ghl_m1_newton_notify(
              callbacks, ghl_m1_solver_stage_completed_solve, error, U, residual,
              diagnostics);
        return error;
      }

      if(error == ghl_success) {
        const double trial_norm = ghl_m1_newton_residual_max_norm(trial_residual);
        const double trial_merit = ghl_m1_newton_weighted_merit(
              m1_params, metric, trial_U, U_base, trial_residual);
        if(trial_merit < merit || trial_merit <= 1.0) {
          for(int i = 0; i < 4; i++) {
            U[i] = trial_U[i];
            residual[i] = trial_residual[i];
          }
          residual_is_current = true;
          ghl_m1_newton_set_diagnostics(
                iter + 1, backtracks, used_projection, trial_norm, trial_merit,
                diagnostics);
          accepted = true;
          if(trial_merit <= 1.0) {
            for(int i = 0; i < 4; i++) {
              U_out[i] = U[i];
            }
            ghl_m1_newton_notify(
                  callbacks, ghl_m1_solver_stage_newton_step, ghl_success, U_out,
                  trial_residual, diagnostics);
            ghl_m1_newton_notify(
                  callbacks, ghl_m1_solver_stage_completed_solve, ghl_success, U_out,
                  trial_residual, diagnostics);
            return ghl_success;
          }
          ghl_m1_newton_notify(
                callbacks, ghl_m1_solver_stage_newton_step, ghl_success, U,
                trial_residual, diagnostics);
          break;
        }
      }

      step_scale *= 0.5;
      backtracks++;
    }

    if(!accepted) {
      ghl_m1_newton_notify(
            callbacks, ghl_m1_solver_stage_newton_step,
            ghl_error_m1_implicit_solve_failure, U, residual, diagnostics);
      ghl_m1_newton_notify(
            callbacks, ghl_m1_solver_stage_completed_solve,
            ghl_error_m1_implicit_solve_failure, U, residual, diagnostics);
      return ghl_error_m1_implicit_solve_failure;
    }
  }

  ghl_m1_newton_notify(
        callbacks, ghl_m1_solver_stage_completed_solve,
        ghl_error_m1_implicit_solve_failure, U, NULL, diagnostics);
  return ghl_error_m1_implicit_solve_failure;
}

ghl_error_codes_t ghl_m1_newton_solve_4d(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_newton_callbacks *restrict callbacks,
      const void *restrict context,
      const double U_base[4],
      double U_out[4],
      ghl_m1_newton_diagnostics *restrict diagnostics) {
  return ghl_m1_newton_solve_4d_with_initial_guess(
        m1_params, metric, callbacks, context, U_base, U_base, U_out, diagnostics);
}
