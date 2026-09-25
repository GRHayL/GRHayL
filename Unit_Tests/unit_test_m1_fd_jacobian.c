#include "../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit.h"
#include "../GRHayL/Radiation/ghl_m1_utils.h"
#include "ghl_m1.h"
#include "m1_test_utils.h"

#include <float.h>
#include <math.h>
#include <string.h>

/*
 * Regression coverage for the shared E/F finite-difference scale used by the
 * neutrino implicit Jacobian.  The tests deliberately keep the radiation
 * energy many orders above E_floor while one or more flux components are
 * exactly zero.  Every independently evaluated difference below uses the
 * coupled energy-scale step expected by the private Jacobian helper.
 */

static void require_condition(
      const bool condition,
      const char *restrict message,
      const int case_index) {
  if(!condition) {
    ghl_error("M1 FD Jacobian case %d: %s\n", case_index, message);
  }
}

static void require_error(
      const ghl_error_codes_t actual,
      const ghl_error_codes_t expected,
      const char *restrict operation,
      const int case_index) {
  if(actual != expected) {
    ghl_error(
          "M1 FD Jacobian case %d: %s returned %d, expected %d\n", case_index, operation,
          (int)actual, (int)expected);
  }
}

static bool close_value(
      const double actual,
      const double expected,
      const double relative_tolerance,
      const double absolute_tolerance) {
  return isfinite(actual) && isfinite(expected)
         && fabs(actual - expected)
                  <= absolute_tolerance
                           + relative_tolerance * fmax(fabs(actual), fabs(expected));
}

static void require_close(
      const double actual,
      const double expected,
      const double relative_tolerance,
      const double absolute_tolerance,
      const char *restrict quantity,
      const int case_index) {
  if(!close_value(actual, expected, relative_tolerance, absolute_tolerance)) {
    ghl_error(
          "M1 FD Jacobian case %d: %s mismatch (got %.17e, expected %.17e)\n",
          case_index, quantity, actual, expected);
  }
}

static void test_scaled_positive_argument_boundaries(void) {
  ghl_m1_scaled_positive value = { .mantissa = 0.25, .exponent = 2 };
  require_condition(
        !ghl_m1_scaled_positive_from_double(1.0, NULL),
        "scaled-positive conversion accepted a NULL output", 378);
  require_condition(
        !ghl_m1_scaled_positive_from_double(NAN, &value),
        "scaled-positive conversion accepted NaN", 378);
  require_condition(
        !ghl_m1_scaled_positive_from_double(INFINITY, &value),
        "scaled-positive conversion accepted infinity", 378);
  require_condition(
        !ghl_m1_scaled_positive_from_double(-1.0, &value),
        "scaled-positive conversion accepted a negative value", 378);
  require_condition(
        ghl_m1_scaled_positive_from_double(0.0, &value) && value.mantissa == 0.0
              && value.exponent == 0,
        "scaled-positive conversion mishandled zero", 378);

  const double factor = 2.0;
  require_condition(
        !ghl_m1_scaled_positive_product(NULL, 1, &value),
        "scaled-positive product accepted a NULL input", 379);
  require_condition(
        !ghl_m1_scaled_positive_product(&factor, 0, &value),
        "scaled-positive product accepted an empty input", 379);
  require_condition(
        !ghl_m1_scaled_positive_product(&factor, 1, NULL),
        "scaled-positive product accepted a NULL output", 379);
  require_condition(
        ghl_m1_scaled_positive_product(&factor, 1, &value)
              && value.mantissa == 0.5 && value.exponent == 2,
        "scaled-positive product returned the wrong value", 379);
}

static bool check_observer_contract;
typedef struct {
  int completed;
  ghl_error_codes_t result;
} completion_record;

static void check_completion_event(
      void *context, const ghl_m1_solver_stage_t stage,
      const ghl_error_codes_t status, const double U[4], const double residual[4],
      const ghl_m1_newton_diagnostics *diagnostics) {
  (void)residual;
  (void)diagnostics;
  completion_record *record = context;
  require_condition(U != NULL, "observer lost current iterate", 780);
  if(stage == ghl_m1_solver_stage_completed_solve) {
    record->completed++;
    record->result = status;
  }
}

static ghl_error_codes_t checked_newton_with_initial_guess(
      const ghl_m1_parameters *params, const ghl_metric_quantities *metric,
      const ghl_m1_newton_callbacks *callbacks, const void *context,
      const double base[4], const double initial[4], double output[4],
      ghl_m1_newton_diagnostics *diagnostics) {
  completion_record record = {0};
  ghl_m1_newton_callbacks observed;
  const bool instrument = check_observer_contract && callbacks != NULL
                          && callbacks->observer == NULL;
  if(instrument) {
    observed = *callbacks;
    observed.observer = check_completion_event;
    observed.observer_context = &record;
    callbacks = &observed;
  }
  const ghl_error_codes_t result = ghl_m1_newton_solve_4d_with_initial_guess(
      params, metric, callbacks, context, base, initial, output, diagnostics);
  if(instrument && params && metric && callbacks->residual && callbacks->jacobian
     && base && initial && output) {
    require_condition(record.completed == 1 && record.result == result,
                      "observer completion disagrees with solver result", 780);
  }
  return result;
}

static ghl_error_codes_t checked_newton(
      const ghl_m1_parameters *params, const ghl_metric_quantities *metric,
      const ghl_m1_newton_callbacks *callbacks, const void *context,
      const double base[4], double output[4], ghl_m1_newton_diagnostics *diagnostics) {
  if(!check_observer_contract) {
    return ghl_m1_newton_solve_4d(
        params, metric, callbacks, context, base, output, diagnostics);
  }
  return checked_newton_with_initial_guess(
      params, metric, callbacks, context, base, base, output, diagnostics);
}

static void make_primitives(ghl_primitive_quantities *restrict prims) {
  *prims = (ghl_primitive_quantities){ 0 };
  prims->rho = 1.0;
  prims->eps = 0.1;
  prims->press = 0.1;
  prims->Y_e = 0.5;
  prims->temperature = 1.0;
  prims->entropy = 1.0;
  prims->vU[1] = 0.2;
}

static void make_neutrino_parameters(ghl_m1_neutrino_parameters *restrict nu_params) {
  *nu_params = (ghl_m1_neutrino_parameters){
    .N_floor = 1.0e-12,
    .mean_energy_min = 0.0,
    .mean_energy_max = 0.0,
    .enforce_mean_energy_bounds = 0,
    .terminal_fallback_policy = ghl_m1_neutrino_terminal_fallback_no_update_all,
    .J_floor = 1.0e-14,
    .Gamma_N_floor = 1.0e-12
  };
}

static void
make_rates(const double target_scale, ghl_m1_neutrino_rates *restrict rates) {
  *rates = (ghl_m1_neutrino_rates){ 0 };
  rates->species = ghl_m1_neutrino_nux;
  rates->kappa_a_N = 1.0;
  rates->kappa_a_E = 1.0;
  rates->kappa_s = 1.0;
  rates->kappa_tr = 2.0;
  rates->n_eq = target_scale;
  rates->mean_energy = 1.0;
  rates->J_eq = target_scale;
  rates->eta_N = target_scale;
  rates->eta_E = target_scale;
  rates->lepton_weight = 0.0;
}

static void make_trigger_state(
      const double scale,
      const bool full_zero_flux,
      ghl_m1_neutrino_state *restrict state) {
  *state = (ghl_m1_neutrino_state){
    .N = scale, .E = scale, .F = { full_zero_flux ? 0.0 : 0.1 * scale, 0.0, 0.0 }
  };
}

static void make_U(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      double U[4]) {
  U[0] = state->E * metric->sqrt_detgamma;
  for(int direction = 0; direction < 3; ++direction) {
    U[direction + 1] = state->F[direction] * metric->sqrt_detgamma;
  }
}

static double coupled_fd_delta(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const double U_base[4],
      const double U[4]) {
  const double floor_scale = metric->sqrt_detgamma * m1_params->E_floor;
  const double energy_scale = fmax(fmax(fabs(U[0]), fabs(U_base[0])), floor_scale);
  return m1_params->fd_epsilon_rel * energy_scale
         + m1_params->fd_epsilon_abs * floor_scale;
}

static void compute_residual(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in,
      const double dt,
      const double U[4],
      double residual[4],
      const int case_index) {
  const ghl_error_codes_t error = ghl_m1_neutrino_compute_implicit_residual(
        m1_params, nu_params, metric, prims, rates, state_in, dt, U, residual);
  require_error(error, ghl_success, "implicit residual", case_index);
}

static void compute_jacobian(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in,
      const double dt,
      const double U[4],
      const double residual[4],
      double jacobian[4][4],
      const int case_index) {
  const ghl_error_codes_t error = ghl_m1_neutrino_compute_implicit_jacobian(
        m1_params, nu_params, metric, prims, rates, state_in, dt, U, residual, jacobian);
  require_error(error, ghl_success, "implicit Jacobian", case_index);
}

static void check_jacobian_against_independent_difference(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in,
      const double dt,
      const int column,
      const bool expect_one_sided,
      const int case_index) {
  double U_base[4];
  make_U(metric, state_in, U_base);
  const double U[4] = { U_base[0], U_base[1], U_base[2], U_base[3] };
  double residual_0[4];
  compute_residual(
        m1_params, nu_params, metric, prims, rates, state_in, dt, U, residual_0,
        case_index);

  double jacobian[4][4] = { { 0.0 } };
  compute_jacobian(
        m1_params, nu_params, metric, prims, rates, state_in, dt, U, residual_0,
        jacobian, case_index);

  const double delta = coupled_fd_delta(m1_params, metric, U_base, U);
  require_condition(
        isfinite(delta) && delta > 0.0, "coupled finite-difference step is invalid",
        case_index);

  double U_perturbed[4] = { U[0], U[1], U[2], U[3] };
  U_perturbed[column] += delta;
  double residual_perturbed[4] = { 0.0, 0.0, 0.0, 0.0 };
  ghl_error_codes_t error = ghl_m1_neutrino_compute_implicit_residual(
        m1_params, nu_params, metric, prims, rates, state_in, dt, U_perturbed,
        residual_perturbed);

  double used_delta = delta;
  if(expect_one_sided) {
    require_error(
          error, ghl_error_m1_implicit_admissibility, "forward boundary residual",
          case_index);
    U_perturbed[column] = U[column] - delta;
    compute_residual(
          m1_params, nu_params, metric, prims, rates, state_in, dt, U_perturbed,
          residual_perturbed, case_index);
    used_delta = -delta;

    const ghl_m1_neutrino_implicit_context context = {
      .m1_params = m1_params, .metric = metric, .prims_frozen = prims, .rates = rates
    };
    double validated_jacobian[4][4] = { { 0.0 } };
    require_error(
          ghl_m1_neutrino_compute_implicit_jacobian_validated(
                &context, dt, U_base, U, residual_0, validated_jacobian),
          ghl_success, "validated one-sided implicit Jacobian", case_index);
    for(int row = 0; row < 4; ++row) {
      for(int component = 0; component < 4; ++component) {
        require_close(
              validated_jacobian[row][component], jacobian[row][component], 0.0, 0.0,
              "checked/validated one-sided Jacobian", case_index);
      }
    }
  }
  else {
    require_error(error, ghl_success, "forward interior residual", case_index);
  }

  const double independent_derivative
        = (residual_perturbed[column] - residual_0[column]) / used_delta;
  require_condition(
        isfinite(independent_derivative) && fabs(independent_derivative) > 1.0e-6,
        "independent transverse derivative is not resolvable", case_index);
  require_condition(
        isfinite(jacobian[column][column]) && fabs(jacobian[column][column]) > 1.0e-6,
        "Jacobian transverse column is zero or nonfinite", case_index);
  require_close(
        jacobian[column][column], independent_derivative, 2.0e-8, 2.0e-10,
        "transverse Jacobian derivative", case_index);

  if(!expect_one_sided) {
    require_condition(
          fabs(jacobian[column][column] - 1.0) > 1.0e-3,
          "transverse source derivative was not retained", case_index);
  }
}

static void check_centered_full_fy_column(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in) {
  double U_base[4];
  make_U(metric, state_in, U_base);
  const double U[4] = { U_base[0], U_base[1], U_base[2], U_base[3] };
  double residual_0[4];
  compute_residual(
        m1_params, nu_params, metric, prims, rates, state_in, 1.0, U, residual_0, 120);
  double jacobian[4][4] = { { 0.0 } };
  compute_jacobian(
        m1_params, nu_params, metric, prims, rates, state_in, 1.0, U, residual_0,
        jacobian, 120);

  const double energy_scale = fmax(fabs(U[0]), fabs(U_base[0]));
  const double step_factors[2] = { 1.0e-6, 1.0e-7 };
  double centered_columns[2][4] = { { 0.0 } };
  for(int step_index = 0; step_index < 2; ++step_index) {
    const double half_step = step_factors[step_index] * energy_scale;
    double U_plus[4] = { U[0], U[1], U[2], U[3] };
    double U_minus[4] = { U[0], U[1], U[2], U[3] };
    U_plus[2] += half_step;
    U_minus[2] -= half_step;
    double residual_plus[4];
    double residual_minus[4];
    compute_residual(
          m1_params, nu_params, metric, prims, rates, state_in, 1.0, U_plus,
          residual_plus, 121 + step_index * 10);
    compute_residual(
          m1_params, nu_params, metric, prims, rates, state_in, 1.0, U_minus,
          residual_minus, 122 + step_index * 10);
    for(int row = 0; row < 4; ++row) {
      centered_columns[step_index][row]
            = (residual_plus[row] - residual_minus[row]) / (2.0 * half_step);
      require_condition(
            isfinite(centered_columns[step_index][row]),
            "centered transverse derivative is nonfinite", 123 + step_index * 10);
      require_close(
            jacobian[row][2], centered_columns[step_index][row], 5.0e-5,
            step_index == 0 ? 2.0e-5 : 3.0e-6, "full transverse Jacobian column",
            124 + step_index * 10 + row);
    }
  }

  for(int row = 0; row < 4; ++row) {
    require_close(
          centered_columns[0][row], centered_columns[1][row], 5.0e-5, 2.0e-5,
          "independent centered transverse columns", 145 + row);
  }
  require_condition(
        fabs(centered_columns[0][2]) > 1.0e-6,
        "centered transverse derivative is not resolvable", 150);
}

static void check_full_zero_flux_columns(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in,
      const double dt) {
  double U_base[4];
  make_U(metric, state_in, U_base);
  const double U[4] = { U_base[0], U_base[1], U_base[2], U_base[3] };
  double residual[4];
  compute_residual(
        m1_params, nu_params, metric, prims, rates, state_in, dt, U, residual, 200);
  double jacobian[4][4] = { { 0.0 } };
  compute_jacobian(
        m1_params, nu_params, metric, prims, rates, state_in, dt, U, residual, jacobian,
        200);

  for(int column = 1; column <= 3; ++column) {
    require_condition(
          isfinite(jacobian[column][column]) && fabs(jacobian[column][column]) > 1.0e-6,
          "full-zero flux diagonal Jacobian entry is not resolved", 200 + column);
  }
}

static void check_scaled_case(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in,
      const double dt,
      const double reference_derivative) {
  double U_base[4];
  make_U(metric, state_in, U_base);
  const double U[4] = { U_base[0], U_base[1], U_base[2], U_base[3] };
  double residual[4];
  compute_residual(
        m1_params, nu_params, metric, prims, rates, state_in, dt, U, residual, 300);
  double jacobian[4][4] = { { 0.0 } };
  compute_jacobian(
        m1_params, nu_params, metric, prims, rates, state_in, dt, U, residual, jacobian,
        300);
  require_close(
        jacobian[2][2], reference_derivative, 3.0e-8, 3.0e-10,
        "common-rescaling transverse derivative", 300);
}

static void check_validated_context_equivalence(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in) {
  double U_base[4];
  make_U(metric, state_in, U_base);
  const double U[4] = { U_base[0], U_base[1], U_base[2], U_base[3] };
  double checked_residual[4];
  compute_residual(
        m1_params, nu_params, metric, prims, rates, state_in, 0.75, U, checked_residual,
        350);

  const ghl_m1_neutrino_implicit_context context = {
    .m1_params = m1_params, .metric = metric, .prims_frozen = prims, .rates = rates
  };
  double validated_residual[4] = { 0.0, 0.0, 0.0, 0.0 };
  bool closure_fallback_observed = false;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_validated(
              &context, 0.75, U_base, U, &closure_fallback_observed, validated_residual),
        ghl_success, "validated implicit residual", 350);
  for(int component = 0; component < 4; ++component) {
    require_close(
          validated_residual[component], checked_residual[component], 0.0, 0.0,
          "checked/validated residual", 351 + component);
  }

  double checked_jacobian[4][4] = { { 0.0 } };
  compute_jacobian(
        m1_params, nu_params, metric, prims, rates, state_in, 0.75, U, checked_residual,
        checked_jacobian, 356);
  double validated_jacobian[4][4] = { { 0.0 } };
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_validated(
              &context, 0.75, U_base, U, checked_residual, validated_jacobian),
        ghl_success, "validated implicit Jacobian", 357);
  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      require_close(
            validated_jacobian[row][column], checked_jacobian[row][column], 0.0, 0.0,
            "checked/validated Jacobian", 358 + row * 4 + column);
    }
  }

  require_error(
        ghl_m1_neutrino_compute_implicit_residual_validated(
              NULL, 0.75, U_base, U, NULL, validated_residual),
        ghl_error_m1_null_pointer, "null validated residual context", 375);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_validated(
              NULL, 0.75, U_base, U, checked_residual, validated_jacobian),
        ghl_error_m1_null_pointer, "null validated Jacobian context", 376);
}

static void check_invalid_fd_step_rejection(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in) {
  double U[4];
  make_U(metric, state_in, U);
  double residual[4];
  compute_residual(
        m1_params, nu_params, metric, prims, rates, state_in, 0.5, U, residual, 380);
  double nonfinite_residual[4] = { residual[0], residual[1], residual[2], residual[3] };
  nonfinite_residual[0] = NAN;
  double jacobian[4][4] = { { 0.0 } };
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, metric, prims, rates, 0.5, U, U, nonfinite_residual, jacobian),
        ghl_error_m1_invalid_state, "nonfinite reference residual", 381);
  ghl_m1_parameters zero_step_params = *m1_params;
  zero_step_params.fd_epsilon_rel = 0.0;
  zero_step_params.fd_epsilon_abs = 0.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              &zero_step_params, metric, prims, rates, 0.5, U, U, residual, jacobian),
        ghl_error_m1_invalid_state, "zero finite-difference step", 380);

  ghl_m1_parameters nonfinite_step_params = *m1_params;
  nonfinite_step_params.fd_epsilon_rel = INFINITY;
  nonfinite_step_params.fd_epsilon_abs = 0.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              &nonfinite_step_params, metric, prims, rates, 0.5, U, U, residual, jacobian),
        ghl_error_m1_invalid_state, "nonfinite finite-difference step", 380);

  ghl_metric_quantities overflow_metric = *metric;
  overflow_metric.lapse = 2.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &overflow_metric, prims, rates, DBL_MAX, U, U, residual),
        ghl_error_m1_invalid_state, "overflowing dt-lapse product", 380);

  ghl_metric_quantities invalid_volume_metric = *metric;
  invalid_volume_metric.sqrt_detgamma = 0.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual(
              m1_params, nu_params, &invalid_volume_metric, prims, rates, state_in, 0.5, U,
              residual),
        ghl_error_m1_invalid_metric, "zero public residual volume", 380);
  invalid_volume_metric.sqrt_detgamma = NAN;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual(
              m1_params, nu_params, &invalid_volume_metric, prims, rates, state_in, 0.5, U,
              residual),
        ghl_error_m1_invalid_metric, "nonfinite public residual volume", 380);

  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, metric, prims, rates, -1.0, U, U, residual, jacobian),
        ghl_error_m1_invalid_state, "forward residual hard-error propagation", 382);

  /* An overflowing forward energy perturbation permits a backward retry.
   * The backward trial is finite, but its invalid frozen rates must retain
   * their microphysics error instead of becoming a Jacobian-domain error. */
  ghl_m1_parameters large_step_params = *m1_params;
  large_step_params.fd_epsilon_rel = 0.5;
  ghl_m1_neutrino_rates invalid_rates = *rates;
  invalid_rates.eta_E = NAN;
  const double large_U[4] = { DBL_MAX, 0.0, 0.0, 0.0 };
  const double zero_residual[4] = { 0.0, 0.0, 0.0, 0.0 };
  double forward_U[4] = { INFINITY, 0.0, 0.0, 0.0 };
  double backward_U[4] = { DBL_MAX / 2.0, 0.0, 0.0, 0.0 };
  double trial_residual[4];
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              &large_step_params, metric, prims, &invalid_rates, 0.0, large_U, forward_U,
              trial_residual),
        ghl_error_m1_implicit_admissibility, "overflowing forward energy trial", 383);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              &large_step_params, metric, prims, &invalid_rates, 0.0, large_U,
              backward_U, trial_residual),
        ghl_error_m1_microphysics_failure, "backward trial invalid frozen rates", 384);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              &large_step_params, metric, prims, &invalid_rates, 0.0, large_U, large_U,
              zero_residual, jacobian),
        ghl_error_m1_microphysics_failure, "backward residual hard-error propagation",
        385);

  /* A positive step below the spacing at DBL_MAX must advance by one ULP.
   * The forward ULP overflows; the backward ULP remains finite and preserves
   * the hard error from the frozen rates. */
  ghl_m1_parameters sub_ulp_params = *m1_params;
  sub_ulp_params.fd_epsilon_rel = nextafter(0.0, 1.0);
  sub_ulp_params.fd_epsilon_abs = 0.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              &sub_ulp_params, metric, prims, &invalid_rates, 0.0, large_U, large_U,
              zero_residual, jacobian),
        ghl_error_m1_microphysics_failure, "sub-ULP one-sided hard-error propagation",
        386);
}

static void check_public_convergence(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in) {
  ghl_m1_neutrino_state state_out = { .N = -1.0, .E = -2.0, .F = { -3.0, -4.0, -5.0 } };
  ghl_m1_neutrino_exchange exchange = { .dE_rad = -1.0 };
  ghl_m1_implicit_solve_diagnostics solve_diagnostics;
  ghl_m1_neutrino_diagnostics neutrino_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);

  const ghl_error_codes_t error = ghl_m1_solve_neutrino_implicit_homogeneous_update(
        m1_params, nu_params, metric, prims, rates, 1.0, 1.0, state_in, &state_out,
        &exchange, &solve_diagnostics, &neutrino_diagnostics);
  require_error(error, ghl_success, "public implicit source solve", 400);
  require_condition(
        solve_diagnostics.residual_scaled_norm <= 1.0 + 1.0e-10,
        "solver reported a residual above its requested tolerance", 400);

  double U_base[4];
  make_U(metric, state_in, U_base);
  double U_final[4];
  make_U(metric, &state_out, U_final);
  double residual[4];
  compute_residual(
        m1_params, nu_params, metric, prims, rates, state_in, 1.0, U_final, residual,
        401);

  const double floor_scale = metric->sqrt_detgamma * m1_params->E_floor;
  const double energy_scale = fmax(fmax(fabs(U_final[0]), fabs(U_base[0])), floor_scale);
  for(int i = 0; i < 4; ++i) {
    const double component_scale
          = i == 0 ? fmax(fabs(U_final[i]), fmax(fabs(U_base[i]), floor_scale))
                   : fmax(fabs(U_final[i]), fmax(fabs(U_base[i]), energy_scale));
    const double requested_tolerance
          = metric->sqrt_detgamma * m1_params->newton_absolute_tolerance
            + m1_params->newton_tolerance * component_scale;
    require_condition(
          isfinite(residual[i])
                && fabs(residual[i]) <= requested_tolerance * (1.0 + 1.0e-8),
          "final residual exceeds the requested component tolerance", 401 + i);
  }
  require_condition(
        state_out.F[1] > 1.0e-3 && fabs(state_out.E - state_in->E) > 1.0e-6,
        "successful solve did not publish the interaction update", 405);
}

typedef struct {
  double target[4];
} affine_newton_context;

static ghl_error_codes_t affine_newton_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  const affine_newton_context *const affine = (const affine_newton_context *)context;
  for(int component = 0; component < 4; ++component) {
    residual[component] = U[component] - affine->target[component];
  }
  return ghl_success;
}

static ghl_error_codes_t affine_newton_jacobian(
      const void *restrict context,
      const double U[4],
      const double residual[4],
      double jacobian[4][4]) {
  (void)context;
  (void)U;
  (void)residual;
  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      jacobian[row][column] = row == column ? 1.0 : 0.0;
    }
  }
  return ghl_success;
}

static void set_identity_matrix(double matrix[4][4]);

static ghl_error_codes_t permuted_affine_newton_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  const affine_newton_context *const affine = (const affine_newton_context *)context;
  residual[0] = U[1] - affine->target[1];
  residual[1] = U[0] - affine->target[0];
  residual[2] = U[2] - affine->target[2];
  residual[3] = U[3] - affine->target[3];
  return ghl_success;
}

static ghl_error_codes_t permuted_affine_newton_jacobian(
      const void *restrict context,
      const double U[4],
      const double residual[4],
      double jacobian[4][4]) {
  (void)context;
  (void)U;
  (void)residual;
  set_identity_matrix(jacobian);
  jacobian[0][0] = 0.0;
  jacobian[0][1] = 1.0;
  jacobian[1][0] = 1.0;
  jacobian[1][1] = 0.0;
  return ghl_success;
}

typedef struct {
  double residual[4];
  double jacobian[4][4];
  ghl_error_codes_t residual_status;
  ghl_error_codes_t jacobian_status;
} fixed_newton_context;

static void set_identity_matrix(double matrix[4][4]) {
  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      matrix[row][column] = row == column ? 1.0 : 0.0;
    }
  }
}

static ghl_error_codes_t fixed_newton_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  (void)U;
  const fixed_newton_context *const fixed = (const fixed_newton_context *)context;
  for(int component = 0; component < 4; ++component) {
    residual[component] = fixed->residual[component];
  }
  return fixed->residual_status;
}

static ghl_error_codes_t fixed_newton_jacobian(
      const void *restrict context,
      const double U[4],
      const double residual[4],
      double jacobian[4][4]) {
  (void)U;
  (void)residual;
  const fixed_newton_context *const fixed = (const fixed_newton_context *)context;
  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      jacobian[row][column] = fixed->jacobian[row][column];
    }
  }
  return fixed->jacobian_status;
}

static void expect_solver_failure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const fixed_newton_context *restrict context,
      const char *restrict operation,
      const int case_index) {
  const ghl_m1_newton_callbacks callbacks = { .residual = fixed_newton_residual,
                                              .jacobian = fixed_newton_jacobian,
                                              .observer = NULL,
                                              .observer_context = NULL };
  const double U_base[4] = { 0.0, 0.0, 0.0, 0.0 };
  double U_out[4] = { -1.0, -2.0, -3.0, -4.0 };
  ghl_m1_newton_diagnostics diagnostics;
  require_error(
        checked_newton(
              m1_params, metric, &callbacks, context, U_base, U_out, &diagnostics),
        ghl_error_m1_implicit_solve_failure, operation, case_index);
}

static void expect_solver_error(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const fixed_newton_context *restrict context,
      const ghl_error_codes_t expected,
      const char *restrict operation,
      const int case_index) {
  const ghl_m1_newton_callbacks callbacks = { .residual = fixed_newton_residual,
                                              .jacobian = fixed_newton_jacobian,
                                              .observer = NULL,
                                              .observer_context = NULL };
  const double U_base[4] = { 0.0, 0.0, 0.0, 0.0 };
  double U_out[4] = { -1.0, -2.0, -3.0, -4.0 };
  ghl_m1_newton_diagnostics diagnostics;
  require_error(
        checked_newton(
              m1_params, metric, &callbacks, context, U_base, U_out, &diagnostics),
        expected, operation, case_index);
}

static void test_weighted_merit_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {
  const double U[4] = { 2.0, 3.0, 0.0, 0.0 };
  const double U_base[4] = { 1.0, 0.0, 0.0, 0.0 };
  const double residual[4] = { 0.0, 0.0, 0.0, 0.0 };
  const double merit
        = ghl_m1_newton_weighted_merit(m1_params, metric, U, U_base, residual);
  require_condition(
        isfinite(merit) && merit == 0.0, "finite weighted merit boundary failed", 720);

  ghl_m1_parameters floor_params = *m1_params;
  floor_params.E_floor = 1.0;
  const double zero[4] = { 0.0, 0.0, 0.0, 0.0 };
  require_condition(
        ghl_m1_newton_weighted_merit(&floor_params, metric, zero, zero, zero) == 0.0,
        "energy-floor weighted merit boundary failed", 721);

  ghl_m1_parameters nonfinite_tolerance_params = *m1_params;
  nonfinite_tolerance_params.newton_absolute_tolerance = INFINITY;
  require_condition(
        isinf(ghl_m1_newton_weighted_merit(
              &nonfinite_tolerance_params, metric, U, U_base, residual)),
        "nonfinite weighted-merit denominator was not rejected", 722);

  ghl_m1_parameters negative_tolerance_params = *m1_params;
  negative_tolerance_params.newton_absolute_tolerance = -1.0;
  negative_tolerance_params.newton_tolerance = 0.0;
  require_condition(
        isinf(ghl_m1_newton_weighted_merit(
              &negative_tolerance_params, metric, U, U_base, residual)),
        "nonpositive weighted-merit denominator was not rejected", 723);

  double nonfinite_residual[4] = { 0.0, 0.0, NAN, 0.0 };
  require_condition(
        isinf(ghl_m1_newton_weighted_merit(
              m1_params, metric, U, U_base, nonfinite_residual)),
        "nonfinite weighted-merit residual was not rejected", 724);

  /* This is an exported entry point and must report "no admissible merit"
   * rather than dereferencing a missing operand. */
  require_condition(
        isinf(ghl_m1_newton_weighted_merit(NULL, metric, U, U_base, residual))
              && isinf(
                    ghl_m1_newton_weighted_merit(m1_params, NULL, U, U_base, residual))
              && isinf(ghl_m1_newton_weighted_merit(
                    m1_params, metric, NULL, U_base, residual))
              && isinf(
                    ghl_m1_newton_weighted_merit(m1_params, metric, U, NULL, residual))
              && isinf(ghl_m1_newton_weighted_merit(m1_params, metric, U, U_base, NULL)),
        "NULL weighted-merit operands were not rejected", 725);
}

static void test_projection_failure_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {
  const double original[4] = { 1.0, 2.0, 3.0, 4.0 };
  for(int component = 0; component < 4; ++component) {
    double projected[4] = { original[0], original[1], original[2], original[3] };
    projected[component] = NAN;
    const double before[4] = { projected[0], projected[1], projected[2], projected[3] };
    require_error(
          ghl_m1_newton_project_admissible(m1_params, metric, projected),
          ghl_error_m1_implicit_admissibility,
          "nonfinite admissibility projection input", 730 + component);
    require_condition(
          memcmp(projected, before, sizeof(projected)) == 0,
          "nonfinite admissibility projection changed its input", 730 + component);
  }

  ghl_metric_quantities invalid_tensor_metric = *metric;
  invalid_tensor_metric.gammaUU[0][0] = NAN;
  double invalid_tensor_projection[4]
        = { original[0], original[1], original[2], original[3] };
  require_error(
        ghl_m1_newton_project_admissible(
              m1_params, &invalid_tensor_metric, invalid_tensor_projection),
        ghl_error_m1_invalid_metric, "invalid repair-metric projection", 734);
  require_condition(
        memcmp(invalid_tensor_projection, original, sizeof(invalid_tensor_projection))
              == 0,
        "invalid repair-metric projection changed its input", 734);

  ghl_metric_quantities nonfinite_metric = *metric;
  nonfinite_metric.sqrt_detgamma = NAN;
  double nonfinite_metric_projection[4]
        = { original[0], original[1], original[2], original[3] };
  require_error(
        ghl_m1_newton_project_admissible(
              m1_params, &nonfinite_metric, nonfinite_metric_projection),
        ghl_error_m1_invalid_metric, "nonfinite projection metric", 736);
  require_condition(
        memcmp(
              nonfinite_metric_projection, original, sizeof(nonfinite_metric_projection))
              == 0,
        "nonfinite projection metric changed its input", 736);

  ghl_m1_parameters overflow_params = *m1_params;
  overflow_params.E_floor = DBL_MAX;
  ghl_metric_quantities overflow_metric = *metric;
  overflow_metric.detgamma = 2.0;
  overflow_metric.sqrt_detgamma = sqrt(2.0);
  overflow_metric.gammaDD[0][0] = 2.0;
  overflow_metric.gammaUU[0][0] = 0.5;
  double overflow_projection[4] = { 1.0, 0.0, 0.0, 0.0 };
  const double overflow_before[4] = { overflow_projection[0], overflow_projection[1],
                                      overflow_projection[2], overflow_projection[3] };
  require_error(
        ghl_m1_newton_project_admissible(
              &overflow_params, &overflow_metric, overflow_projection),
        ghl_error_m1_implicit_admissibility, "overflowing admissibility projection",
        735);
  require_condition(
        memcmp(overflow_projection, overflow_before, sizeof(overflow_projection)) == 0,
        "overflowing admissibility projection published partial output", 735);
}

static void test_linear_rejection_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {
  fixed_newton_context context = { .residual = { 1.0, 0.0, 0.0, 0.0 },
                                   .residual_status = ghl_success,
                                   .jacobian_status = ghl_success };

  set_identity_matrix(context.jacobian);
  context.residual[0] = NAN;
  expect_solver_failure(m1_params, metric, &context, "nonfinite Newton RHS", 740);

  context.residual[0] = 1.0;
  context.jacobian[0][0] = NAN;
  expect_solver_failure(m1_params, metric, &context, "nonfinite Newton matrix", 741);

  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      context.jacobian[row][column] = 0.0;
    }
  }
  expect_solver_failure(m1_params, metric, &context, "zero Newton matrix row", 742);

  for(int row = 0; row < 4; ++row) {
    context.jacobian[row][0] = 1.0;
  }
  expect_solver_failure(m1_params, metric, &context, "rank-deficient Newton pivot", 743);

  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      context.jacobian[row][column] = row == column ? 1.0 : 0.0;
    }
  }
  context.jacobian[3][3] = nextafter(DBL_MIN, 0.0);
  expect_solver_failure(m1_params, metric, &context, "subnormal Newton pivot", 744);

  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      context.jacobian[row][column] = 0.0;
    }
  }
  /* Finite input entries can overflow during elimination. The scaled-pivot
   * guard must reject the resulting infinite pivot before back substitution. */
  context.jacobian[0][0] = 1.0;
  context.jacobian[0][1] = 1.0;
  context.jacobian[1][0] = -0x1p971;
  context.jacobian[1][1] = DBL_MAX;
  context.jacobian[2][2] = 1.0;
  context.jacobian[3][3] = 1.0;
  expect_solver_failure(
        m1_params, metric, &context, "finite matrix with overflowing Newton pivot", 746);

  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      context.jacobian[row][column] = 0.0;
    }
  }
  affine_newton_context permuted_affine = { .target = { 1.25, 0.10, -0.20, 0.05 } };
  const ghl_m1_newton_callbacks permuted_callbacks
        = { .residual = permuted_affine_newton_residual,
            .jacobian = permuted_affine_newton_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  const double permuted_base[4] = { 1.0, 0.0, 0.0, 0.0 };
  double permuted_output[4] = { -1.0, -2.0, -3.0, -4.0 };
  ghl_m1_newton_diagnostics permuted_diagnostics;
  require_error(
        checked_newton(
              m1_params, metric, &permuted_callbacks, &permuted_affine, permuted_base,
              permuted_output, &permuted_diagnostics),
        ghl_success, "row-swapped Newton matrix solve", 745);
  for(int component = 0; component < 4; ++component) {
    require_close(
          permuted_output[component], permuted_affine.target[component], 0.0, 0.0,
          "row-swapped Newton solution", 745 + component);
  }

  context.residual_status = ghl_error_m1_invalid_state;
  expect_solver_error(
        m1_params, metric, &context, ghl_error_m1_invalid_state,
        "initial Newton residual callback failure", 749);
  context.residual_status = ghl_success;

  context.jacobian_status = ghl_error_m1_invalid_state;
  expect_solver_error(
        m1_params, metric, &context, ghl_error_m1_invalid_state,
        "Newton Jacobian callback failure", 750);
  context.jacobian_status = ghl_success;

  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      context.jacobian[row][column] = 0.0;
    }
  }
  context.jacobian[0][0] = DBL_MIN;
  context.jacobian[1][0] = DBL_MAX;
  context.jacobian[1][1] = 1.0;
  context.jacobian[2][2] = 1.0;
  context.jacobian[3][3] = 1.0;
  expect_solver_failure(
        m1_params, metric, &context, "overflowing Newton elimination factor", 746);

  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      context.jacobian[row][column] = 0.0;
    }
  }
  context.jacobian[0][0] = 1.0e-200;
  context.jacobian[0][1] = 1.0e300;
  context.jacobian[1][0] = 1.0e-100;
  context.jacobian[1][1] = DBL_MAX;
  context.jacobian[2][2] = 1.0;
  context.jacobian[3][3] = 1.0;
  expect_solver_failure(
        m1_params, metric, &context, "nonfinite scaled Newton pivot after elimination",
        749);

  for(int row = 0; row < 4; ++row) {
    for(int column = 0; column < 4; ++column) {
      context.jacobian[row][column] = row == column ? 1.0 : 0.0;
    }
  }
  context.jacobian[0][0] = DBL_MAX;
  context.jacobian[0][1] = 2.0;
  context.residual[1] = -DBL_MAX;
  expect_solver_failure(
        m1_params, metric, &context, "overflowing Newton back-substitution", 747);

  context.residual[1] = 0.0;
  context.jacobian[3][3] = DBL_MIN;
  context.residual[3] = -DBL_MAX;
  expect_solver_failure(m1_params, metric, &context, "overflowing Newton solution", 748);
}

static void test_newton_input_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {
  const affine_newton_context affine = { .target = { 1.25, 0.10, -0.20, 0.05 } };
  const ghl_m1_newton_callbacks callbacks = { .residual = affine_newton_residual,
                                              .jacobian = affine_newton_jacobian,
                                              .observer = NULL,
                                              .observer_context = NULL };
  const double U_base[4] = { 1.0, 0.0, 0.0, 0.0 };
  const double U_initial[4] = { 1.0, 0.0, 0.0, 0.0 };
  double U_out[4] = { -1.0, -2.0, -3.0, -4.0 };
  ghl_m1_newton_diagnostics diagnostics;

  require_error(
        checked_newton_with_initial_guess(
              NULL, metric, &callbacks, &affine, U_base, U_initial, U_out, &diagnostics),
        ghl_error_m1_null_pointer, "NULL Newton parameters", 750);
  require_error(
        checked_newton_with_initial_guess(
              m1_params, NULL, &callbacks, &affine, U_base, U_initial, U_out,
              &diagnostics),
        ghl_error_m1_null_pointer, "NULL Newton metric", 751);
  require_error(
        checked_newton_with_initial_guess(
              m1_params, metric, NULL, &affine, U_base, U_initial, U_out, &diagnostics),
        ghl_error_m1_null_pointer, "NULL Newton callbacks", 752);

  ghl_m1_newton_callbacks null_residual_callbacks = callbacks;
  null_residual_callbacks.residual = NULL;
  require_error(
        checked_newton_with_initial_guess(
              m1_params, metric, &null_residual_callbacks, &affine, U_base, U_initial,
              U_out, &diagnostics),
        ghl_error_m1_null_pointer, "NULL Newton residual callback", 753);

  ghl_m1_newton_callbacks null_jacobian_callbacks = callbacks;
  null_jacobian_callbacks.jacobian = NULL;
  require_error(
        checked_newton_with_initial_guess(
              m1_params, metric, &null_jacobian_callbacks, &affine, U_base, U_initial,
              U_out, &diagnostics),
        ghl_error_m1_null_pointer, "NULL Newton Jacobian callback", 754);
  require_error(
        checked_newton_with_initial_guess(
              m1_params, metric, &callbacks, &affine, NULL, U_initial, U_out,
              &diagnostics),
        ghl_error_m1_null_pointer, "NULL Newton base state", 755);
  require_error(
        checked_newton_with_initial_guess(
              m1_params, metric, &callbacks, &affine, U_base, NULL, U_out, &diagnostics),
        ghl_error_m1_null_pointer, "NULL Newton initial state", 756);
  require_error(
        checked_newton_with_initial_guess(
              m1_params, metric, &callbacks, &affine, U_base, U_initial, NULL,
              &diagnostics),
        ghl_error_m1_null_pointer, "NULL Newton output state", 757);
}

typedef struct {
  int calls;
  int completed_calls;
  int null_residual_calls;
} newton_observer_context;

static void record_newton_event(
      void *restrict context,
      const ghl_m1_solver_stage_t stage,
      const ghl_error_codes_t status,
      const double U[4],
      const double residual[4],
      const ghl_m1_newton_diagnostics *restrict diagnostics) {
  (void)status;
  (void)U;
  (void)diagnostics;
  newton_observer_context *const observer = (newton_observer_context *)context;
  observer->calls++;
  if(stage == ghl_m1_solver_stage_completed_solve) {
    observer->completed_calls++;
  }
  if(residual == NULL) {
    observer->null_residual_calls++;
  }
}

static void test_newton_optional_paths(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {
  const affine_newton_context affine = { .target = { 1.25, 0.10, -0.20, 0.05 } };
  const ghl_m1_newton_callbacks callbacks_without_observer
        = { .residual = affine_newton_residual,
            .jacobian = affine_newton_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  const double U_base[4] = { 1.0, 0.0, 0.0, 0.0 };
  double U_out[4] = { -1.0, -2.0, -3.0, -4.0 };
  require_error(
        checked_newton(
              m1_params, metric, &callbacks_without_observer, &affine, U_base, U_out,
              NULL),
        ghl_success, "Newton solve without diagnostics", 760);
  for(int component = 0; component < 4; ++component) {
    require_close(
          U_out[component], affine.target[component], 0.0, 0.0,
          "Newton solution without diagnostics", 760);
  }

  newton_observer_context observer_context = { 0 };
  const ghl_m1_newton_callbacks callbacks_with_observer
        = { .residual = affine_newton_residual,
            .jacobian = affine_newton_jacobian,
            .observer = record_newton_event,
            .observer_context = &observer_context };
  ghl_m1_newton_diagnostics diagnostics;
  require_error(
        checked_newton(
              m1_params, metric, &callbacks_with_observer, &affine, U_base, U_out,
              &diagnostics),
        ghl_success, "Newton solve with observer", 761);
  for(int component = 0; component < 4; ++component) {
    require_close(
          U_out[component], affine.target[component], 0.0, 0.0,
          "Newton solution with observer", 761);
  }
  require_condition(
        observer_context.calls > 0 && observer_context.completed_calls > 0,
        "Newton observer did not receive solve events", 761);
}

typedef struct {
  int calls;
} scripted_newton_context;

static ghl_error_codes_t zero_residual_for_any_state(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  (void)U;
  scripted_newton_context *const scripted = (scripted_newton_context *)context;
  scripted->calls++;
  for(int component = 0; component < 4; ++component) {
    residual[component] = 0.0;
  }
  return ghl_success;
}

static ghl_error_codes_t scripted_newton_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  scripted_newton_context *const scripted = (scripted_newton_context *)context;
  scripted->calls++;
  residual[0] = 0.0;
  residual[1] = residual[2] = residual[3] = 0.0;
  if(scripted->calls == 1) {
    residual[0] = 1.0;
    return ghl_success;
  }
  if(U[0] >= 1.75) {
    residual[0] = 2.0;
  }
  else {
    residual[0] = 0.1;
  }
  return ghl_success;
}

static ghl_error_codes_t always_admissible_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  (void)U;
  scripted_newton_context *const scripted = (scripted_newton_context *)context;
  scripted->calls++;
  for(int component = 0; component < 4; ++component) {
    residual[component] = 0.0;
  }
  if(scripted->calls == 1) {
    residual[0] = 1.0;
    return ghl_success;
  }
  return ghl_error_m1_implicit_admissibility;
}

static ghl_error_codes_t post_projection_error_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  (void)U;
  scripted_newton_context *const scripted = (scripted_newton_context *)context;
  scripted->calls++;
  for(int component = 0; component < 4; ++component) {
    residual[component] = 0.0;
  }
  if(scripted->calls == 1) {
    residual[0] = 1.0;
    return ghl_success;
  }
  /* Calls 2-12 are the backtracking ladder (step scales 1 .. 1/1024); only
   * after every shortened step fails does the driver project the full step. */
  if(scripted->calls <= 12) {
    return ghl_error_m1_implicit_admissibility;
  }
  return ghl_error_m1_invalid_state;
}

static ghl_error_codes_t post_projection_success_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  (void)U;
  scripted_newton_context *const scripted = (scripted_newton_context *)context;
  scripted->calls++;
  for(int component = 0; component < 4; ++component) {
    residual[component] = 0.0;
  }
  if(scripted->calls == 1) {
    residual[0] = 1.0;
    return ghl_success;
  }
  /* Calls 2-12 are the backtracking ladder (step scales 1 .. 1/1024); only
   * after every shortened step fails does the driver project the full step. */
  if(scripted->calls <= 12) {
    return ghl_error_m1_implicit_admissibility;
  }
  return ghl_success;
}

static ghl_error_codes_t overflow_trial_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  (void)U;
  scripted_newton_context *const scripted = (scripted_newton_context *)context;
  scripted->calls++;
  for(int component = 0; component < 4; ++component) {
    residual[component] = 0.0;
  }
  if(scripted->calls == 1) {
    residual[0] = -DBL_MAX;
    return ghl_success;
  }
  return ghl_error_m1_implicit_admissibility;
}

static ghl_error_codes_t trial_callback_error_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  (void)U;
  scripted_newton_context *const scripted = (scripted_newton_context *)context;
  scripted->calls++;
  for(int component = 0; component < 4; ++component) {
    residual[component] = 0.0;
  }
  if(scripted->calls == 1) {
    residual[0] = 1.0;
    return ghl_success;
  }
  return ghl_error_m1_invalid_state;
}

static ghl_error_codes_t scripted_identity_jacobian(
      const void *restrict context,
      const double U[4],
      const double residual[4],
      double jacobian[4][4]) {
  (void)context;
  (void)U;
  (void)residual;
  set_identity_matrix(jacobian);
  return ghl_success;
}

static ghl_error_codes_t backtracking_jacobian(
      const void *restrict context,
      const double U[4],
      const double residual[4],
      double jacobian[4][4]) {
  (void)context;
  (void)U;
  (void)residual;
  set_identity_matrix(jacobian);
  jacobian[0][0] = -1.0;
  return ghl_success;
}

static ghl_error_codes_t small_step_jacobian(
      const void *restrict context,
      const double U[4],
      const double residual[4],
      double jacobian[4][4]) {
  (void)context;
  (void)U;
  (void)residual;
  set_identity_matrix(jacobian);
  /* Residual 1 needs a 1e-12 correction, inside the configured 1e-10 tolerance. */
  jacobian[0][0] = 1.0e12;
  return ghl_success;
}

static ghl_error_codes_t no_root_small_step_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  (void)context;
  const double t = 1.0e12 * (U[0] - 1.0);
  residual[0] = 1.0 + t + t * t;
  for(int i = 1; i < 4; ++i) {
    residual[i] = U[i];
  }
  return ghl_success;
}

static ghl_error_codes_t no_root_small_step_jacobian(
      const void *restrict context,
      const double U[4],
      const double residual[4],
      double jacobian[4][4]) {
  (void)context;
  (void)residual;
  set_identity_matrix(jacobian);
  const double t = 1.0e12 * (U[0] - 1.0);
  jacobian[0][0] = 1.0e12 * (1.0 + 2.0 * t);
  return ghl_success;
}

static ghl_error_codes_t floor_limited_small_step_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {
  (void)context;
  for(int i = 0; i < 4; ++i) {
    residual[i] = 0.0;
  }
  if(U[0] < 1.0) {
    return ghl_error_m1_implicit_admissibility;
  }
  residual[0] = 1.0 + 1.0e12 * (U[0] - 1.0);
  for(int i = 1; i < 4; ++i) {
    residual[i] = U[i];
  }
  return ghl_success;
}

static void test_newton_retry_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {
  const double base_state[4] = { 1.0, 0.0, 0.0, 0.0 };
  double output[4] = { -1.0, -2.0, -3.0, -4.0 };
  ghl_m1_newton_diagnostics diagnostics;

  scripted_newton_context backtracking_context = { 0 };
  ghl_m1_parameters one_iteration_params = *m1_params;
  one_iteration_params.newton_max_iterations = 1;
  const ghl_m1_newton_callbacks backtracking_callbacks
        = { .residual = scripted_newton_residual,
            .jacobian = backtracking_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  require_error(
        checked_newton(
              &one_iteration_params, metric, &backtracking_callbacks,
              &backtracking_context, base_state, output, &diagnostics),
        ghl_error_m1_implicit_solve_failure, "Newton backtracking and exhaustion", 770);
  require_condition(
        diagnostics.backtracks > 0, "Newton backtracking case did not backtrack", 770);

  scripted_newton_context admissible_context = { 0 };
  const ghl_m1_newton_callbacks admissible_callbacks
        = { .residual = always_admissible_residual,
            .jacobian = scripted_identity_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  require_error(
        checked_newton(
              m1_params, metric, &admissible_callbacks, &admissible_context, base_state,
              output, &diagnostics),
        ghl_error_m1_implicit_solve_failure, "Newton admissibility retry exhaustion",
        771);

  ghl_metric_quantities invalid_projection_metric = *metric;
  invalid_projection_metric.sqrt_detgamma = 0.0;
  scripted_newton_context invalid_projection_context = { 0 };
  require_error(
        checked_newton(
              m1_params, &invalid_projection_metric, &admissible_callbacks,
              &invalid_projection_context, base_state, output, &diagnostics),
        ghl_error_m1_invalid_metric, "Newton projection metric failure", 772);

  scripted_newton_context post_projection_context = { 0 };
  const ghl_m1_newton_callbacks post_projection_callbacks
        = { .residual = post_projection_error_residual,
            .jacobian = scripted_identity_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  require_error(
        checked_newton(
              m1_params, metric, &post_projection_callbacks, &post_projection_context,
              base_state, output, &diagnostics),
        ghl_error_m1_invalid_state, "Newton post-projection callback failure", 773);

  scripted_newton_context post_projection_success_context = { 0 };
  const ghl_m1_newton_callbacks post_projection_success_callbacks
        = { .residual = post_projection_success_residual,
            .jacobian = scripted_identity_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  require_error(
        checked_newton(
              m1_params, metric, &post_projection_success_callbacks,
              &post_projection_success_context, base_state, output, &diagnostics),
        ghl_success, "Newton successful projected trial", 774);
  require_condition(
        diagnostics.used_projection, "successful projected trial was not recorded", 774);

  scripted_newton_context overflow_context = { 0 };
  const ghl_m1_newton_callbacks overflow_callbacks
        = { .residual = overflow_trial_residual,
            .jacobian = scripted_identity_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  const double overflowing_base_state[4] = { DBL_MAX, 0.0, 0.0, 0.0 };
  require_error(
        checked_newton(
              m1_params, metric, &overflow_callbacks, &overflow_context,
              overflowing_base_state, output, &diagnostics),
        ghl_error_m1_implicit_solve_failure, "Newton nonfinite projected trial retry",
        775);

  scripted_newton_context trial_error_context = { 0 };
  const ghl_m1_newton_callbacks trial_error_callbacks
        = { .residual = trial_callback_error_residual,
            .jacobian = scripted_identity_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  require_error(
        checked_newton(
              m1_params, metric, &trial_error_callbacks, &trial_error_context,
              base_state, output, &diagnostics),
        ghl_error_m1_invalid_state, "Newton trial residual callback failure", 776);

  scripted_newton_context small_step_error_context = { 0 };
  const ghl_m1_newton_callbacks small_step_error_callbacks
        = { .residual = trial_callback_error_residual,
            .jacobian = small_step_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  const double unchanged_output[4] = { -1.0, -2.0, -3.0, -4.0 };
  memcpy(output, unchanged_output, sizeof(output));
  require_error(
        checked_newton(
              m1_params, metric, &small_step_error_callbacks, &small_step_error_context,
              base_state, output, &diagnostics),
        ghl_error_m1_invalid_state, "Newton small-step residual callback failure", 777);
  require_condition(
        small_step_error_context.calls == 2
              && memcmp(output, unchanged_output, sizeof(output)) == 0,
        "Newton small-step callback failure changed the output", 777);

  /* This residual is at least 3/4 for every E, despite its tiny first
   * Newton correction. A correction-only rule would falsely report success. */
  const ghl_m1_newton_callbacks no_root_callbacks
        = { .residual = no_root_small_step_residual,
            .jacobian = no_root_small_step_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  memcpy(output, unchanged_output, sizeof(output));
  require_error(
        checked_newton(
              &one_iteration_params, metric, &no_root_callbacks, NULL, base_state,
              output, &diagnostics),
        ghl_error_m1_implicit_solve_failure, "Newton small-step no-root residual", 778);
  require_condition(
        memcmp(output, unchanged_output, sizeof(output)) == 0
              && diagnostics.residual_weighted_merit > 1.0,
        "Newton small-step no-root case published false convergence", 778);

  /* The linear residual has its root below the energy floor. Its full Newton
   * correction is small, but the trial is inadmissible and projection cannot
   * lower the residual at the floor. */
  ghl_m1_parameters floor_params = *m1_params;
  floor_params.E_floor = 1.0;
  const ghl_m1_newton_callbacks floor_callbacks
        = { .residual = floor_limited_small_step_residual,
            .jacobian = small_step_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  memcpy(output, unchanged_output, sizeof(output));
  require_error(
        checked_newton(
              &floor_params, metric, &floor_callbacks, NULL, base_state, output,
              &diagnostics),
        ghl_error_m1_implicit_solve_failure, "Newton inadmissible small trial", 779);
  require_condition(
        memcmp(output, unchanged_output, sizeof(output)) == 0
              && diagnostics.residual_weighted_merit > 1.0,
        "Newton inadmissible small trial published false convergence", 779);

  /* A callback may succeed without checking the physical state domain. The
   * driver must not publish a zero-merit iterate below E_floor or one whose
   * components are nonfinite. */
  const ghl_m1_newton_callbacks zero_residual_callbacks
        = { .residual = zero_residual_for_any_state,
            .jacobian = scripted_identity_jacobian,
            .observer = NULL,
            .observer_context = NULL };
  const double valid_base[4] = { 1.0, 0.0, 0.0, 0.0 };
  const double invalid_initials[2][4]
        = { { 0.0, 0.0, 0.0, 0.0 }, { NAN, 0.0, 0.0, 0.0 } };
  const char *const invalid_initial_names[2] = { "Newton zero-energy converged iterate",
                                                 "Newton nonfinite converged iterate" };
  for(int invalid_case = 0; invalid_case < 2; ++invalid_case) {
    scripted_newton_context invalid_iterate_context = { 0 };
    memcpy(output, unchanged_output, sizeof(output));
    require_error(
          checked_newton_with_initial_guess(
                m1_params, metric, &zero_residual_callbacks, &invalid_iterate_context,
                valid_base, invalid_initials[invalid_case], output, &diagnostics),
          ghl_error_m1_implicit_admissibility, invalid_initial_names[invalid_case],
          780 + invalid_case);
    require_condition(
          invalid_iterate_context.calls == 1
                && memcmp(output, unchanged_output, sizeof(output)) == 0,
          "Newton invalid converged iterate was published", 780 + invalid_case);
  }

  ghl_m1_parameters roundtrip_params = *m1_params;
  roundtrip_params.E_floor = 1000.0;
  ghl_metric_quantities roundtrip_metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 3.0, 0.0, 1.0, &roundtrip_metric);
  volatile double densitized_floor
        = roundtrip_params.E_floor * roundtrip_metric.sqrt_detgamma;
  const double roundtrip_energy = densitized_floor / roundtrip_metric.sqrt_detgamma;
  require_condition(
        roundtrip_energy < roundtrip_params.E_floor,
        "Newton floor round-trip fixture did not cross E_floor", 782);
  const double floor_base[4] = { densitized_floor, 0.0, 0.0, 0.0 };
  const double floor_initial[4] = { densitized_floor, 0.0, 0.0, 0.0 };
  scripted_newton_context roundtrip_context = { 0 };
  memcpy(output, unchanged_output, sizeof(output));
  require_error(
        checked_newton_with_initial_guess(
              &roundtrip_params, &roundtrip_metric, &zero_residual_callbacks,
              &roundtrip_context, floor_base, floor_initial, output, &diagnostics),
        ghl_success, "Newton densitized E_floor round-trip", 782);
  require_condition(
        memcmp(output, floor_initial, sizeof(output)) == 0
              && roundtrip_context.calls == 1,
        "Newton rejected or changed the densitized floor round-trip", 782);

  const double subfloor_state[4] = { nextafter(densitized_floor, 0.0), 0.0, 0.0, 0.0 };
  scripted_newton_context subfloor_context = { 0 };
  memcpy(output, unchanged_output, sizeof(output));
  require_error(
        checked_newton_with_initial_guess(
              &roundtrip_params, &roundtrip_metric, &zero_residual_callbacks,
              &subfloor_context, floor_base, subfloor_state, output, &diagnostics),
        ghl_error_m1_implicit_admissibility, "Newton genuinely subfloor iterate", 783);
  require_condition(
        memcmp(output, unchanged_output, sizeof(output)) == 0
              && subfloor_context.calls == 1,
        "Newton published the genuinely subfloor iterate", 783);
}

/* The public solver does not validate the configuration before iterating. An
 * invalid repair policy therefore first surfaces at the trial admissibility
 * check and must abort without publishing. */
static void test_newton_trial_configuration_failure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {
  ghl_m1_parameters invalid_params = *m1_params;
  invalid_params.repair_policy = (ghl_m1_repair_policy_t)12345;
  const affine_newton_context affine = { .target = { 2.0, 0.1, 0.0, 0.0 } };
  const ghl_m1_newton_callbacks callbacks = { .residual = affine_newton_residual,
                                              .jacobian = affine_newton_jacobian,
                                              .observer = NULL,
                                              .observer_context = NULL };
  const double U_base[4] = { 1.0, 0.0, 0.0, 0.0 };
  const double unchanged[4] = { -1.0, -2.0, -3.0, -4.0 };
  double U_out[4] = { unchanged[0], unchanged[1], unchanged[2], unchanged[3] };
  ghl_m1_newton_diagnostics diagnostics;
  require_error(
        checked_newton(
              &invalid_params, metric, &callbacks, &affine, U_base, U_out, &diagnostics),
        ghl_error_m1_invalid_repair_policy, "Newton trial configuration failure", 776);
  for(int component = 0; component < 4; ++component) {
    require_close(
          U_out[component], unchanged[component], 0.0, 0.0,
          "Newton trial configuration failure transaction", 776);
  }
}

static void test_newton_extreme_admissibility(
      const ghl_m1_parameters *params, const ghl_metric_quantities *metric) {
  const affine_newton_context affine = { .target = {1.0, 2.0, 0.0, 0.0} };
  const ghl_m1_newton_callbacks affine_callbacks = {
    .residual = affine_newton_residual, .jacobian = affine_newton_jacobian
  };
  const double base[4] = {1.0, 0.0, 0.0, 0.0};
  double output[4] = {-1.0, -2.0, -3.0, -4.0};
  const double unchanged[4] = {-1.0, -2.0, -3.0, -4.0};
  require_error(checked_newton(params, metric, &affine_callbacks, &affine,
                  base, output, NULL),
                ghl_error_m1_implicit_solve_failure,
                "unrealizable affine target", 781);
  require_condition(memcmp(output, unchanged, sizeof(output)) == 0,
                    "unrealizable Newton target published output", 781);

  /* Densitizing a positive floor can lose more than one ulp when the product
   * is subnormal. Such a round trip must not be accepted as a floor tie. */
  const double volumes[] = {0.5, 1.0e-20};
  const double floors[] = {0x1p-1074, 1.0e-300};
  for(size_t i = 0; i < sizeof(volumes)/sizeof(volumes[0]); ++i) {
    ghl_metric_quantities small_metric = *metric;
    small_metric.sqrt_detgamma = volumes[i];
    small_metric.detgamma = volumes[i]*volumes[i];
    small_metric.gammaDD[0][0] = small_metric.detgamma;
    small_metric.gammaUU[0][0] = 1.0/small_metric.detgamma;
    ghl_m1_parameters small_params = *params;
    small_params.E_floor = floors[i];
    const double initial[4] = {floors[i]*volumes[i], 0.0, 0.0, 0.0};
    scripted_newton_context context = {0};
    const ghl_m1_newton_callbacks callbacks = {
      .residual = zero_residual_for_any_state, .jacobian = scripted_identity_jacobian
    };
    require_error(checked_newton(&small_params, &small_metric, &callbacks, &context,
                                initial, output, NULL),
                  ghl_error_m1_implicit_admissibility,
                  "subnormal densitized floor round trip", 782);
    require_condition(memcmp(output, unchanged, sizeof(output)) == 0,
                      "subnormal floor rejection published output", 782);
  }

  ghl_m1_parameters huge_floor = *params;
  huge_floor.E_floor = DBL_MAX;
  ghl_metric_quantities volume_metric = *metric;
  volume_metric.gammaDD[0][0] = volume_metric.detgamma = 4.0;
  volume_metric.gammaUU[0][0] = 0.25;
  volume_metric.sqrt_detgamma = 2.0;
  fixed_newton_context fixed = {
    .residual = {1.0, 0.0, 0.0, 0.0}, .residual_status = ghl_success,
    .jacobian_status = ghl_success
  };
  set_identity_matrix(fixed.jacobian);
  expect_solver_failure(&huge_floor, &volume_metric, &fixed,
                        "overflowing densitized Newton floor", 783);
}

static void test_newton_coverage_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {
  test_newton_extreme_admissibility(m1_params, metric);
  test_newton_trial_configuration_failure(m1_params, metric);
  test_weighted_merit_boundaries(m1_params, metric);
  test_projection_failure_boundaries(m1_params, metric);
  test_linear_rejection_boundaries(m1_params, metric);
  test_newton_input_boundaries(m1_params, metric);
  test_newton_optional_paths(m1_params, metric);
  test_newton_retry_boundaries(m1_params, metric);
}

static void test_public_newton_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric) {
  const affine_newton_context affine = { .target = { 1.25, 0.10, -0.20, 0.05 } };
  const ghl_m1_newton_callbacks callbacks = { .residual = affine_newton_residual,
                                              .jacobian = affine_newton_jacobian,
                                              .observer = NULL,
                                              .observer_context = NULL };
  const double U_base[4] = { 1.0, 0.0, 0.0, 0.0 };
  double U_out[4] = { -1.0, -2.0, -3.0, -4.0 };
  ghl_m1_newton_diagnostics diagnostics;
  require_error(
        checked_newton(
              m1_params, metric, &callbacks, &affine, U_base, U_out, &diagnostics),
        ghl_success, "public four-dimensional Newton solve", 700);
  for(int component = 0; component < 4; ++component) {
    require_close(
          U_out[component], affine.target[component], 0.0, 0.0, "public Newton solution",
          700 + component);
  }
  require_condition(
        diagnostics.iterations > 0 && diagnostics.residual_weighted_merit <= 1.0,
        "public Newton diagnostics did not report convergence", 705);

  double projected[4] = { 1.0, 2.0, 0.0, 0.0 };
  require_error(
        ghl_m1_newton_project_admissible(m1_params, metric, projected), ghl_success,
        "public admissibility projection", 706);
  require_condition(
        isfinite(projected[0]) && isfinite(projected[1]) && isfinite(projected[2])
              && isfinite(projected[3]) && fabs(projected[1]) < projected[0],
        "public admissibility projection published an invalid state", 706);

  const double unchanged[4] = { 1.0, 2.0, 3.0, 4.0 };
  double null_metric_projection[4]
        = { unchanged[0], unchanged[1], unchanged[2], unchanged[3] };
  require_error(
        ghl_m1_newton_project_admissible(m1_params, NULL, null_metric_projection),
        ghl_error_m1_null_pointer, "NULL-metric admissibility projection", 707);
  for(int component = 0; component < 4; ++component) {
    require_close(
          null_metric_projection[component], unchanged[component], 0.0, 0.0,
          "NULL-metric projection transaction", 707);
  }

  double null_parameters_projection[4]
        = { unchanged[0], unchanged[1], unchanged[2], unchanged[3] };
  require_error(
        ghl_m1_newton_project_admissible(NULL, metric, null_parameters_projection),
        ghl_error_m1_null_pointer, "NULL-parameter admissibility projection", 708);
  for(int component = 0; component < 4; ++component) {
    require_close(
          null_parameters_projection[component], unchanged[component], 0.0, 0.0,
          "NULL-parameter projection transaction", 708);
  }

  require_error(
        ghl_m1_newton_project_admissible(m1_params, metric, NULL),
        ghl_error_m1_null_pointer, "NULL-output admissibility projection", 709);

  ghl_metric_quantities invalid_metric = *metric;
  invalid_metric.sqrt_detgamma = 0.0;
  projected[0] = unchanged[0];
  projected[1] = unchanged[1];
  projected[2] = unchanged[2];
  projected[3] = unchanged[3];
  require_error(
        ghl_m1_newton_project_admissible(m1_params, &invalid_metric, projected),
        ghl_error_m1_invalid_metric, "invalid-metric admissibility projection", 710);
  for(int component = 0; component < 4; ++component) {
    require_close(
          projected[component], unchanged[component], 0.0, 0.0,
          "transactional admissibility projection", 710);
  }
}

int main(void) {
  ghl_m1_parameters m1_params;
  ghl_error_codes_t error = ghl_m1_initialize_with_newton_tolerances(
        1.0e-10, 1.0e-30, 1.0e-6, 1.0e-6, 1.0e-8, 40, 1.0e-10, 1.0e-12, &m1_params);
  if(error != ghl_success) {
    ghl_error("ghl_m1_initialize_with_newton_tolerances returned %d\n", (int)error);
  }

  ghl_metric_quantities metric;
  m1_setup_flat_metric(&metric);
  test_public_newton_boundaries(&m1_params, &metric);
  test_newton_coverage_boundaries(&m1_params, &metric);
  /* Repeat independent scenarios to establish the observer's success and
   * failure notification contract without changing solver outcomes. */
  check_observer_contract = true;
  test_public_newton_boundaries(&m1_params, &metric);
  test_newton_coverage_boundaries(&m1_params, &metric);
  check_observer_contract = false;
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(1.0, &rates);

  ghl_m1_neutrino_state partial_zero_state;
  make_trigger_state(1.0, false, &partial_zero_state);
  double partial_U[4];
  make_U(&metric, &partial_zero_state, partial_U);
  double partial_residual[4];
  compute_residual(
        &m1_params, &nu_params, &metric, &prims, &rates, &partial_zero_state, 1.0,
        partial_U, partial_residual, 100);
  require_close(
        partial_residual[2], -0.535304066609, 2.0e-10, 2.0e-12,
        "partial-zero transverse residual", 100);
  check_jacobian_against_independent_difference(
        &m1_params, &nu_params, &metric, &prims, &rates, &partial_zero_state, 1.0, 2,
        false, 101);
  check_centered_full_fy_column(
        &m1_params, &nu_params, &metric, &prims, &rates, &partial_zero_state);
  check_validated_context_equivalence(
        &m1_params, &nu_params, &metric, &prims, &rates, &partial_zero_state);
  check_invalid_fd_step_rejection(
        &m1_params, &nu_params, &metric, &prims, &rates, &partial_zero_state);

  ghl_m1_neutrino_state full_zero_state;
  make_trigger_state(1.0, true, &full_zero_state);
  check_full_zero_flux_columns(
        &m1_params, &nu_params, &metric, &prims, &rates, &full_zero_state, 1.0);

  /* Keep this explicit reference check separate so a zero column cannot be
   * hidden by a common-rescaling assertion. */
  double reference_U[4];
  make_U(&metric, &partial_zero_state, reference_U);
  double reference_residual[4];
  compute_residual(
        &m1_params, &nu_params, &metric, &prims, &rates, &partial_zero_state, 1.0,
        reference_U, reference_residual, 302);
  double reference_jacobian[4][4] = { { 0.0 } };
  compute_jacobian(
        &m1_params, &nu_params, &metric, &prims, &rates, &partial_zero_state, 1.0,
        reference_U, reference_residual, reference_jacobian, 302);

  const double scales[3] = { 1.0e-4, 1.0, 1.0e4 };
  for(int scale_index = 0; scale_index < 3; ++scale_index) {
    ghl_m1_neutrino_rates scaled_rates;
    make_rates(scales[scale_index], &scaled_rates);
    ghl_m1_neutrino_state scaled_state;
    make_trigger_state(scales[scale_index], false, &scaled_state);
    check_scaled_case(
          &m1_params, &nu_params, &metric, &prims, &scaled_rates, &scaled_state, 1.0,
          reference_jacobian[2][2]);
  }

  ghl_m1_neutrino_state boundary_state = partial_zero_state;
  boundary_state.F[0] = sqrt(m1_params.one_minus_epsilon_c_sq) - 1.0e-11;
  check_jacobian_against_independent_difference(
        &m1_params, &nu_params, &metric, &prims, &rates, &boundary_state, 1.0, 1, true,
        500);

  check_public_convergence(
        &m1_params, &nu_params, &metric, &prims, &rates, &partial_zero_state);

  ghl_metric_quantities nonunit_metric;
  m1_setup_flat_metric(&nonunit_metric);
  nonunit_metric.gammaDD[0][0] = 4.0;
  nonunit_metric.gammaUU[0][0] = 0.25;
  nonunit_metric.detgamma = 4.0;
  nonunit_metric.sqrt_detgamma = 2.0;
  check_jacobian_against_independent_difference(
        &m1_params, &nu_params, &nonunit_metric, &prims, &rates, &partial_zero_state,
        1.0, 2, false, 600);
  ghl_m1_neutrino_state nonunit_full_zero_state;
  make_trigger_state(1.0, true, &nonunit_full_zero_state);
  check_public_convergence(
        &m1_params, &nu_params, &nonunit_metric, &prims, &rates,
        &nonunit_full_zero_state);

  ghl_info("M1 finite-difference Jacobian regression passed\n");
  return 0;
}
