#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"
#include "../ghl_m1_utils.h"

/*
 * Local homogeneous source solve for one neutrino species. It solves
 * E/F_i first, derives the endpoint Gamma_N from the fresh shared current,
 * then updates N and constructs one signed charged-current exchange packet in
 * temporary storage.
 *
 * The solve follows the shared Newton/line-search/substep pattern but uses
 * frozen rates and primitives instead of Con2Prim or opacity callbacks. After
 * required pointers are validated, a failure leaves state_in and a zero packet
 * published; terminal recovery never publishes a partial source update.
 */

static ghl_error_codes_t ghl_m1_neutrino_publish_hard_failure(
      const ghl_error_codes_t error,
      ghl_m1_neutrino_diagnostics *restrict diagnostics) {

  diagnostics->source_failures++;
  return error;
}

ghl_error_codes_t ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double n_b_cons,
      const ghl_m1_neutrino_state *restrict state_in,
      int *restrict thin_inequalities_hold,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      ghl_m1_neutrino_diagnostics *restrict diagnostics) {

  if(m1_params == NULL || nu_params == NULL || metric == NULL || prims == NULL ||
     rates == NULL || state_in == NULL || thin_inequalities_hold == NULL ||
     state_out == NULL || exchange == NULL)
    return ghl_error_m1_null_pointer;

  *state_out = *state_in;
  *exchange = (ghl_m1_neutrino_exchange){ 0 };
  *thin_inequalities_hold = 0;

  ghl_m1_neutrino_diagnostics candidate_diagnostics;
  ghl_m1_neutrino_diagnostics *candidate_diagnostics_ptr = NULL;
  if(diagnostics != NULL) {
    candidate_diagnostics = *diagnostics;
    candidate_diagnostics_ptr = &candidate_diagnostics;
  }

  if(!isfinite(dt) || dt < 0.0 || !isfinite(n_b_cons) || n_b_cons <= 0.0)
    return ghl_error_m1_invalid_state;
  const ghl_error_codes_t rates_error =
        ghl_m1_neutrino_validate_single_species_rates(rates, NULL);
  if(rates_error != ghl_success)
    return rates_error;

  const double dt_alpha = metric->lapse * dt;
  if(!isfinite(dt_alpha) || dt_alpha < 0.0)
    return ghl_error_m1_invalid_state;

  /* This is the two-inequality four-point reference branch selector. A thick packet is a
   * successful non-selection so the validation host can retain its implicit
   * source solve without changing the default GRHayL algorithm. */
  if(!(dt_alpha * rates->kappa_a_E < 1.0) ||
     !(dt_alpha * rates->kappa_s < 1.0))
    return ghl_success;
  *thin_inequalities_hold = 1;

  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state_in);
  ghl_m1_closure closure;
  ghl_error_codes_t error = ghl_m1_compute_closure_with_primitives(
        m1_params, metric, prims, &rad_state, &closure);
  if(error != ghl_success)
    return error;

  ghl_m1_sources EF_sources = { 0 };
  double ignored_N_source = 0.0;
  error = ghl_m1_compute_neutrino_interaction_sources_from_closure(
        m1_params, nu_params, metric, prims, state_in, &closure, rates,
        &EF_sources, &ignored_N_source);
  if(error != ghl_success)
    return error;

  ghl_m1_neutrino_state candidate = *state_in;
  candidate.E += dt_alpha * EF_sources.S_E;
  for(int direction = 0; direction < 3; ++direction)
    candidate.F[direction] += dt_alpha * EF_sources.S[direction];

  /* The explicit source increment can move a near-streaming state a few
   * ulps outside the realizability cone even when the input was repaired.
   * Repair the complete candidate transactionally so the diagnostics-aware
   * API accounts for both E/F and N mutations. */
  error = ghl_m1_repair_neutrino_state(
      m1_params, nu_params, metric, &candidate, candidate_diagnostics_ptr);
  if(error != ghl_success)
    return error;

  /* The implicit source route updates N after the explicit E/F endpoint.
   * Gamma_N therefore belongs to the repaired endpoint current, and the
   * number update is backward Euler with the same frozen rates used for E/F. */
  ghl_m1_neutrino_current endpoint_current;
  error = ghl_m1_neutrino_derive_current(
        m1_params, nu_params, metric, prims, &candidate, &endpoint_current);
  if(error != ghl_success)
    return error;

  error = ghl_m1_update_neutrino_number_backward_euler(
        nu_params, rates, dt_alpha, endpoint_current.Gamma_N,
        state_in->N, &candidate.N);
  if(error != ghl_success)
    return error;

  /* Apply the configured number floor before constructing the exchange packet.
   * E/F are already repaired above; this second repair is consequently an
   * N-only operation for the successful thin branch. */
  error = ghl_m1_repair_neutrino_state(
        m1_params, nu_params, metric, &candidate, candidate_diagnostics_ptr);
  if(error != ghl_success)
    return error;

  error = ghl_m1_neutrino_derive_current(
      m1_params, nu_params, metric, prims, &candidate, &endpoint_current);
  if(error != ghl_success)
    return error;
  error = ghl_m1_neutrino_check_EN_bounds(
      &candidate, nu_params, &endpoint_current);
  if(error != ghl_success)
    return error;

  /* The charged-current exchange packet remains endpoint-based, as required
   * by the public GRHayL exchange contract. */
  const double dN_cc = dt_alpha *
        (rates->eta_N_cc - rates->kappa_a_N_cc * candidate.N /
         endpoint_current.Gamma_N);
  const double dL_rad_cc = rates->lepton_weight * dN_cc;
  error = ghl_m1_neutrino_assemble_exchange(
        state_in, &candidate, rates, dL_rad_cc, metric->sqrt_detgamma,
        n_b_cons, exchange);
  if(error != ghl_success)
    return error;

  *state_out = candidate;
  if(diagnostics != NULL)
    *diagnostics = candidate_diagnostics;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_try_neutrino_explicit_thin_update(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double n_b_cons,
      const ghl_m1_neutrino_state *restrict state_in,
      int *restrict thin_inequalities_hold,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange) {
  return ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
      m1_params, nu_params, metric, prims, rates, dt, n_b_cons, state_in,
      thin_inequalities_hold, state_out, exchange, NULL);
}

typedef struct {
  ghl_m1_neutrino_implicit_context validated;
  double dt_sub;
  double U_base[4];
  bool closure_fallback_observed;
} ghl_m1_neutrino_newton_context;

static ghl_error_codes_t ghl_m1_neutrino_newton_residual(
      const void *restrict context,
      const double U[4],
      double residual[4]) {

  ghl_m1_neutrino_newton_context *restrict solve_context
        = (ghl_m1_neutrino_newton_context *)context;
  return ghl_m1_neutrino_compute_implicit_residual_validated(
        &solve_context->validated, solve_context->dt_sub, solve_context->U_base, U,
        &solve_context->closure_fallback_observed, residual);
}

static ghl_error_codes_t ghl_m1_neutrino_newton_jacobian(
      const void *restrict context,
      const double U[4],
      const double residual[4],
      double jacobian[4][4]) {

  const ghl_m1_neutrino_newton_context *restrict solve_context = context;
  return ghl_m1_neutrino_compute_implicit_jacobian_validated(
        &solve_context->validated, solve_context->dt_sub, solve_context->U_base,
        U, residual, jacobian);
}

static ghl_error_codes_t ghl_m1_neutrino_build_EF_initial_guess(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt_sub,
      const double U_in[4],
      double U_initial[4],
      bool *restrict closure_fallback_observed) {
  if(m1_params == NULL || metric == NULL || prims_frozen == NULL || rates == NULL
     || U_in == NULL || U_initial == NULL || closure_fallback_observed == NULL) {
    return ghl_error_m1_null_pointer;
  }
  *closure_fallback_observed = false;
  for(int i = 0; i < 4; ++i) {
    U_initial[i] = U_in[i];
  }

  if(!isfinite(dt_sub) || dt_sub <= 0.0 || U_in[1] != 0.0 || U_in[2] != 0.0
     || U_in[3] != 0.0) {
    return ghl_success;
  }

  const double dt_alpha_sqrt_detgamma = metric->lapse * dt_sub * metric->sqrt_detgamma;
  if(!isfinite(dt_alpha_sqrt_detgamma) || dt_alpha_sqrt_detgamma < 0.0) {
    return ghl_error_m1_invalid_state;
  }

  ghl_m1_rad_state rad_state;
  ghl_error_codes_t error = ghl_m1_neutrino_build_trial_state(metric, U_in, &rad_state);
  if(error != ghl_success) {
    return error;
  }
  bool source_closure_fallback_observed = false;
  ghl_m1_sources sources;
  error = ghl_m1_neutrino_compute_EF_interaction_sources_validated(
        m1_params, metric, prims_frozen, &rad_state, rates,
        &source_closure_fallback_observed, &sources);
  if(error != ghl_success) {
    return error;
  }
  *closure_fallback_observed = source_closure_fallback_observed;
  U_initial[0] += dt_alpha_sqrt_detgamma * sources.S_E;
  for(int i = 0; i < 3; ++i) {
    U_initial[i + 1] += dt_alpha_sqrt_detgamma * sources.S[i];
  }
  for(int i = 0; i < 4; ++i) {
    if(!isfinite(U_initial[i])) {
      return ghl_error_m1_invalid_state;
    }
  }

  ghl_m1_rad_state predicted_state;
  error = ghl_m1_neutrino_build_trial_state(metric, U_initial, &predicted_state);
  if(error != ghl_success
     || ghl_m1_neutrino_check_trial_admissibility(m1_params, metric, &predicted_state)
              != ghl_success) {
    for(int i = 0; i < 4; ++i) {
      U_initial[i] = U_in[i];
    }
  }
  return ghl_success;
}

static ghl_error_codes_t ghl_m1_neutrino_attempt_newton_step_with_initial_guess(
      const ghl_m1_neutrino_implicit_context *restrict validated_context,
      const double dt_sub,
      const double U_in[4],
      const double U_initial[4],
      double U_out[4],
      ghl_m1_newton_diagnostics *restrict diagnostics,
      bool *restrict closure_fallback_observed) {
  ghl_m1_neutrino_newton_context solve_context
        = { .validated = *validated_context,
            .dt_sub = dt_sub,
            .U_base = { U_in[0], U_in[1], U_in[2], U_in[3] },
            .closure_fallback_observed = false };
  const ghl_m1_newton_callbacks callbacks
        = { .residual = ghl_m1_neutrino_newton_residual,
            .jacobian = ghl_m1_neutrino_newton_jacobian };
  const ghl_error_codes_t error = ghl_m1_newton_solve_4d_with_initial_guess(
        validated_context->m1_params, validated_context->metric, &callbacks,
        &solve_context, U_in, U_initial, U_out, diagnostics);

  *closure_fallback_observed = solve_context.closure_fallback_observed;
  return error;
}

static ghl_error_codes_t ghl_m1_neutrino_attempt_newton_step(
      const ghl_m1_neutrino_implicit_context *restrict validated_context,
      const double dt_sub,
      const double U_in[4],
      double U_out[4],
      ghl_m1_newton_diagnostics *restrict diagnostics,
      bool *restrict closure_fallback_observed) {
  double U_initial[4];
  bool predictor_fallback_observed = false;
  const ghl_error_codes_t guess_error = ghl_m1_neutrino_build_EF_initial_guess(
        validated_context->m1_params, validated_context->metric,
        validated_context->prims_frozen, validated_context->rates, dt_sub, U_in,
        U_initial, &predictor_fallback_observed);
  if(guess_error != ghl_success) {
    for(int i = 0; i < 4; ++i) {
      U_initial[i] = U_in[i];
    }
  }
  bool newton_fallback_observed = false;
  const ghl_error_codes_t error = ghl_m1_neutrino_attempt_newton_step_with_initial_guess(
        validated_context, dt_sub, U_in, U_initial, U_out, diagnostics,
        &newton_fallback_observed);
  *closure_fallback_observed = predictor_fallback_observed || newton_fallback_observed;
  return error;
}

ghl_error_codes_t ghl_m1_neutrino_attempt_EF_newton_step(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt_sub,
      const double U_in[4],
      double U_out[4],
      ghl_m1_newton_diagnostics *restrict diagnostics,
      bool *restrict closure_fallback_observed) {
  if(m1_params == NULL || metric == NULL || prims_frozen == NULL || rates == NULL
     || U_in == NULL || U_out == NULL || diagnostics == NULL
     || closure_fallback_observed == NULL) {
    return ghl_error_m1_null_pointer;
  }
  const ghl_m1_neutrino_implicit_context context = { .m1_params = m1_params,
                                                     .metric = metric,
                                                     .prims_frozen = prims_frozen,
                                                     .rates = rates };
  return ghl_m1_neutrino_attempt_newton_step(
        &context, dt_sub, U_in, U_out, diagnostics, closure_fallback_observed);
}

/*
 * Populate the comoving mean-energy diagnostic fields of
 * ghl_m1_neutrino_diagnostics. Called once per solve in both the
 * success and terminal-fallback paths.
 *
 * The comoving mean energy is J * Gamma_N / N. Consistency checks use
 * the fixed relative tolerance |a - b| / b < 1e-3. The diagnostic is
 * invalid when N is at or below the number floor, where this ratio
 * is not meaningful.
 */
void ghl_m1_neutrino_populate_mean_energy_diagnostics(
      const ghl_m1_neutrino_state *restrict state_out,
      const ghl_m1_neutrino_current *restrict current,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics) {

  if(neutrino_diagnostics == NULL || state_out == NULL || current == NULL
     || nu_params == NULL || rates == NULL) {
    return;
  }

  neutrino_diagnostics->mean_energy_diag = 0.0;
  neutrino_diagnostics->mean_energy_consistent = 0;
  neutrino_diagnostics->Jeq_over_neq_consistent = 0;
  neutrino_diagnostics->mean_energy_diag_invalid = 1;
  const double tol_meandiff = 1.0e-3;

  /* The post-update comoving grey mean energy is J Gamma_N / N. */
  if(!isfinite(state_out->N) || state_out->N <= nu_params->N_floor ||
     !isfinite(current->J) || !isfinite(current->Gamma_N))
    return;

  const double mean_energy_diag = current->J * current->Gamma_N / state_out->N;
  if(!isfinite(mean_energy_diag))
    return;
  neutrino_diagnostics->mean_energy_diag = mean_energy_diag;
  neutrino_diagnostics->mean_energy_diag_invalid = 0;

  if(!isfinite(rates->mean_energy) || rates->mean_energy <= 0.0)
    return;

  const double ref = rates->mean_energy;
  const double rel = fabs(mean_energy_diag - ref) / ref;
  neutrino_diagnostics->mean_energy_consistent = rel < tol_meandiff;

  if(rates->n_eq > 0.0 && rates->J_eq > 0.0) {
    const double Jeq_over_neq = rates->J_eq / rates->n_eq;
    if(isfinite(Jeq_over_neq)) {
      const double rel_eq = fabs(Jeq_over_neq - ref) / ref;
      neutrino_diagnostics->Jeq_over_neq_consistent = rel_eq < tol_meandiff;
    }
  }
}

ghl_error_codes_t ghl_m1_solve_neutrino_implicit_homogeneous_update(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double n_b_cons,
      const ghl_m1_neutrino_state *restrict state_in,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      ghl_m1_implicit_solve_diagnostics *restrict solve_diagnostics,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics) {

  if(m1_params == NULL || nu_params == NULL || metric == NULL || prims_frozen == NULL
     || rates == NULL || state_in == NULL || state_out == NULL || exchange == NULL
     || solve_diagnostics == NULL || neutrino_diagnostics == NULL) {
    return ghl_error_m1_null_pointer;
  }

  /* Once all output pointers are known, state and exchange are transactional.
   * Diagnostics are observability state: every later hard failure increments
   * source_failures, while exhausted retry schedules increment
   * source_terminal_fallbacks. */
  *state_out = *state_in;
  *exchange = (ghl_m1_neutrino_exchange){ 0 };

  ghl_m1_implicit_solve_diagnostics candidate_solve_diagnostics;
  ghl_m1_initialize_implicit_solve_diagnostics(&candidate_solve_diagnostics);
  ghl_m1_neutrino_diagnostics candidate_neutrino_diagnostics = *neutrino_diagnostics;

  /* Validate basic inputs. */
  if(!isfinite(dt) || dt < 0.0 || !isfinite(n_b_cons) || n_b_cons <= 0.0) {
    return ghl_m1_neutrino_publish_hard_failure(
          ghl_error_m1_invalid_state, neutrino_diagnostics);
  }
  if(nu_params->terminal_fallback_policy
     != ghl_m1_neutrino_terminal_fallback_no_update_all) {
    return ghl_m1_neutrino_publish_hard_failure(
          ghl_error_m1_invalid_state, neutrino_diagnostics);
  }

  if(!ghl_m1_metric_is_symmetric_spd(metric)) {
    return ghl_m1_neutrino_publish_hard_failure(
          ghl_error_m1_invalid_metric, neutrino_diagnostics);
  }

  const double dt_alpha = metric->lapse * dt;
  if(!isfinite(dt_alpha) || dt_alpha < 0.0) {
    return ghl_m1_neutrino_publish_hard_failure(
          ghl_error_m1_invalid_state, neutrino_diagnostics);
  }

  if(!isfinite(state_in->N) || !isfinite(state_in->E)) {
    return ghl_m1_neutrino_publish_hard_failure(
          ghl_error_m1_invalid_state, neutrino_diagnostics);
  }
  for(int i = 0; i < 3; i++) {
    if(!isfinite(state_in->F[i])) {
      return ghl_m1_neutrino_publish_hard_failure(
            ghl_error_m1_invalid_state, neutrino_diagnostics);
    }
  }

  /* Validate frozen rates. */
  ghl_error_codes_t rates_error = ghl_m1_neutrino_validate_single_species_rates(rates, NULL);
  if(rates_error != ghl_success) {
    return ghl_m1_neutrino_publish_hard_failure(rates_error, neutrino_diagnostics);
  }

  /* The public solver consumes a full N/E/F state. Validate that state with
   * the shared fresh-current path before entering Newton so a nonrepairable
   * input failure is deterministic and cannot be converted into a fallback
   * solve failure. The endpoint current is recomputed after E/F converges. */
  ghl_m1_neutrino_current input_current;
  const ghl_error_codes_t input_current_error = ghl_m1_neutrino_derive_current(
        m1_params, nu_params, metric, prims_frozen, state_in, &input_current);
  if(input_current_error != ghl_success) {
    return ghl_m1_neutrino_publish_hard_failure(
          input_current_error, neutrino_diagnostics);
  }

  const ghl_m1_neutrino_implicit_context validated_context = {
        .m1_params = m1_params,
        .metric = metric,
        .prims_frozen = prims_frozen,
        .rates = rates };

  const double sqrt_detgamma = metric->sqrt_detgamma;

  /* Densitized input unknowns U_in = (tildeE, tildeF_i). */
  double U_in[4] = { 0.0, 0.0, 0.0, 0.0 };
  U_in[0] = state_in->E * sqrt_detgamma;
  for(int i = 0; i < 3; i++) {
    U_in[i + 1] = state_in->F[i] * sqrt_detgamma;
  }

  const int fallback_schedule[] = { 1, 2, 4, 8, 16 };
  const int num_schedules
        = (int)(sizeof(fallback_schedule) / sizeof(fallback_schedule[0]));

  double U_final[4] = { U_in[0], U_in[1], U_in[2], U_in[3] };
  int total_iterations = 0;
  int total_backtracks = 0;
  double residual_norm = INFINITY;
  double residual_scaled_norm = INFINITY;
  bool any_projection = false;
  bool any_closure_fallback = false;
  bool any_backtracking = false;
  int successful_substeps = 0;

  for(int schedule_idx = 0; schedule_idx < num_schedules; schedule_idx++) {
    const int num_substeps = fallback_schedule[schedule_idx];
    const double dt_sub = dt / (double)num_substeps;

    double U_current[4] = { U_in[0], U_in[1], U_in[2], U_in[3] };
    int step_iterations = 0;
    int step_backtracks = 0;
    double step_residual_norm = INFINITY;
    double step_residual_scaled_norm = INFINITY;
    bool solve_failed = false;

    for(int step = 0; step < num_substeps; step++) {
      double U_next[4] = { 0.0, 0.0, 0.0, 0.0 };
      ghl_m1_newton_diagnostics step_diagnostics;
      bool closure_fallback_used = false;
      const ghl_error_codes_t error = ghl_m1_neutrino_attempt_newton_step(
            &validated_context, dt_sub, U_current, U_next, &step_diagnostics,
            &closure_fallback_used);
      any_backtracking = any_backtracking || step_diagnostics.backtracks > 0;
      any_projection = any_projection || step_diagnostics.used_projection;
      any_closure_fallback = any_closure_fallback || closure_fallback_used;
      if(error != ghl_success) {
        if(!ghl_m1_schedule_error_allows_retry(error))
          return ghl_m1_neutrino_publish_hard_failure(
              error, neutrino_diagnostics);
        solve_failed = true;
        break;
      }

      step_iterations += step_diagnostics.iterations;
      step_backtracks += step_diagnostics.backtracks;
      step_residual_norm = step_diagnostics.residual_max_norm;
      step_residual_scaled_norm = step_diagnostics.residual_weighted_merit;

      for(int i = 0; i < 4; i++) {
        U_current[i] = U_next[i];
      }
    }

    if(!solve_failed) {
      for(int i = 0; i < 4; i++) {
        U_final[i] = U_current[i];
      }
      total_iterations = step_iterations;
      total_backtracks = step_backtracks;
      residual_norm = step_residual_norm;
      residual_scaled_norm = step_residual_scaled_norm;
      successful_substeps = num_substeps;
      break;
    }
  }

  if(successful_substeps > 0) {
    ghl_m1_neutrino_state candidate = *state_in;
    ghl_m1_neutrino_exchange candidate_exchange = { 0 };
    /* Undensitize the converged E/F_i. */
    const double inv_sqrt_detgamma = 1.0 / sqrt_detgamma;
    candidate.E = U_final[0] * inv_sqrt_detgamma;
    for(int i = 0; i < 3; i++) {
      candidate.F[i] = U_final[i + 1] * inv_sqrt_detgamma;
    }

    /* Backward-Euler N update with dt_alpha = dt*lapse and the fresh endpoint
     * current normalization. Gamma_N is not a fluid Lorentz factor. */
    double N_out = 0.0;
    ghl_m1_neutrino_current endpoint_current;
    ghl_error_codes_t n_error = ghl_m1_neutrino_derive_current(
          m1_params, nu_params, metric, prims_frozen, &candidate, &endpoint_current);
    if(n_error != ghl_success) {
      return ghl_m1_neutrino_publish_hard_failure(n_error, neutrino_diagnostics);
    }
    n_error = ghl_m1_update_neutrino_number_backward_euler(
          nu_params, rates, dt_alpha, endpoint_current.Gamma_N, state_in->N, &N_out);
    if(n_error != ghl_success) {
      return ghl_m1_neutrino_publish_hard_failure(n_error, neutrino_diagnostics);
    }
    candidate.N = N_out;

    /* The Newton solve only evolves E/F_i.  Apply the same transactional
     * post-update repair used by the explicit compatibility branches so the
     * general implicit route cannot publish N below the configured floor. */
    n_error = ghl_m1_repair_neutrino_state(
          m1_params, nu_params, metric, &candidate,
          &candidate_neutrino_diagnostics);
    if(n_error != ghl_success) {
      return ghl_m1_neutrino_publish_hard_failure(n_error, neutrino_diagnostics);
    }

    /* Repair may have changed N (and, for a marginal E/F endpoint, the
     * realizability projection may change E/F as well).  All endpoint-based
     * diagnostics and charged-current exchange must use the published
     * candidate. */
    n_error = ghl_m1_neutrino_derive_current(
          m1_params, nu_params, metric, prims_frozen, &candidate,
          &endpoint_current);
    if(n_error != ghl_success) {
      return ghl_m1_neutrino_publish_hard_failure(n_error, neutrino_diagnostics);
    }

    /* The bound applies to the final repaired endpoint.  Check it before
     * constructing any exchange packet so a rejected endpoint publishes no
     * derived source data. */
    ghl_m1_neutrino_populate_mean_energy_diagnostics(
          &candidate, &endpoint_current, nu_params, rates,
          &candidate_neutrino_diagnostics);
    const ghl_error_codes_t en_error =
          ghl_m1_neutrino_check_EN_bounds(&candidate, nu_params, &endpoint_current);
    if(en_error != ghl_success) {
      return ghl_m1_neutrino_publish_hard_failure(en_error, neutrino_diagnostics);
    }

    /* Populate dYe_matter from the signed charged-current packet only. The
     * explicit n_b_cons argument fixes the host's baryon-number convention. */
    const double dN_cc
          = dt_alpha
            * (rates->eta_N_cc
               - rates->kappa_a_N_cc * candidate.N / endpoint_current.Gamma_N);
    const double dL_rad_cc = rates->lepton_weight * dN_cc;
    const ghl_error_codes_t ye_error = ghl_m1_neutrino_assemble_exchange(
          state_in, &candidate, rates, dL_rad_cc, sqrt_detgamma, n_b_cons,
          &candidate_exchange);
    if(ye_error != ghl_success) {
      return ghl_m1_neutrino_publish_hard_failure(ye_error, neutrino_diagnostics);
    }

    /* Update diagnostics. */
    candidate_neutrino_diagnostics.source_converged++;

    candidate_solve_diagnostics.newton_iterations = total_iterations;
    candidate_solve_diagnostics.line_search_backtracks = total_backtracks;
    candidate_solve_diagnostics.fallback_substeps = successful_substeps;
    candidate_solve_diagnostics.used_fallback_substepping = successful_substeps > 1;
    candidate_solve_diagnostics.residual_max_norm = residual_norm;
    candidate_solve_diagnostics.residual_scaled_norm = residual_scaled_norm;
    candidate_solve_diagnostics.solution_path_flags
          = (successful_substeps == 1 ? ghl_m1_solution_path_primary_convergence
                                      : ghl_m1_solution_path_substepping)
            | (any_backtracking ? ghl_m1_solution_path_line_search_backtracking : 0u)
            | (any_projection ? ghl_m1_solution_path_projection : 0u)
            | (any_closure_fallback ? ghl_m1_solution_path_closure_fallback : 0u)
            | ghl_m1_solution_path_endpoint_acceptance;

    *state_out = candidate;
    *exchange = candidate_exchange;
    *solve_diagnostics = candidate_solve_diagnostics;
    *neutrino_diagnostics = candidate_neutrino_diagnostics;
    return ghl_success;
  }

  candidate_neutrino_diagnostics.source_terminal_fallbacks++;
  ghl_m1_neutrino_populate_mean_energy_diagnostics(
        state_in, &input_current, nu_params, rates, &candidate_neutrino_diagnostics);
  candidate_solve_diagnostics.fallback_substeps = fallback_schedule[num_schedules - 1];
  candidate_solve_diagnostics.used_fallback_substepping = true;
  candidate_solve_diagnostics.solution_path_flags
        = ghl_m1_solution_path_substepping | ghl_m1_solution_path_terminal_failure
          | (any_backtracking ? ghl_m1_solution_path_line_search_backtracking : 0u)
          | (any_projection ? ghl_m1_solution_path_projection : 0u)
          | (any_closure_fallback ? ghl_m1_solution_path_closure_fallback : 0u);
  *solve_diagnostics = candidate_solve_diagnostics;
  *neutrino_diagnostics = candidate_neutrino_diagnostics;
  return ghl_error_m1_implicit_terminal_fallback;
}
