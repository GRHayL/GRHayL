#include "../ghl_m1_utils.h"
#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"

#include <math.h>

/*
 * Coupled grey electron-neutrino pair source.
 *
 * The number reaction is solved analytically as one scalar backward-Euler
 * extent shared by nu_e and anti-nu_e.  E/F is then solved for both species
 * with the same substep schedule as the ordinary implicit source route.  The
 * effective E/F rates are private to this file: the public rate validator is
 * intentionally applied only to the provider bundle, before this spectral
 * projection is constructed.
 */

enum { ghl_m1_pair_species_count = 2 };

static void ghl_m1_pair_zero_source_diagnostics(
      ghl_m1_neutrino_source_diagnostics *restrict diagnostics) {
  *diagnostics = (ghl_m1_neutrino_source_diagnostics){ 0 };
  ghl_m1_initialize_implicit_solve_diagnostics(&diagnostics->implicit);
}

static bool ghl_m1_pair_state_is_finite(const ghl_m1_neutrino_state *restrict state) {
  if(state == NULL || !isfinite(state->N) || !isfinite(state->E)) {
    return false;
  }
  for(int i = 0; i < 3; ++i) {
    if(!isfinite(state->F[i])) {
      return false;
    }
  }
  return true;
}

static bool ghl_m1_pair_rates_are_active(
      const ghl_m1_neutrino_rates rates[ghl_m1_pair_species_count]) {
  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
      if(rates[species].eta_N_pair[process] != 0.0
         || rates[species].eta_E_pair[process] != 0.0) {
        return true;
      }
    }
  }
  return false;
}

static void ghl_m1_pair_zero_fields(ghl_m1_neutrino_rates *restrict rates) {
  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    rates->eta_N_pair[process] = 0.0;
    rates->eta_E_pair[process] = 0.0;
  }
}

/* Solve the shared extent equation by evolving the smaller endpoint number.
 * Forming N+d directly loses all digits when a large population undergoes a
 * nearly complete inverse reaction.  With m=min(N_e,N_a), M=max(...),
 * delta=M-m, H=h*q, and D=Gamma_e*Gamma_a*n_eq_e*n_eq_a, the smaller new
 * population y satisfies
 *
 *   y^2 + (delta + D/H)y - D(1 + m/H) = 0.
 *
 * The positive root is evaluated in scaled form and the original difference
 * is restored exactly in the other endpoint. */
static ghl_error_codes_t ghl_m1_pair_number_extent(
      const ghl_m1_neutrino_rates rates[ghl_m1_pair_species_count],
      const ghl_m1_neutrino_current current[ghl_m1_pair_species_count],
      const ghl_m1_neutrino_state state[ghl_m1_pair_species_count],
      const double h,
      double new_number[ghl_m1_pair_species_count]) {
  if(new_number == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(h) || h < 0.0) {
    return ghl_error_m1_invalid_state;
  }

  double q = 0.0;
  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    if(rates[0].eta_N_pair[process] != rates[1].eta_N_pair[process]) {
      return ghl_error_m1_microphysics_failure;
    }
    q += rates[0].eta_N_pair[process];
  }
  if(!isfinite(q) || q < 0.0) {
    return ghl_error_m1_microphysics_failure;
  }
  if(q == 0.0 || h == 0.0) {
    new_number[0] = state[0].N;
    new_number[1] = state[1].N;
    return ghl_success;
  }

  const double gamma_e = current[0].Gamma_N;
  const double gamma_a = current[1].Gamma_N;
  const double n_eq_e = rates[0].n_eq;
  const double n_eq_a = rates[1].n_eq;
  const double N_e = state[0].N;
  const double N_a = state[1].N;
  if(!isfinite(gamma_e) || !isfinite(gamma_a) || gamma_e <= 0.0 || gamma_a <= 0.0
     || !isfinite(n_eq_e) || !isfinite(n_eq_a) || n_eq_e <= 0.0 || n_eq_a <= 0.0
     || !isfinite(N_e) || !isfinite(N_a)) {
    return ghl_error_m1_invalid_state;
  }

  if(N_e < 0.0 || N_a < 0.0) {
    return ghl_error_m1_invalid_state;
  }

  const long double H = (long double)h * (long double)q;
  const long double D = (long double)gamma_e * (long double)gamma_a * (long double)n_eq_e
                        * (long double)n_eq_a;
  const long double smaller = N_e < N_a ? (long double)N_e : (long double)N_a;
  const long double difference = fabsl((long double)N_e - (long double)N_a);
  if(!isfinite(H) || H <= 0.0L || !isfinite(D) || D <= 0.0L || !isfinite(smaller)
     || !isfinite(difference)) {
    return ghl_error_m1_implicit_admissibility;
  }

  const long double D_over_H = D / H;
  const long double C = D * (1.0L + smaller / H);
  const long double B = difference + D_over_H;
  if(!isfinite(D_over_H) || !isfinite(B) || !isfinite(C) || B < 0.0L || C <= 0.0L) {
    return ghl_error_m1_implicit_admissibility;
  }
  const long double scale = fmaxl(B, sqrtl(C));
  if(!isfinite(scale) || scale <= 0.0L) {
    return ghl_error_m1_implicit_admissibility;
  }
  const long double b_scaled = B / scale;
  const long double c_over_scale = C / scale;
  const long double discriminant_scaled
        = b_scaled * b_scaled + 4.0L * c_over_scale / scale;
  if(!isfinite(discriminant_scaled) || discriminant_scaled < 0.0L) {
    return ghl_error_m1_implicit_admissibility;
  }
  const long double y = 2.0L * c_over_scale / (b_scaled + sqrtl(discriminant_scaled));
  if(!isfinite(y) || y < 0.0L || !isfinite((double)y)) {
    return ghl_error_m1_implicit_admissibility;
  }
  const long double larger = y + difference;
  if(!isfinite(larger) || !isfinite((double)larger)) {
    return ghl_error_m1_implicit_admissibility;
  }
  if(N_e <= N_a) {
    new_number[0] = (double)y;
    new_number[1] = (double)larger;
  }
  else {
    new_number[0] = (double)larger;
    new_number[1] = (double)y;
  }
  return ghl_success;
}

static ghl_error_codes_t ghl_m1_pair_make_effective_rates(
      const ghl_m1_neutrino_rates *restrict pair_rates,
      const double partner_n_com,
      const double partner_n_eq,
      ghl_m1_neutrino_rates *restrict effective_rates) {
  if(pair_rates == NULL || effective_rates == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(partner_n_com) || partner_n_com < 0.0 || !isfinite(partner_n_eq)
     || partner_n_eq <= 0.0) {
    return ghl_error_m1_invalid_state;
  }

  double eta_E = 0.0;
  double kappa_a_E = 0.0;
  if(pair_rates->J_eq <= 0.0 || !isfinite(pair_rates->J_eq)) {
    return ghl_error_m1_invalid_state;
  }
  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    eta_E += pair_rates->eta_E_pair[process];
    if(pair_rates->eta_E_pair[process] != 0.0) {
      kappa_a_E += pair_rates->eta_E_pair[process] / pair_rates->J_eq * partner_n_com
                   / partner_n_eq;
    }
  }
  if(!isfinite(eta_E) || !isfinite(kappa_a_E) || eta_E < 0.0 || kappa_a_E < 0.0) {
    return ghl_error_m1_invalid_state;
  }

  *effective_rates = *pair_rates;
  effective_rates->eta_N = 0.0;
  effective_rates->kappa_a_N = 0.0;
  effective_rates->eta_E = eta_E;
  effective_rates->kappa_a_E = kappa_a_E;
  effective_rates->kappa_s = 0.0;
  effective_rates->kappa_tr = kappa_a_E;
  effective_rates->eta_N_cc = 0.0;
  effective_rates->kappa_a_N_cc = 0.0;
  ghl_m1_pair_zero_fields(effective_rates);
  return ghl_success;
}

static void ghl_m1_pair_record_ef_repair(
      const ghl_m1_neutrino_state *restrict before,
      const ghl_m1_neutrino_state *restrict after,
      ghl_m1_neutrino_diagnostics *restrict diagnostics) {
  if(diagnostics == NULL) {
    return;
  }
  bool changed = before->E != after->E;
  diagnostics->repair_dE += fabs(after->E - before->E);
  for(int i = 0; i < 3; ++i) {
    changed |= before->F[i] != after->F[i];
    diagnostics->repair_dF[i] += fabs(after->F[i] - before->F[i]);
  }
  if(changed) {
    diagnostics->EF_repairs++;
  }
}

static ghl_error_codes_t ghl_m1_pair_attempt_schedule(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters nu_params[ghl_m1_pair_species_count],
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates rates[ghl_m1_pair_species_count],
      const ghl_m1_neutrino_state base_state[ghl_m1_pair_species_count],
      const ghl_m1_neutrino_diagnostics base_diagnostics[ghl_m1_pair_species_count],
      const double dt,
      const int num_substeps,
      ghl_m1_neutrino_state output_state[ghl_m1_pair_species_count],
      ghl_m1_neutrino_diagnostics output_diagnostics[ghl_m1_pair_species_count],
      ghl_m1_implicit_solve_diagnostics
            pair_solve_diagnostics[ghl_m1_pair_species_count]) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL || prims_frozen == NULL
     || rates == NULL || base_state == NULL || base_diagnostics == NULL
     || output_state == NULL || output_diagnostics == NULL
     || pair_solve_diagnostics == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(dt) || dt < 0.0 || num_substeps <= 0) {
    return ghl_error_m1_invalid_state;
  }

  ghl_m1_neutrino_state current_state[ghl_m1_pair_species_count]
        = { base_state[0], base_state[1] };
  ghl_m1_neutrino_diagnostics current_diagnostics[ghl_m1_pair_species_count]
        = { base_diagnostics[0], base_diagnostics[1] };
  ghl_m1_implicit_solve_diagnostics totals[ghl_m1_pair_species_count];
  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    ghl_m1_initialize_implicit_solve_diagnostics(&totals[species]);
  }

  const double dt_sub = dt / (double)num_substeps;
  const double h = metric->lapse * dt_sub;
  if(!isfinite(dt_sub) || !isfinite(h) || h < 0.0) {
    return ghl_error_m1_invalid_state;
  }

  for(int step = 0; step < num_substeps; ++step) {
    ghl_m1_neutrino_current current[ghl_m1_pair_species_count];
    for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
      const ghl_error_codes_t error = ghl_m1_neutrino_derive_current(
            m1_params, &nu_params[species], metric, prims_frozen,
            &current_state[species], &current[species]);
      if(error != ghl_success) {
        return error;
      }
    }

    double new_number[ghl_m1_pair_species_count];
    ghl_error_codes_t error
          = ghl_m1_pair_number_extent(rates, current, current_state, h, new_number);
    if(error != ghl_success) {
      return error;
    }

    ghl_m1_neutrino_state candidate_state[ghl_m1_pair_species_count]
          = { current_state[0], current_state[1] };
    candidate_state[0].N = new_number[0];
    candidate_state[1].N = new_number[1];
    for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
      if(!isfinite(candidate_state[species].N)
         || candidate_state[species].N < nu_params[species].N_floor) {
        return ghl_error_m1_implicit_admissibility;
      }
    }

    for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
      const int partner = 1 - species;
      const double partner_n_com = candidate_state[partner].N / current[partner].Gamma_N;
      ghl_m1_neutrino_rates effective_rates;
      error = ghl_m1_pair_make_effective_rates(
            &rates[species], partner_n_com, rates[partner].n_eq, &effective_rates);
      if(error != ghl_success) {
        return error;
      }

      const double sqrt_detgamma = metric->sqrt_detgamma;
      if(!isfinite(sqrt_detgamma) || sqrt_detgamma <= 0.0) {
        return ghl_error_m1_invalid_metric;
      }
      const double U_in[4] = { current_state[species].E * sqrt_detgamma,
                               current_state[species].F[0] * sqrt_detgamma,
                               current_state[species].F[1] * sqrt_detgamma,
                               current_state[species].F[2] * sqrt_detgamma };
      double U_out[4];
      ghl_m1_newton_diagnostics newton_diagnostics;
      bool closure_fallback_observed = false;
      error = ghl_m1_neutrino_attempt_EF_newton_step(
            m1_params, metric, prims_frozen, &effective_rates, dt_sub, U_in, U_out,
            &newton_diagnostics, &closure_fallback_observed);
      if(error != ghl_success) {
        return error;
      }

      totals[species].newton_iterations += newton_diagnostics.iterations;
      totals[species].line_search_backtracks += newton_diagnostics.backtracks;
      totals[species].residual_max_norm = newton_diagnostics.residual_max_norm;
      totals[species].residual_scaled_norm = newton_diagnostics.residual_weighted_merit;
      totals[species].used_fallback_substepping |= num_substeps > 1;
      totals[species].solution_path_flags
            |= num_substeps > 1 ? ghl_m1_solution_path_substepping
                                : ghl_m1_solution_path_primary_convergence;
      if(newton_diagnostics.used_projection) {
        totals[species].solution_path_flags |= ghl_m1_solution_path_projection;
      }
      if(newton_diagnostics.backtracks > 0) {
        totals[species].solution_path_flags
              |= ghl_m1_solution_path_line_search_backtracking;
      }
      if(closure_fallback_observed) {
        totals[species].solution_path_flags |= ghl_m1_solution_path_closure_fallback;
      }

      ghl_m1_neutrino_state repaired_state = candidate_state[species];
      repaired_state.E = U_out[0] / sqrt_detgamma;
      for(int i = 0; i < 3; ++i) {
        repaired_state.F[i] = U_out[i + 1] / sqrt_detgamma;
      }
      if(!ghl_m1_pair_state_is_finite(&repaired_state)) {
        return ghl_error_m1_invalid_state;
      }
      const ghl_m1_neutrino_state before_repair = repaired_state;
      ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(&repaired_state);
      error = ghl_m1_realizability_repair(m1_params, metric, &rad_state);
      if(error != ghl_success) {
        return error;
      }
      repaired_state.E = rad_state.E;
      for(int i = 0; i < 3; ++i) {
        repaired_state.F[i] = rad_state.F[i];
      }
      if(!ghl_m1_pair_state_is_finite(&repaired_state)) {
        return ghl_error_m1_invalid_state;
      }
      ghl_m1_pair_record_ef_repair(
            &before_repair, &repaired_state, &current_diagnostics[species]);
      candidate_state[species] = repaired_state;
    }

    current_state[0] = candidate_state[0];
    current_state[1] = candidate_state[1];
  }

  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    output_state[species] = current_state[species];
    output_diagnostics[species] = current_diagnostics[species];
    pair_solve_diagnostics[species] = totals[species];
    pair_solve_diagnostics[species].fallback_substeps = num_substeps;
  }
  return ghl_success;
}

static ghl_error_codes_t ghl_m1_pair_publish_failure(
      const ghl_error_codes_t error,
      const ghl_m1_neutrino_state state_transport[ghl_m1_pair_species_count],
      ghl_m1_neutrino_state state_out[ghl_m1_pair_species_count],
      ghl_m1_neutrino_exchange exchange[ghl_m1_pair_species_count],
      ghl_m1_neutrino_source_diagnostics diagnostics[ghl_m1_pair_species_count]) {
  if(state_transport != NULL && state_out != NULL) {
    for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
      state_out[species] = state_transport[species];
    }
  }
  if(exchange != NULL) {
    for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
      exchange[species] = (ghl_m1_neutrino_exchange){ 0 };
    }
  }
  if(diagnostics != NULL) {
    for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
      diagnostics[species].path = ghl_m1_neutrino_source_path_hard_failure;
      diagnostics[species].terminal_no_update = false;
    }
  }
  return error;
}

static ghl_error_codes_t ghl_m1_pair_publish_terminal(
      const ghl_error_codes_t error,
      const ghl_m1_neutrino_state state_transport[ghl_m1_pair_species_count],
      ghl_m1_neutrino_state state_out[ghl_m1_pair_species_count],
      ghl_m1_neutrino_exchange exchange[ghl_m1_pair_species_count],
      ghl_m1_neutrino_source_diagnostics diagnostics[ghl_m1_pair_species_count]) {
  const ghl_error_codes_t published = ghl_m1_pair_publish_failure(
        error, state_transport, state_out, exchange, diagnostics);
  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    diagnostics[species].path = ghl_m1_neutrino_source_path_terminal_no_update;
    diagnostics[species].terminal_no_update = true;
  }
  return published;
}

ghl_error_codes_t ghl_m1_solve_neutrino_pair_source_update(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters nu_params[ghl_m1_pair_species_count],
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates rates[ghl_m1_pair_species_count],
      const ghl_m1_neutrino_state state_input[ghl_m1_pair_species_count],
      const ghl_m1_neutrino_state state_transport[ghl_m1_pair_species_count],
      const double dt,
      const double n_b_cons,
      ghl_m1_neutrino_state state_out[ghl_m1_pair_species_count],
      ghl_m1_neutrino_exchange exchange[ghl_m1_pair_species_count],
      ghl_m1_neutrino_source_diagnostics diagnostics[ghl_m1_pair_species_count],
      ghl_m1_neutrino_diagnostics neutrino_diagnostics[ghl_m1_pair_species_count]) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL || prims_frozen == NULL
     || rates == NULL || state_input == NULL || state_transport == NULL
     || state_out == NULL || exchange == NULL || diagnostics == NULL
     || neutrino_diagnostics == NULL) {
    return ghl_error_m1_null_pointer;
  }

  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    state_out[species] = state_transport[species];
    exchange[species] = (ghl_m1_neutrino_exchange){ 0 };
    ghl_m1_pair_zero_source_diagnostics(&diagnostics[species]);
  }

  if(!isfinite(dt) || dt < 0.0 || !isfinite(n_b_cons) || n_b_cons <= 0.0) {
    return ghl_m1_pair_publish_failure(
          ghl_error_m1_invalid_state, state_transport, state_out, exchange, diagnostics);
  }
  if(rates[0].species != ghl_m1_neutrino_nue
     || rates[1].species != ghl_m1_neutrino_anue) {
    return ghl_m1_pair_publish_failure(
          ghl_error_m1_microphysics_failure, state_transport, state_out, exchange,
          diagnostics);
  }

  ghl_error_codes_t error = ghl_m1_validate_configuration(m1_params, metric);
  if(error != ghl_success) {
    return ghl_m1_pair_publish_failure(
          error, state_transport, state_out, exchange, diagnostics);
  }
  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    error = ghl_m1_validate_neutrino_rates(&rates[species], NULL);
    if(error != ghl_success) {
      return ghl_m1_pair_publish_failure(
            error, state_transport, state_out, exchange, diagnostics);
    }
  }
  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    if(rates[0].eta_N_pair[process] != rates[1].eta_N_pair[process]) {
      return ghl_m1_pair_publish_failure(
            ghl_error_m1_microphysics_failure, state_transport, state_out, exchange,
            diagnostics);
    }
  }

  /* Remove pair fields before crossing the established single-species source
   * boundary. This is the independent CC/scattering stage required by the
   * canonical paired contract. */
  ghl_m1_neutrino_rates base_rates[ghl_m1_pair_species_count] = { rates[0], rates[1] };
  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    ghl_m1_pair_zero_fields(&base_rates[species]);
  }
  ghl_m1_neutrino_state nonpair_state[ghl_m1_pair_species_count];
  ghl_m1_neutrino_exchange nonpair_exchange[ghl_m1_pair_species_count];
  ghl_m1_neutrino_diagnostics stage_diagnostics[ghl_m1_pair_species_count]
        = { neutrino_diagnostics[0], neutrino_diagnostics[1] };
  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    error = ghl_m1_solve_neutrino_source_update(
          NULL, m1_params, &nu_params[species], metric, prims_frozen,
          &base_rates[species], &state_input[species], &state_transport[species], dt,
          n_b_cons, &nonpair_state[species], &nonpair_exchange[species],
          &diagnostics[species], &stage_diagnostics[species]);
    if(error != ghl_success) {
      return ghl_m1_pair_publish_failure(
            error, state_transport, state_out, exchange, diagnostics);
    }
    if(stage_diagnostics[species].N_floor_repairs
       != neutrino_diagnostics[species].N_floor_repairs) {
      return ghl_m1_pair_publish_failure(
            ghl_error_m1_invalid_state, state_transport, state_out, exchange,
            diagnostics);
    }
  }

  if(!ghl_m1_pair_rates_are_active(rates)) {
    ghl_m1_neutrino_state final_state[ghl_m1_pair_species_count]
          = { nonpair_state[0], nonpair_state[1] };
    ghl_m1_neutrino_exchange final_exchange[ghl_m1_pair_species_count];
    for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
      error = ghl_m1_neutrino_assemble_exchange(
            &state_transport[species], &final_state[species], &base_rates[species],
            nonpair_exchange[species].dL_rad_cc, metric->sqrt_detgamma, n_b_cons,
            &final_exchange[species]);
      if(error != ghl_success) {
        return ghl_m1_pair_publish_failure(
              error, state_transport, state_out, exchange, diagnostics);
      }
    }
    for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
      state_out[species] = final_state[species];
      exchange[species] = final_exchange[species];
      diagnostics[species].path = ghl_m1_neutrino_source_path_general_implicit;
      neutrino_diagnostics[species] = stage_diagnostics[species];
    }
    return ghl_success;
  }

  const int fallback_schedule[] = { 1, 2, 4, 8, 16 };
  const int num_schedules
        = (int)(sizeof(fallback_schedule) / sizeof(fallback_schedule[0]));
  ghl_m1_neutrino_state pair_state[ghl_m1_pair_species_count];
  ghl_m1_neutrino_diagnostics pair_diagnostics[ghl_m1_pair_species_count];
  ghl_m1_implicit_solve_diagnostics pair_solve_diagnostics[ghl_m1_pair_species_count];
  int successful_substeps = 0;
  for(int schedule = 0; schedule < num_schedules; ++schedule) {
    error = ghl_m1_pair_attempt_schedule(
          m1_params, nu_params, metric, prims_frozen, rates, nonpair_state,
          stage_diagnostics, dt, fallback_schedule[schedule], pair_state,
          pair_diagnostics, pair_solve_diagnostics);
    if(error == ghl_success) {
      successful_substeps = fallback_schedule[schedule];
      break;
    }
    if(!ghl_m1_schedule_error_allows_retry(error)) {
      return ghl_m1_pair_publish_failure(
            error, state_transport, state_out, exchange, diagnostics);
    }
  }
  if(successful_substeps == 0) {
    return ghl_m1_pair_publish_terminal(
          ghl_error_m1_implicit_terminal_fallback, state_transport, state_out, exchange,
          diagnostics);
  }

  ghl_m1_neutrino_state final_state[ghl_m1_pair_species_count]
        = { pair_state[0], pair_state[1] };
  ghl_m1_neutrino_exchange final_exchange[ghl_m1_pair_species_count];
  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    ghl_m1_neutrino_current endpoint_current;
    error = ghl_m1_neutrino_derive_current(
          m1_params, &nu_params[species], metric, prims_frozen, &final_state[species],
          &endpoint_current);
    if(error != ghl_success) {
      return ghl_m1_pair_publish_failure(
            error, state_transport, state_out, exchange, diagnostics);
    }
    error = ghl_m1_neutrino_check_EN_bounds(
          &final_state[species], &nu_params[species], &endpoint_current);
    if(error != ghl_success) {
      return ghl_m1_pair_publish_failure(
            error, state_transport, state_out, exchange, diagnostics);
    }
    ghl_m1_neutrino_populate_mean_energy_diagnostics(
          &final_state[species], &endpoint_current, &nu_params[species], &rates[species],
          &pair_diagnostics[species]);

    error = ghl_m1_neutrino_assemble_exchange(
          &state_transport[species], &final_state[species], &base_rates[species],
          nonpair_exchange[species].dL_rad_cc, metric->sqrt_detgamma, n_b_cons,
          &final_exchange[species]);
    if(error != ghl_success) {
      return ghl_m1_pair_publish_failure(
            error, state_transport, state_out, exchange, diagnostics);
    }
  }
  for(int species = 0; species < ghl_m1_pair_species_count; ++species) {
    state_out[species] = final_state[species];
    exchange[species] = final_exchange[species];
    diagnostics[species].path = ghl_m1_neutrino_source_path_general_implicit;
    diagnostics[species].closure_fallback_used
          = diagnostics[species].closure_fallback_used
            || ((pair_solve_diagnostics[species].solution_path_flags
                 & ghl_m1_solution_path_closure_fallback)
                != 0u);
    diagnostics[species].implicit = pair_solve_diagnostics[species];
    diagnostics[species].implicit.solution_path_flags
          |= ghl_m1_solution_path_endpoint_acceptance;
    neutrino_diagnostics[species] = pair_diagnostics[species];
  }
  return ghl_success;
}
