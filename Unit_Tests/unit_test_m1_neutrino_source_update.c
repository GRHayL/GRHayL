#if defined(__linux__)
#define _GNU_SOURCE
#endif
#include "../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit.h"
#include "../GRHayL/Radiation/ghl_m1_utils.h"
#include "ghl_m1.h"
#include "m1_test_prng.h"

#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/*
 * Deterministic property-style coverage for the public neutrino source
 * boundaries.  The source update is pointwise and transactional, so this
 * test deliberately builds its own admissible states and frozen rate bundles
 * rather than depending on a table or a stale shared fixture.
 */

enum { SOURCE_RANDOM_CASES = 96, PAIR_RANDOM_CASES = 24 };

static void require_condition(
      const bool condition,
      const char *restrict message,
      const int case_index) {
  if(!condition) {
    ghl_error("M1 source-update case %d: %s\n", case_index, message);
  }
}

static void require_error(
      const ghl_error_codes_t actual,
      const ghl_error_codes_t expected,
      const char *restrict operation,
      const int case_index) {
  if(actual != expected) {
    ghl_error(
          "M1 source-update case %d: %s returned %d, expected %d\n", case_index,
          operation, (int)actual, (int)expected);
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
          "M1 source-update case %d: %s mismatch (got %.17e, expected %.17e)\n",
          case_index, quantity, actual, expected);
  }
}

static void require_zero_exchange(
      const ghl_m1_neutrino_exchange *restrict exchange,
      const char *restrict operation,
      const int case_index) {
  require_close(exchange->dN_rad_total, 0.0, 0.0, 0.0, operation, case_index);
  require_close(exchange->dL_rad_cc, 0.0, 0.0, 0.0, operation, case_index);
  require_close(exchange->dE_rad, 0.0, 0.0, 0.0, operation, case_index);
  require_close(exchange->dTau_matter, 0.0, 0.0, 0.0, operation, case_index);
  require_close(exchange->dYe_matter, 0.0, 0.0, 0.0, operation, case_index);
  for(int direction = 0; direction < 3; ++direction) {
    require_close(exchange->dF_rad[direction], 0.0, 0.0, 0.0, operation, case_index);
    require_close(exchange->dS_matter[direction], 0.0, 0.0, 0.0, operation, case_index);
  }
}

static void make_metric(
      m1_test_rng *restrict rng,
      const double lapse,
      ghl_metric_quantities *restrict metric) {
  /* Powers of two keep the diagonal inverse and unit determinant exact while
   * still varying the spatial geometry and all direction components. */
  const double diagonal_choices[3] = { 0.5, 1.0, 2.0 };
  const double g0 = diagonal_choices[m1_test_rng_next_u64(rng) % 3];
  const double g1 = diagonal_choices[m1_test_rng_next_u64(rng) % 3];
  const double g2 = 1.0 / (g0 * g1);

  *metric = (ghl_metric_quantities){ 0 };
  metric->lapse = lapse;
  metric->lapseinv = 1.0 / lapse;
  metric->lapseinv2 = metric->lapseinv * metric->lapseinv;
  metric->detgamma = 1.0;
  metric->sqrt_detgamma = 1.0;
  metric->gammaDD[0][0] = g0;
  metric->gammaDD[1][1] = g1;
  metric->gammaDD[2][2] = g2;
  metric->gammaUU[0][0] = 1.0 / g0;
  metric->gammaUU[1][1] = 1.0 / g1;
  metric->gammaUU[2][2] = 1.0 / g2;
}

static void make_primitives(ghl_primitive_quantities *restrict prims) {
  *prims = (ghl_primitive_quantities){ 0 };
  prims->rho = 1.0;
  prims->press = 0.1;
  prims->eps = 0.1;
  prims->Y_e = 0.5;
  prims->temperature = 1.0;
  prims->entropy = 1.0;
  /* vU is the coordinate three-velocity plus shift in the M1 convention. */
  prims->vU[0] = 0.0;
  prims->vU[1] = 0.0;
  prims->vU[2] = 0.0;
}

static void make_state(
      m1_test_rng *restrict rng,
      const ghl_metric_quantities *restrict metric,
      ghl_m1_neutrino_state *restrict state) {
  double direction[3]
        = { m1_test_rng_between(rng, -1.0, 1.0), m1_test_rng_between(rng, -1.0, 1.0),
            m1_test_rng_between(rng, -1.0, 1.0) };
  double direction_norm = sqrt(
        direction[0] * direction[0] + direction[1] * direction[1]
        + direction[2] * direction[2]);
  if(direction_norm == 0.0) {
    direction[0] = direction_norm = 1.0;
  }
  for(int i = 0; i < 3; ++i) {
    direction[i] /= direction_norm;
  }

  state->N = m1_test_rng_between(rng, 0.25, 1.75);
  state->E = m1_test_rng_between(rng, 0.5, 2.5);
  const double flux_factor = m1_test_rng_between(rng, 0.02, 0.68);
  const double covector_scale = state->E * flux_factor;
  for(int i = 0; i < 3; ++i) {
    state->F[i] = covector_scale * sqrt(metric->gammaDD[i][i]) * direction[i];
  }
}

static void make_neutrino_parameters(ghl_m1_neutrino_parameters *restrict nu_params) {
  *nu_params = (ghl_m1_neutrino_parameters){ .N_floor = 1.0e-12,
                                             .mean_energy_min = 0.0,
                                             .mean_energy_max = 0.0,
                                             .enforce_mean_energy_bounds = 0,
                                             .J_floor = 1.0e-14,
                                             .Gamma_N_floor = 1.0e-12 };
}

static void make_rates(
      const ghl_m1_neutrino_species_t species,
      const double kappa_a_N,
      const double kappa_a_E,
      const double kappa_s,
      const double n_eq,
      const double mean_energy,
      const double pair_number,
      const double pair_energy,
      ghl_m1_neutrino_rates *restrict rates) {
  *rates = (ghl_m1_neutrino_rates){ 0 };
  rates->species = species;
  rates->kappa_a_N = kappa_a_N;
  rates->kappa_a_E = kappa_a_E;
  rates->kappa_s = kappa_s;
  rates->kappa_tr = kappa_a_E + kappa_s;
  rates->n_eq = n_eq;
  rates->mean_energy = mean_energy;
  rates->J_eq = n_eq * mean_energy;
  rates->eta_N = kappa_a_N * n_eq;
  rates->eta_E = kappa_a_E * rates->J_eq;
  if(species == ghl_m1_neutrino_nue) {
    rates->lepton_weight = 1.0;
    rates->kappa_a_N_cc = kappa_a_N;
    rates->eta_N_cc = rates->eta_N;
  }
  else if(species == ghl_m1_neutrino_anue) {
    rates->lepton_weight = -1.0;
    rates->kappa_a_N_cc = kappa_a_N;
    rates->eta_N_cc = rates->eta_N;
  }
  else {
    rates->lepton_weight = 0.0;
  }
  for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
    rates->eta_N_pair[process] = pair_number;
    rates->eta_E_pair[process] = pair_energy;
  }
}

static void require_state_finite_and_admissible(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      const int case_index) {
  require_condition(
        isfinite(state->N) && isfinite(state->E), "published state has nonfinite N or E",
        case_index);
  require_condition(
        state->N >= nu_params->N_floor && state->E >= m1_params->E_floor,
        "published state violates configured floors", case_index);
  double flux_squared = 0.0;
  for(int i = 0; i < 3; ++i) {
    require_condition(
          isfinite(state->F[i]), "published state has a nonfinite flux", case_index);
    for(int j = 0; j < 3; ++j) {
      flux_squared += metric->gammaUU[i][j] * state->F[i] * state->F[j];
    }
  }
  require_condition(
        isfinite(flux_squared) && flux_squared >= 0.0, "published flux norm is invalid",
        case_index);
  require_condition(
        flux_squared <= state->E * state->E * (1.0 - m1_params->epsilon_c) + 1.0e-12,
        "published state violates the realizability cone", case_index);
}

static void require_exchange_contract(
      const ghl_m1_neutrino_state *restrict state_transport,
      const ghl_m1_neutrino_state *restrict state_out,
      const ghl_m1_neutrino_exchange *restrict exchange,
      const ghl_metric_quantities *restrict metric,
      const double n_b_cons,
      const int case_index) {
  require_close(
        exchange->dN_rad_total, state_out->N - state_transport->N, 1.0e-11, 1.0e-13,
        "total number exchange", case_index);
  require_close(
        exchange->dE_rad, state_out->E - state_transport->E, 1.0e-11, 1.0e-13,
        "energy exchange", case_index);
  require_close(
        exchange->dTau_matter, -metric->sqrt_detgamma * exchange->dE_rad, 1.0e-11,
        1.0e-13, "matter energy exchange", case_index);
  for(int direction = 0; direction < 3; ++direction) {
    require_close(
          exchange->dF_rad[direction],
          state_out->F[direction] - state_transport->F[direction], 1.0e-11, 1.0e-13,
          "flux exchange", case_index);
    require_close(
          exchange->dS_matter[direction],
          -metric->sqrt_detgamma * exchange->dF_rad[direction], 1.0e-11, 1.0e-13,
          "matter momentum exchange", case_index);
  }
  require_close(
        exchange->dYe_matter, -exchange->dL_rad_cc / n_b_cons, 1.0e-11, 1.0e-13,
        "electron-fraction exchange", case_index);
}

static ghl_m1_neutrino_source_options branched_options(void) {
  return (ghl_m1_neutrino_source_options){
    .policy = ghl_m1_neutrino_source_branched_compatibility,
    .thick_equilibrium_threshold = 0.5,
    .scattering_threshold = 0.5,
    .thermalized_number_threshold = -1.0,
    .allow_closure_fallback = true,
    .interaction_sources_already_applied = false,
    .ye_policy = ghl_m1_neutrino_ye_from_charged_current
  };
}

static void test_source_regimes(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  static const ghl_m1_neutrino_source_path_t expected_paths[4]
        = { ghl_m1_neutrino_source_path_thin_explicit,
            ghl_m1_neutrino_source_path_thick_equilibrium,
            ghl_m1_neutrino_source_path_scattering_dominated,
            ghl_m1_neutrino_source_path_general_implicit };
  static const double kappa_a_E[4] = { 0.10, 4.0, 0.0, 2.0 };
  static const double kappa_s[4] = { 0.10, 1.0, 4.0, 0.10 };

  for(int case_index = 0; case_index < SOURCE_RANDOM_CASES; ++case_index) {
    const int regime = case_index % 4;
    ghl_metric_quantities metric;
    make_metric(rng, m1_test_rng_between(rng, 0.55, 1.45), &metric);
    ghl_primitive_quantities prims;
    make_primitives(&prims);
    ghl_m1_neutrino_parameters nu_params;
    make_neutrino_parameters(&nu_params);
    ghl_m1_neutrino_rates rates;
    make_rates(
          ghl_m1_neutrino_nue, 0.08, kappa_a_E[regime], kappa_s[regime], 1.0, 2.0, 0.0,
          0.0, &rates);
    ghl_m1_neutrino_state state_input, state_transport, state_out;
    make_state(rng, &metric, &state_transport);
    state_input = state_transport;
    /* Deliberately make the pre-transport state distinct. The source base is
     * state_transport by contract, not state_input. */
    state_input.E *= 0.73;
    state_input.N *= 1.17;
    for(int i = 0; i < 3; ++i) {
      state_input.F[i] *= 0.73;
    }

    ghl_m1_neutrino_exchange exchange;
    ghl_m1_neutrino_source_diagnostics diagnostics;
    ghl_m1_neutrino_diagnostics neutrino_diagnostics;
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
    ghl_m1_neutrino_source_options options = branched_options();
    if(regime == 3) {
      options.thick_equilibrium_threshold = 0.0;
      options.scattering_threshold = 0.0;
    }
    ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
          &options, m1_params, &nu_params, &metric, &prims, &rates, &state_input,
          &state_transport, 1.0, 3.0, &state_out, &exchange, &diagnostics,
          &neutrino_diagnostics);
    require_error(error, ghl_success, "source regime update", case_index);
    require_condition(
          diagnostics.path == expected_paths[regime], "unexpected source-regime path",
          case_index);
    require_condition(
          !diagnostics.terminal_no_update, "source regime unexpectedly terminal",
          case_index);
    require_state_finite_and_admissible(
          m1_params, &nu_params, &metric, &state_out, case_index);
    require_exchange_contract(
          &state_transport, &state_out, &exchange, &metric, 3.0, case_index);
    require_condition(
          neutrino_diagnostics.source_converged > 0,
          "successful source path did not record convergence", case_index);
  }
}

static void test_source_base_and_lapse_scaling(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric_half, metric_one;
  make_metric(rng, 0.5, &metric_half);
  metric_one = metric_half;
  metric_one.lapse = 1.0;
  metric_one.lapseinv = 1.0;
  metric_one.lapseinv2 = 1.0;
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 0.12, 0.20, 0.10, 1.0, 2.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_state state_transport;
  make_state(rng, &metric_half, &state_transport);
  ghl_m1_neutrino_state state_input = state_transport;
  state_input.E *= 0.4;
  state_input.N *= 1.6;
  for(int i = 0; i < 3; ++i) {
    state_input.F[i] *= 0.4;
  }
  ghl_m1_neutrino_source_options options = branched_options();
  options.thick_equilibrium_threshold = 0.0;
  options.scattering_threshold = 0.0;
  ghl_m1_neutrino_state out_half, out_one, out_zero;
  ghl_m1_neutrino_exchange exchange_half, exchange_one, exchange_zero;
  ghl_m1_neutrino_source_diagnostics diagnostics_half, diagnostics_one, diagnostics_zero;
  ghl_m1_neutrino_diagnostics nd_half, nd_one, nd_zero;
  ghl_m1_neutrino_diagnostics_initialize(&nd_half);
  ghl_m1_neutrino_diagnostics_initialize(&nd_one);
  ghl_m1_neutrino_diagnostics_initialize(&nd_zero);

  ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
        &options, m1_params, &nu_params, &metric_half, &prims, &rates, &state_input,
        &state_transport, 0.4, 3.0, &out_half, &exchange_half, &diagnostics_half,
        &nd_half);
  require_error(error, ghl_success, "half-lapse source update", 1000);
  error = ghl_m1_solve_neutrino_source_update(
        &options, m1_params, &nu_params, &metric_one, &prims, &rates, &state_input,
        &state_transport, 0.2, 3.0, &out_one, &exchange_one, &diagnostics_one, &nd_one);
  require_error(error, ghl_success, "unit-lapse source update", 1001);
  for(int component = 0; component < 3; ++component) {
    require_close(
          out_half.F[component], out_one.F[component], 2.0e-10, 2.0e-12,
          "lapse-scaled flux", 1002);
  }
  require_close(out_half.N, out_one.N, 2.0e-10, 2.0e-12, "lapse-scaled number", 1002);
  require_close(out_half.E, out_one.E, 2.0e-10, 2.0e-12, "lapse-scaled energy", 1002);

  error = ghl_m1_solve_neutrino_source_update(
        &options, m1_params, &nu_params, &metric_half, &prims, &rates, &state_input,
        &state_transport, 0.0, 3.0, &out_zero, &exchange_zero, &diagnostics_zero,
        &nd_zero);
  require_error(error, ghl_success, "zero-timestep source update", 1003);
  require_condition(
        memcmp(&out_zero, &state_transport, sizeof(out_zero)) == 0,
        "zero-timestep update changed the source base", 1003);
  require_zero_exchange(&exchange_zero, "zero-timestep exchange", 1003);
}

static void test_scaled_ratio_underflow_consumers(void) {
  const double true_min = ldexp(DBL_MIN, -52);
  ghl_m1_neutrino_parameters nu_params = { 0 };
  ghl_m1_neutrino_rates rates = { 0 };
  rates.kappa_a_N = 1.0;
  rates.n_eq = 1.0;
  rates.J_eq = true_min;
  rates.mean_energy = true_min;
  const ghl_m1_neutrino_state state_base = { .N = 1.0 };
  double N_out = -1.0;

  /* If the final quotient itself is not representable, the endpoint must
   * remain unpublished rather than silently becoming zero. */
  rates.mean_energy = DBL_MAX;
  const ghl_m1_neutrino_current nonrepresentable_current
        = { .J = DBL_MIN, .Gamma_N = DBL_MIN };
  N_out = -2.0;
  require_error(
        ghl_m1_neutrino_update_endpoint_number_with_policy(
              &nu_params, &rates, 0.0, 0.0, 0.0, &state_base, &nonrepresentable_current,
              &N_out, NULL),
        ghl_error_m1_invalid_state, "nonrepresentable projected endpoint", 1211);
  require_condition(N_out == -2.0, "failed projected endpoint changed its output", 1211);
}

static void
test_thermalized_number_projection(const ghl_m1_parameters *restrict m1_params) {
  m1_test_rng rng = { .state = UINT64_C(0x4d315f544845524d) };
  ghl_metric_quantities metric;
  make_metric(&rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, &rates);
  const ghl_m1_neutrino_state state_transport
        = { .N = 0.25, .E = 2.0, .F = { 0.0, 0.0, 0.0 } };
  const ghl_m1_neutrino_state state_input = state_transport;
  ghl_m1_neutrino_source_options ordinary_options = branched_options();
  ghl_m1_neutrino_source_options projection_options = ordinary_options;
  projection_options.thermalized_number_threshold = 0.0;
  ghl_m1_neutrino_state ordinary_out, projection_out;
  ghl_m1_neutrino_exchange ordinary_exchange, projection_exchange;
  ghl_m1_neutrino_source_diagnostics ordinary_diagnostics, projection_diagnostics;
  ghl_m1_neutrino_diagnostics ordinary_nd, projection_nd;
  ghl_m1_neutrino_diagnostics_initialize(&ordinary_nd);
  ghl_m1_neutrino_diagnostics_initialize(&projection_nd);

  ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
        &ordinary_options, m1_params, &nu_params, &metric, &prims, &rates, &state_input,
        &state_transport, 0.0, 3.0, &ordinary_out, &ordinary_exchange,
        &ordinary_diagnostics, &ordinary_nd);
  require_error(error, ghl_success, "ordinary number endpoint", 1100);
  error = ghl_m1_solve_neutrino_source_update(
        &projection_options, m1_params, &nu_params, &metric, &prims, &rates,
        &state_input, &state_transport, 0.0, 3.0, &projection_out, &projection_exchange,
        &projection_diagnostics, &projection_nd);
  require_error(error, ghl_success, "thermalized number endpoint", 1101);
  require_condition(
        ordinary_diagnostics.path == ghl_m1_neutrino_source_path_thin_explicit
              && projection_diagnostics.path
                       == ghl_m1_neutrino_source_path_thin_explicit,
        "thermalized number case did not use the thin endpoint", 1101);
  /* A branch that runs no Newton iteration must still publish a coherent
   * implicit record: ghl_m1_implicit_solve_diagnostics documents that success
   * implies residual_scaled_norm <= 1, so the initializer's INFINITY sentinels
   * must not survive beside the convergence flags. */
  require_condition(
        isfinite(ordinary_diagnostics.implicit.residual_scaled_norm)
              && ordinary_diagnostics.implicit.residual_scaled_norm <= 1.0
              && isfinite(ordinary_diagnostics.implicit.residual_max_norm)
              && ordinary_diagnostics.implicit.newton_iterations == 0
              && ordinary_diagnostics.implicit.line_search_backtracks == 0
              && ordinary_diagnostics.implicit.fallback_substeps == 1
              && !ordinary_diagnostics.implicit.used_fallback_substepping,
        "thin branch published an incoherent implicit solve record", 1102);
  require_condition(
        isfinite(projection_diagnostics.implicit.residual_scaled_norm)
              && projection_diagnostics.implicit.residual_scaled_norm <= 1.0
              && projection_diagnostics.implicit.newton_iterations == 0,
        "projected thin branch published an incoherent implicit solve record", 1102);

  require_close(
        ordinary_out.N, state_transport.N, 0.0, 0.0, "ordinary number endpoint", 1100);
  require_close(
        projection_out.N, 1.0, 1.0e-12, 1.0e-12, "thermalized number projection", 1101);
  require_condition(
        projection_out.N != ordinary_out.N, "thermalized number threshold had no effect",
        1101);
  require_exchange_contract(
        &state_transport, &projection_out, &projection_exchange, &metric, 3.0, 1101);
}

static void test_branched_general_thermalized_number_projection(
      const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 0.0, &rates);

  const ghl_m1_neutrino_state state_transport
        = { .N = 0.5, .E = 4.0, .F = { 0.0, 0.0, 0.0 } };
  const ghl_m1_neutrino_state state_input = state_transport;

  /* Make both stiff shortcuts ineligible without testing their private
   * selector order. The flat/rest-frame E/F endpoint is independently
   * (4 + 1*4)/(1 + 1*2) = 8/3. */
  ghl_m1_neutrino_source_options projected_options = branched_options();
  projected_options.thick_equilibrium_threshold = DBL_MAX;
  projected_options.scattering_threshold = DBL_MAX;
  projected_options.thermalized_number_threshold = 0.0;
  ghl_m1_neutrino_source_options disabled_options = projected_options;
  disabled_options.thermalized_number_threshold = -1.0;
  const ghl_m1_neutrino_source_options *const options[3]
        = { &disabled_options, NULL, &projected_options };
  const double expected_N[3] = { 0.75, 0.75, 4.0 / 3.0 };
  const double expected_dL_rad_cc[3] = { 0.25, 0.25, -1.0 / 3.0 };
  const double expected_dYe_matter[3] = { -1.0 / 12.0, -1.0 / 12.0, 1.0 / 9.0 };
  const char *const labels[3] = { "branched general disabled thermalized endpoint",
                                  "default implicit disabled thermalized endpoint",
                                  "branched general thermalized endpoint" };

  for(int case_index = 0; case_index < 3; ++case_index) {
    ghl_m1_neutrino_state state_out;
    ghl_m1_neutrino_exchange exchange;
    ghl_m1_neutrino_source_diagnostics diagnostics;
    ghl_m1_neutrino_diagnostics neutrino_diagnostics;
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
    const ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
          options[case_index], m1_params, &nu_params, &metric, &prims, &rates,
          &state_input, &state_transport, 1.0, 3.0, &state_out, &exchange, &diagnostics,
          &neutrino_diagnostics);
    require_error(error, ghl_success, labels[case_index], 1102 + case_index);
    require_condition(
          diagnostics.path == ghl_m1_neutrino_source_path_general_implicit
                && !diagnostics.terminal_no_update,
          "thermalized general case did not use the general endpoint",
          1102 + case_index);
    require_close(
          state_out.E, 8.0 / 3.0, 2.0e-11, 2.0e-13,
          "thermalized general energy endpoint", 1102 + case_index);
    require_close(
          state_out.N, expected_N[case_index], 2.0e-11, 2.0e-13,
          "thermalized general number endpoint", 1102 + case_index);
    for(int direction = 0; direction < 3; ++direction) {
      require_close(
            state_out.F[direction], 0.0, 0.0, 0.0, "thermalized general flux endpoint",
            1102 + case_index);
    }
    require_exchange_contract(
          &state_transport, &state_out, &exchange, &metric, 3.0, 1102 + case_index);
    require_close(
          exchange.dL_rad_cc, expected_dL_rad_cc[case_index], 2.0e-11, 2.0e-13,
          "thermalized endpoint charged-current exchange", 1102 + case_index);
    require_close(
          exchange.dYe_matter, expected_dYe_matter[case_index], 2.0e-11, 2.0e-13,
          "thermalized endpoint Y_e exchange", 1102 + case_index);
  }
}

static void
test_stiff_branch_arithmetic_boundaries(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  const ghl_m1_neutrino_state state_transport
        = { .N = 1.0, .E = 2.0, .F = { 0.0, 0.0, 0.0 } };
  const ghl_m1_neutrino_state state_input = state_transport;

  /* Powers of two make the intended products exact before the deliberately
   * overflowing/underflowing intermediate. All three rate bundles remain
   * finite and satisfy the public rate identities. */
  const double opacity_large = ldexp(1.0, 600);
  const double dt_large = ldexp(1.0, 600);
  const double dt_small = ldexp(1.0, -600);
  require_condition(
        isfinite(opacity_large) && isfinite(dt_large) && isfinite(dt_small),
        "stiff arithmetic fixture is not finite", 1290);

  ghl_m1_neutrino_rates rates[2];
  make_rates(
        ghl_m1_neutrino_nue, 0.0, opacity_large, 0.0, 1.0, 1.0, 0.0, 0.0, &rates[0]);
  make_rates(
        ghl_m1_neutrino_nue, 0.0, 0.0, opacity_large, 1.0, 1.0, 0.0, 0.0, &rates[1]);
  const double dt[2] = { dt_small, dt_large };
  const double expected_E[2] = { 1.5, 2.0 };
  const ghl_m1_neutrino_source_path_t expected_path[2]
        = { ghl_m1_neutrino_source_path_thick_equilibrium,
            ghl_m1_neutrino_source_path_scattering_dominated };
  const char *const labels[2]
        = { "thick opacity-product overflow", "scattering optical-depth overflow" };

  for(int case_index = 0; case_index < 2; ++case_index) {
    ghl_m1_neutrino_source_options options = branched_options();
    ghl_m1_neutrino_state state_out;
    ghl_m1_neutrino_exchange exchange;
    ghl_m1_neutrino_source_diagnostics diagnostics;
    ghl_m1_neutrino_diagnostics neutrino_diagnostics;
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
    const ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
          &options, m1_params, &nu_params, &metric, &prims, &rates[case_index],
          &state_input, &state_transport, dt[case_index], 3.0, &state_out, &exchange,
          &diagnostics, &neutrino_diagnostics);
    require_error(error, ghl_success, labels[case_index], 1291 + case_index);
    require_condition(
          diagnostics.path == expected_path[case_index]
                && !diagnostics.terminal_no_update,
          "stiff arithmetic selected the wrong public path", 1291 + case_index);
    require_close(
          state_out.N, state_transport.N, 0.0, 0.0, "stiff arithmetic number endpoint",
          1291 + case_index);
    require_close(
          state_out.E, expected_E[case_index], 2.0e-11, 2.0e-13,
          "stiff arithmetic energy endpoint", 1291 + case_index);
    for(int direction = 0; direction < 3; ++direction) {
      require_close(
            state_out.F[direction], 0.0, 0.0, 0.0, "stiff arithmetic flux endpoint",
            1291 + case_index);
    }
    require_exchange_contract(
          &state_transport, &state_out, &exchange, &metric, 3.0, 1291 + case_index);
  }

  /* The endpoint below is genuinely above DBL_MAX, not merely an overflowing
   * intermediate. In flat space with v=3/5 and chi=1/3, the independent
   * isotropic boost oracle gives W^2(1+v^2/3)=7/4,
   * J_transport=(7/8) DBL_MAX, J_endpoint=(67/72) DBL_MAX, and
   * E_endpoint=(469/288) DBL_MAX. The margins are deliberately wide enough
   * that the failure does not depend on compiler expression contraction. */
  prims.vU[0] = 0.6;
  prims.u0 = 1.25;
  const long double velocity = 0.6L;
  const long double W_squared = 1.0L / (1.0L - velocity * velocity);
  const long double isotropic_boost_factor
        = W_squared * (1.0L + velocity * velocity / 3.0L);
  const long double J_transport_over_max = 0.5L * isotropic_boost_factor;
  const long double dtau = 4.0L / 5.0L;
  const long double J_endpoint_over_max = (J_transport_over_max + dtau) / (1.0L + dtau);
  const long double E_endpoint_over_max = isotropic_boost_factor * J_endpoint_over_max;
  require_condition(
        E_endpoint_over_max > 1.5L,
        "nonrepresentable endpoint oracle is not above DBL_MAX", 1295);

  ghl_m1_neutrino_rates nonrepresentable_rates;
  make_rates(
        ghl_m1_neutrino_nue, 0.0, 1.0, 0.0, DBL_MAX, 1.0, 0.0, 0.0,
        &nonrepresentable_rates);
  const ghl_m1_neutrino_state nonrepresentable_state
        = { .N = 1.0, .E = 0.5 * DBL_MAX, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_neutrino_source_options nonrepresentable_options = branched_options();
  ghl_m1_neutrino_state state_out
        = { .N = -41.0, .E = -42.0, .F = { -43.0, -44.0, -45.0 } };
  ghl_m1_neutrino_exchange exchange = { .dN_rad_total = -46.0,
                                        .dL_rad_cc = -47.0,
                                        .dE_rad = -48.0,
                                        .dF_rad = { -49.0, -50.0, -51.0 },
                                        .dTau_matter = -52.0,
                                        .dS_matter = { -53.0, -54.0, -55.0 },
                                        .dYe_matter = -56.0 };
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics neutrino_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &nonrepresentable_options, m1_params, &nu_params, &metric, &prims,
              &nonrepresentable_rates, &nonrepresentable_state, &nonrepresentable_state,
              1.0, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "genuinely nonrepresentable thick endpoint", 1295);
  require_condition(
        memcmp(&state_out, &nonrepresentable_state, sizeof(state_out)) == 0,
        "nonrepresentable thick endpoint changed state", 1295);
  require_zero_exchange(&exchange, "nonrepresentable thick endpoint exchange", 1295);
  require_condition(
        diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && !diagnostics.terminal_no_update
              && neutrino_diagnostics.source_failures == 1,
        "nonrepresentable endpoint failure was not transactional", 1295);

  /* A positive thick threshold with both opacities zero selects no thick
   * branch: the scaled stiffness product is exactly zero and must compare
   * below the threshold without saturating. */
  ghl_m1_neutrino_rates zero_opacity_rates;
  make_rates(
        ghl_m1_neutrino_nue, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, &zero_opacity_rates);
  ghl_m1_neutrino_source_options zero_opacity_options = branched_options();
  zero_opacity_options.thick_equilibrium_threshold = 0.5;
  ghl_m1_neutrino_state zero_opacity_out;
  ghl_m1_neutrino_exchange zero_opacity_exchange;
  ghl_m1_neutrino_source_diagnostics zero_opacity_diagnostics;
  ghl_m1_neutrino_diagnostics zero_opacity_nd;
  ghl_m1_neutrino_diagnostics_initialize(&zero_opacity_nd);
  const ghl_error_codes_t zero_opacity_error = ghl_m1_solve_neutrino_source_update(
        &zero_opacity_options, m1_params, &nu_params, &metric, &prims,
        &zero_opacity_rates, &state_input, &state_transport, 1.0, 3.0, &zero_opacity_out,
        &zero_opacity_exchange, &zero_opacity_diagnostics, &zero_opacity_nd);
  require_error(zero_opacity_error, ghl_success, "zero-opacity thick selection", 1296);
  require_condition(
        zero_opacity_diagnostics.path == ghl_m1_neutrino_source_path_thin_explicit,
        "zero-opacity stiffness did not skip the thick branch", 1296);
  require_state_finite_and_admissible(
        m1_params, &nu_params, &metric, &zero_opacity_out, 1296);
  require_exchange_contract(
        &state_transport, &zero_opacity_out, &zero_opacity_exchange, &metric, 3.0, 1296);

  /* A thick-selected step whose dtau*kappa_a_E overflows takes the ratio's
   * scaled-denominator fallback inside the predictor while the endpoint stays
   * at the representable equilibrium. Restore the rest fluid first: the
   * oracle above deliberately runs at v = 3/5. */
  make_primitives(&prims);
  ghl_m1_neutrino_rates overflow_predictor_rates;
  make_rates(
        ghl_m1_neutrino_nue, 0.0, opacity_large, 0.0, 1.0, 1.0, 0.0, 0.0,
        &overflow_predictor_rates);
  ghl_m1_neutrino_source_options overflow_predictor_options = branched_options();
  overflow_predictor_options.thick_equilibrium_threshold = 0.5;
  ghl_m1_neutrino_state overflow_predictor_out;
  ghl_m1_neutrino_exchange overflow_predictor_exchange;
  ghl_m1_neutrino_source_diagnostics overflow_predictor_diagnostics;
  ghl_m1_neutrino_diagnostics overflow_predictor_nd;
  ghl_m1_neutrino_diagnostics_initialize(&overflow_predictor_nd);
  const ghl_error_codes_t overflow_predictor_error = ghl_m1_solve_neutrino_source_update(
        &overflow_predictor_options, m1_params, &nu_params, &metric, &prims,
        &overflow_predictor_rates, &state_input, &state_transport, dt_large, 3.0,
        &overflow_predictor_out, &overflow_predictor_exchange,
        &overflow_predictor_diagnostics, &overflow_predictor_nd);
  require_error(
        overflow_predictor_error, ghl_success, "thick predictor denominator overflow",
        1297);
  require_condition(
        overflow_predictor_diagnostics.path
              == ghl_m1_neutrino_source_path_thick_equilibrium,
        "thick predictor overflow selected the wrong path", 1297);
  require_state_finite_and_admissible(
        m1_params, &nu_params, &metric, &overflow_predictor_out, 1297);
  require_exchange_contract(
        &state_transport, &overflow_predictor_out, &overflow_predictor_exchange, &metric,
        3.0, 1297);

  /* A predictor whose dtau*kappa overflows while dtau*eta stays finite takes
   * the ratio's denominator fallback and still lands on the equilibrium. */
  make_primitives(&prims);
  ghl_m1_neutrino_rates denominator_overflow_rates;
  make_rates(
        ghl_m1_neutrino_nue, 0.0, opacity_large, 0.0, 1.0, dt_small, 0.0, 0.0,
        &denominator_overflow_rates);
  require_error(
        ghl_m1_validate_neutrino_rates(&denominator_overflow_rates, NULL), ghl_success,
        "denominator-overflow predictor rates", 1298);
  ghl_m1_neutrino_source_options denominator_overflow_options = branched_options();
  denominator_overflow_options.thick_equilibrium_threshold = 0.5;
  ghl_m1_neutrino_state denominator_overflow_out;
  ghl_m1_neutrino_exchange denominator_overflow_exchange;
  ghl_m1_neutrino_source_diagnostics denominator_overflow_diagnostics;
  ghl_m1_neutrino_diagnostics denominator_overflow_nd;
  ghl_m1_neutrino_diagnostics_initialize(&denominator_overflow_nd);
  const ghl_error_codes_t denominator_overflow_error
        = ghl_m1_solve_neutrino_source_update(
              &denominator_overflow_options, m1_params, &nu_params, &metric, &prims,
              &denominator_overflow_rates, &state_input, &state_transport, dt_large, 3.0,
              &denominator_overflow_out, &denominator_overflow_exchange,
              &denominator_overflow_diagnostics, &denominator_overflow_nd);
  require_error(
        denominator_overflow_error, ghl_success, "thick predictor denominator overflow",
        1298);
  require_condition(
        denominator_overflow_diagnostics.path
              == ghl_m1_neutrino_source_path_thick_equilibrium,
        "denominator overflow selected the wrong path", 1298);
  require_state_finite_and_admissible(
        m1_params, &nu_params, &metric, &denominator_overflow_out, 1298);
  require_exchange_contract(
        &state_transport, &denominator_overflow_out, &denominator_overflow_exchange,
        &metric, 3.0, 1298);
}

static void
test_reachable_scaled_ratio_failure(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);

  /* The one-ULP excess in eta_E is within the provider validator's stated
   * representational slack, but makes the exact endpoint eta_E/kappa_a_E
   * exceed DBL_MAX.  The direct BE products overflow, so the public thick
   * branch must use the scaled fallback and reject that nonrepresentable
   * quotient transactionally. */
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nux, 0.0, 0.5, 0.0, DBL_MAX, 1.0, 0.0, 0.0, &rates);
  rates.eta_E = nextafter(0.5 * DBL_MAX, INFINITY);
  ghl_m1_neutrino_diagnostics rate_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&rate_diagnostics);
  require_error(
        ghl_m1_validate_neutrino_rates(&rates, &rate_diagnostics), ghl_success,
        "representational-slack rates were rejected", 1297);

  const ghl_m1_neutrino_state state_transport
        = { .N = 1.0, .E = DBL_MAX, .F = { 0.0, 0.0, 0.0 } };
  const ghl_m1_neutrino_state state_input = state_transport;
  ghl_m1_neutrino_source_options options = branched_options();
  ghl_m1_neutrino_state state_out
        = { .N = -41.0, .E = -42.0, .F = { -43.0, -44.0, -45.0 } };
  ghl_m1_neutrino_exchange exchange = { .dN_rad_total = -46.0,
                                        .dL_rad_cc = -47.0,
                                        .dE_rad = -48.0,
                                        .dF_rad = { -49.0, -50.0, -51.0 },
                                        .dTau_matter = -52.0,
                                        .dS_matter = { -53.0, -54.0, -55.0 },
                                        .dYe_matter = -56.0 };
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics neutrino_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);

  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state_input,
              &state_transport, DBL_MAX, 3.0, &state_out, &exchange, &diagnostics,
              &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "representational-slack thick endpoint was accepted",
        1297);
  require_condition(
        memcmp(&state_out, &state_transport, sizeof(state_out)) == 0,
        "scaled-ratio failure changed the state", 1297);
  require_zero_exchange(&exchange, "scaled-ratio failure exchange", 1297);
  require_condition(
        diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && !diagnostics.terminal_no_update
              && neutrino_diagnostics.source_failures == 1,
        "scaled-ratio failure was not published transactionally", 1297);
}

static void test_public_source_api_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  metric.betaU[1] = 0.07;
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 0.08, 0.10, 0.10, 1.0, 2.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_state state;
  make_state(rng, &metric, &state);
  const ghl_m1_rad_state rad_state
        = { .E = state.E, .F = { state.F[0], state.F[1], state.F[2] } };

  bool closure_fallback_observed = false;
  ghl_m1_sources diagnostic_sources = { 0 };
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              m1_params, &metric, &prims, &rad_state, &rates, &closure_fallback_observed,
              &diagnostic_sources),
        ghl_success, "diagnostic E/F interaction sources", 1200);

  /* A moving fluid with zero Eulerian flux selects the documented closure
   * endpoint fallback. The diagnostic entry point must expose that status to
   * its caller while still publishing finite interaction sources. */
  ghl_primitive_quantities moving_prims = prims;
  moving_prims.vU[0] = 0.4;
  moving_prims.u0 = 1.0 / sqrt(1.0 - moving_prims.vU[0] * moving_prims.vU[0]);
  const ghl_m1_rad_state zero_flux_rad_state = { .E = 1.0, .F = { 0.0, 0.0, 0.0 } };
  closure_fallback_observed = false;
  ghl_m1_sources fallback_sources = { 0 };
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              m1_params, &metric, &moving_prims, &zero_flux_rad_state, &rates,
              &closure_fallback_observed, &fallback_sources),
        ghl_success, "diagnostic closure fallback sources", 1205);
  require_condition(
        closure_fallback_observed, "diagnostic closure fallback was not reported", 1205);
  require_condition(
        isfinite(fallback_sources.S_E),
        "diagnostic closure fallback source was nonfinite", 1205);

  ghl_m1_sources public_sources = { 0 };
  double number_source = NAN;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, &prims, &state, &rates, &public_sources,
              &number_source),
        ghl_success, "public interaction sources", 1201);
  require_close(
        diagnostic_sources.S_E, public_sources.S_E, 2.0e-11, 2.0e-13,
        "diagnostic energy source", 1202);
  for(int direction = 0; direction < 3; ++direction) {
    require_close(
          diagnostic_sources.S[direction], public_sources.S[direction], 2.0e-11, 2.0e-13,
          "diagnostic momentum source", 1203 + direction);
  }

  /* Exercise every public interaction-source output/input guard through the
   * installed API.  Each call is made with otherwise valid operands so the
   * selected short-circuit disjunct is the one under test; the initialized
   * source outputs must remain untouched on failure. */
  ghl_m1_sources source_sentinel = { .S_E = -71.0, .S = { -72.0, -73.0, -74.0 } };
  double number_sentinel = -75.0;
  /* A densitization overflow inside the current walk rejects the whole
   * interaction-source call. */
  {
    ghl_metric_quantities huge_volume_metric;
    ghl_initialize_metric(
          1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &huge_volume_metric);
    huge_volume_metric.sqrt_detgamma = 1.0e154;
    huge_volume_metric.detgamma = 1.0e308;
    for(int i = 0; i < 3; ++i) {
      huge_volume_metric.gammaDD[i][i] = 1.0e308;
      huge_volume_metric.gammaUU[i][i] = 1.0e-308;
    }
    ghl_m1_sources overflow_sources = { .S_E = -81.0, .S = { -82.0, -83.0, -84.0 } };
    double overflow_number = -85.0;
    require_error(
          ghl_m1_compute_neutrino_interaction_sources(
                m1_params, &nu_params, &huge_volume_metric, &prims, &state, &rates,
                &overflow_sources, &overflow_number),
          ghl_error_m1_invalid_state, "interaction current densitization overflow",
          1204);
    ghl_m1_neutrino_state huge_number_state = state;
    huge_number_state.N = DBL_MAX;
    ghl_m1_sources huge_number_sources = { .S_E = -86.0, .S = { -87.0, -88.0, -89.0 } };
    double huge_number_source = -90.0;
    require_error(
          ghl_m1_compute_neutrino_interaction_sources(
                m1_params, &nu_params, &huge_volume_metric, &prims, &huge_number_state,
                &rates, &huge_number_sources, &huge_number_source),
          ghl_error_m1_invalid_state, "interaction number densitization overflow", 1204);
    require_condition(
          huge_number_sources.S_E == -86.0 && huge_number_source == -90.0,
          "interaction number overflow changed outputs", 1204);
    /* A number floor above the comoving number rejects the current walk. */
    ghl_m1_neutrino_parameters huge_floor_nu = nu_params;
    huge_floor_nu.J_floor = DBL_MAX;
    ghl_m1_sources huge_floor_sources = { .S_E = -91.0, .S = { -92.0, -93.0, -94.0 } };
    double huge_floor_number = -95.0;
    require_error(
          ghl_m1_compute_neutrino_interaction_sources(
                m1_params, &huge_floor_nu, &metric, &prims, &state, &rates,
                &huge_floor_sources, &huge_floor_number),
          ghl_error_m1_invalid_state, "interaction number-floor rejection", 1204);
    require_condition(
          huge_floor_sources.S_E == -91.0 && huge_floor_number == -95.0,
          "interaction floor rejection changed outputs", 1204);
    require_condition(
          overflow_sources.S_E == -81.0 && overflow_sources.S[0] == -82.0
                && overflow_number == -85.0,
          "interaction overflow changed outputs", 1204);
  }
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              NULL, &nu_params, &metric, &prims, &state, &rates, &source_sentinel,
              &number_sentinel),
        ghl_error_m1_null_pointer, "null interaction M1 parameters", 1204);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, NULL, &metric, &prims, &state, &rates, &source_sentinel,
              &number_sentinel),
        ghl_error_m1_null_pointer, "null interaction number parameters", 1205);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, NULL, &prims, &state, &rates, &source_sentinel,
              &number_sentinel),
        ghl_error_m1_null_pointer, "null interaction metric", 1206);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, NULL, &state, &rates, &source_sentinel,
              &number_sentinel),
        ghl_error_m1_null_pointer, "null interaction primitives", 1207);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, &prims, NULL, &rates, &source_sentinel,
              &number_sentinel),
        ghl_error_m1_null_pointer, "null interaction state", 1208);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, &prims, &state, NULL, &source_sentinel,
              &number_sentinel),
        ghl_error_m1_null_pointer, "null interaction rates", 1209);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, &prims, &state, &rates, NULL,
              &number_sentinel),
        ghl_error_m1_null_pointer, "null interaction EF output", 1210);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, &prims, &state, &rates, &source_sentinel,
              NULL),
        ghl_error_m1_null_pointer, "null interaction number output", 1211);
  if(source_sentinel.S_E != -71.0 || source_sentinel.S[0] != -72.0
     || number_sentinel != -75.0) {
    ghl_error("M1 source-update case 1212: null interaction guard published outputs\n");
  }

  ghl_m1_closure supplied_closure;
  require_error(
        ghl_m1_compute_neutrino_closure(
              m1_params, &metric, &prims, &state, &supplied_closure),
        ghl_success, "interaction closure construction", 1213);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              NULL, &nu_params, &metric, &prims, &state, &supplied_closure, &rates,
              &source_sentinel, &number_sentinel),
        ghl_error_m1_null_pointer, "null closure-source M1 parameters", 1214);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, NULL, &metric, &prims, &state, &supplied_closure, &rates,
              &source_sentinel, &number_sentinel),
        ghl_error_m1_null_pointer, "null closure-source number parameters", 1215);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, NULL, &prims, &state, &supplied_closure, &rates,
              &source_sentinel, &number_sentinel),
        ghl_error_m1_null_pointer, "null closure-source metric", 1216);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, NULL, &state, &supplied_closure, &rates,
              &source_sentinel, &number_sentinel),
        ghl_error_m1_null_pointer, "null closure-source primitives", 1217);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, &prims, NULL, &supplied_closure, &rates,
              &source_sentinel, &number_sentinel),
        ghl_error_m1_null_pointer, "null closure-source state", 1218);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, &prims, &state, NULL, &rates,
              &source_sentinel, &number_sentinel),
        ghl_error_m1_null_pointer, "null closure-source closure", 1219);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, &prims, &state, &supplied_closure, NULL,
              &source_sentinel, &number_sentinel),
        ghl_error_m1_null_pointer, "null closure-source rates", 1220);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, &prims, &state, &supplied_closure, &rates,
              NULL, &number_sentinel),
        ghl_error_m1_null_pointer, "null closure-source EF output", 1221);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, &prims, &state, &supplied_closure, &rates,
              &source_sentinel, NULL),
        ghl_error_m1_null_pointer, "null closure-source number output", 1222);
  ghl_m1_neutrino_parameters bad_source_nu = nu_params;
  bad_source_nu.N_floor = NAN;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &bad_source_nu, &metric, &prims, &state, &supplied_closure,
              &rates, &source_sentinel, &number_sentinel),
        ghl_error_m1_invalid_state, "invalid closure-source number floor", 1223);
  ghl_m1_neutrino_state bad_source_state = state;
  bad_source_state.N = NAN;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, &prims, &bad_source_state,
              &supplied_closure, &rates, &source_sentinel, &number_sentinel),
        ghl_error_m1_invalid_state, "invalid closure-source state", 1224);

  double number_flux[3], number_transport_velocity[3];
  require_error(
        ghl_m1_compute_neutrino_number_flux(
              m1_params, &nu_params, &metric, &prims, &state, number_flux,
              number_transport_velocity),
        ghl_success, "public number-current flux", 1210);
  for(int direction = 0; direction < 3; ++direction) {
    require_condition(
          isfinite(number_flux[direction])
                && isfinite(number_transport_velocity[direction]),
          "public number-current flux is nonfinite", 1210 + direction);
    double physical_number_flux = NAN;
    require_error(
          ghl_m1_compute_neutrino_physical_number_flux(
                m1_params, &nu_params, &metric, &prims, &state,
                (ghl_m1_direction_t)direction, &physical_number_flux),
          ghl_success, "public physical number flux", 1220 + direction);
    require_close(
          physical_number_flux,
          metric.lapse * number_flux[direction] - metric.betaU[direction] * state.N,
          2.0e-11, 2.0e-13, "physical number flux", 1220 + direction);
  }

  ghl_m1_neutrino_state thin_with_diagnostics;
  ghl_m1_neutrino_state thin_compatibility;
  ghl_m1_neutrino_exchange exchange_with_diagnostics;
  ghl_m1_neutrino_exchange compatibility_exchange;
  ghl_m1_neutrino_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&diagnostics);
  int thin_with_diagnostics_selected = 0;
  int compatibility_selected = 0;
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &rates, 0.1, 3.0, &state,
              &thin_with_diagnostics_selected, &thin_with_diagnostics,
              &exchange_with_diagnostics, &diagnostics),
        ghl_success, "diagnostic explicit thin update", 1230);
  require_condition(
        thin_with_diagnostics_selected == 1,
        "diagnostic explicit thin update was not selected", 1230);
  require_state_finite_and_admissible(
        m1_params, &nu_params, &metric, &thin_with_diagnostics, 1230);
  require_exchange_contract(
        &state, &thin_with_diagnostics, &exchange_with_diagnostics, &metric, 3.0, 1230);

  require_error(
        ghl_m1_try_neutrino_explicit_thin_update(
              m1_params, &nu_params, &metric, &prims, &rates, 0.1, 3.0, &state,
              &compatibility_selected, &thin_compatibility, &compatibility_exchange),
        ghl_success, "compatibility explicit thin update", 1231);
  require_condition(
        compatibility_selected == 1,
        "compatibility explicit thin update was not selected", 1231);
  require_state_finite_and_admissible(
        m1_params, &nu_params, &metric, &thin_compatibility, 1231);
  require_exchange_contract(
        &state, &thin_compatibility, &compatibility_exchange, &metric, 3.0, 1231);
  require_close(
        thin_compatibility.N, thin_with_diagnostics.N, 2.0e-11, 2.0e-13,
        "compatibility thin number", 1232);
  require_close(
        thin_compatibility.E, thin_with_diagnostics.E, 2.0e-11, 2.0e-13,
        "compatibility thin energy", 1232);

  /* Heavy-lepton number exchange may be nonzero, but it must not change
   * electron lepton number or the matter electron fraction. */
  ghl_m1_neutrino_rates heavy_rates;
  make_rates(ghl_m1_neutrino_nux, 0.08, 0.10, 0.10, 2.0, 2.0, 0.0, 0.0, &heavy_rates);
  ghl_m1_neutrino_state heavy_state = state;
  heavy_state.N = 0.5;
  ghl_m1_neutrino_state heavy_output = heavy_state;
  ghl_m1_neutrino_exchange heavy_exchange = { 0 };
  int heavy_selected = 0;
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &heavy_rates, 0.1, 3.0,
              &heavy_state, &heavy_selected, &heavy_output, &heavy_exchange, NULL),
        ghl_success, "heavy-lepton explicit thin update", 1233);
  require_condition(
        heavy_selected == 1, "heavy-lepton thin update was not selected", 1233);
  require_state_finite_and_admissible(
        m1_params, &nu_params, &metric, &heavy_output, 1233);
  require_exchange_contract(
        &heavy_state, &heavy_output, &heavy_exchange, &metric, 3.0, 1233);
  require_condition(
        heavy_exchange.dN_rad_total > 0.0 && heavy_exchange.dL_rad_cc == 0.0
              && heavy_exchange.dYe_matter == 0.0,
        "heavy-lepton number exchange changed electron lepton number", 1233);
}

static void test_source_dispatcher_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 0.08, 0.10, 0.10, 1.0, 2.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_state state;
  make_state(rng, &metric, &state);
  ghl_m1_neutrino_state state_out = state;
  ghl_m1_neutrino_exchange exchange = { 0 };
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics neutrino_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);

  /* Exercise each required output/state pointer independently. The outer
   * dispatcher guard must reject these before it reads any input object. */
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, NULL, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source state output", 1250);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, NULL, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source exchange output", 1251);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, NULL, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source diagnostics output", 1252);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, NULL),
        ghl_error_m1_null_pointer, "NULL source neutrino diagnostics", 1253);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, &state, NULL, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source transport state", 1254);

  /* These pointers are checked after the output transaction has been
   * initialized. Vary one input at a time so each short-circuit operand is
   * independently selected. */
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, NULL, &nu_params, &metric, &prims, &rates, &state, &state, 0.1, 3.0,
              &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source M1 parameters", 1255);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, NULL, &metric, &prims, &rates, &state, &state, 0.1, 3.0,
              &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source neutrino parameters", 1256);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, NULL, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source metric", 1257);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, NULL, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source primitives", 1258);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, NULL, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source rates", 1259);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, NULL, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_null_pointer, "NULL source input state", 1260);

  /* State and floor fields are validated in sequence. Each mutation selects
   * one finite/nonfinite predicate while preserving the preceding terms. */
  ghl_m1_neutrino_state bad_state = state;
  bad_state.N = NAN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, &bad_state, &state,
              0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite source number state", 1261);
  bad_state = state;
  bad_state.E = NAN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, &bad_state, &state,
              0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite source energy state", 1262);
  for(int direction = 0; direction < 3; ++direction) {
    bad_state = state;
    bad_state.F[direction] = NAN;
    require_error(
          ghl_m1_solve_neutrino_source_update(
                NULL, m1_params, &nu_params, &metric, &prims, &rates, &bad_state, &state,
                0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
          ghl_error_m1_invalid_state, "nonfinite source flux state", 1263 + direction);
  }
  bad_state = state;
  bad_state.N = 0.0;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, &bad_state, &state,
              0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "source number below floor", 1266);

  ghl_m1_neutrino_parameters bad_nu = nu_params;
  bad_nu.N_floor = NAN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &bad_nu, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite source N floor", 1267);
  bad_nu = nu_params;
  bad_nu.N_floor = -DBL_MIN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &bad_nu, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "negative source N floor", 1268);
  bad_nu = nu_params;
  bad_nu.J_floor = NAN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &bad_nu, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite source J floor", 1269);
  bad_nu = nu_params;
  bad_nu.J_floor = -DBL_MIN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &bad_nu, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "negative source J floor", 1270);
  bad_nu = nu_params;
  bad_nu.Gamma_N_floor = NAN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &bad_nu, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite source Gamma floor", 1271);
  bad_nu = nu_params;
  bad_nu.Gamma_N_floor = -DBL_MIN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &bad_nu, &metric, &prims, &rates, &state, &state, 0.1,
              3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "negative source Gamma floor", 1272);

  /* Validate every dispatcher option error before branch selection. */
  ghl_m1_neutrino_source_options options = branched_options();
  options.policy = (ghl_m1_neutrino_source_policy_t)99;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "invalid source policy", 1273);
  options = branched_options();
  options.ye_policy = (ghl_m1_neutrino_ye_policy_t)99;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "invalid source Ye policy", 1274);
  options = branched_options();
  options.thick_equilibrium_threshold = NAN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite thick threshold", 1275);
  options = branched_options();
  options.scattering_threshold = NAN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite scattering threshold", 1276);
  options = branched_options();
  options.thermalized_number_threshold = NAN;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite thermalized threshold", 1277);
  options = branched_options();
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              NAN, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite source timestep", 1278);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              -0.1, 3.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "negative source timestep", 1278);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              0.1, 0.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "zero source baryon normalization", 1279);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              0.1, NAN, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite source baryon normalization", 1279);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state,
              0.1, -1.0, &state_out, &exchange, &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "negative source baryon normalization", 1279);

  /* The dispatcher validates dt*alpha after the metric/configuration boundary,
   * and validates the post-transport state separately from state_input. */
  ghl_metric_quantities overflow_dt_metric = metric;
  overflow_dt_metric.lapse = DBL_MAX;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &overflow_dt_metric, &prims, &rates, &state,
              &state, DBL_MAX, 3.0, &state_out, &exchange, &diagnostics,
              &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "overflowed dispatcher dt alpha", 1282);
  require_condition(
        diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && neutrino_diagnostics.source_failures == 1,
        "dispatcher dt-alpha failure was not published", 1282);
  require_zero_exchange(&exchange, "dispatcher dt-alpha exchange", 1282);

  ghl_m1_neutrino_state bad_transport = state;
  bad_transport.F[1] = NAN;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &metric, &prims, &rates, &state,
              &bad_transport, 0.1, 3.0, &state_out, &exchange, &diagnostics,
              &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "invalid dispatcher transport state", 1283);
  require_condition(
        diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && neutrino_diagnostics.source_failures == 1,
        "dispatcher transport-state failure was not published", 1283);
  require_zero_exchange(&exchange, "dispatcher transport-state exchange", 1283);

  /* A frozen-primitives failure propagates from the closure/moments walk into
   * every dispatch route: the thin branch reaches it through
   * evaluate_closure_and_sources, whose own closure failure propagates. */
  {
    ghl_primitive_quantities bad_thin_prims = prims;
    bad_thin_prims.vU[0] = NAN;
    ghl_m1_neutrino_rates thin_rates;
    make_rates(ghl_m1_neutrino_nue, 0.08, 0.10, 0.10, 1.0, 2.0, 0.0, 0.0, &thin_rates);
    ghl_m1_neutrino_source_options thin_options = branched_options();
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
    require_error(
          ghl_m1_solve_neutrino_source_update(
                &thin_options, m1_params, &nu_params, &metric, &bad_thin_prims,
                &thin_rates, &state, &state, 0.1, 3.0, &state_out, &exchange,
                &diagnostics, &neutrino_diagnostics),
          ghl_error_m1_invalid_state, "dispatcher thin closure failure", 1283);
    require_condition(
          diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
                && neutrino_diagnostics.source_failures == 1,
          "dispatcher thin closure failure was not published", 1283);
    require_zero_exchange(&exchange, "dispatcher thin closure exchange", 1283);
  }

  /* Mean-energy bounds on the neutrino parameters reject the thin endpoint
   * through the shared bounds check. */
  {
    ghl_m1_neutrino_parameters bounded_nu = nu_params;
    bounded_nu.enforce_mean_energy_bounds = 1;
    bounded_nu.mean_energy_min = 1.0e6;
    ghl_m1_neutrino_rates thin_rates;
    make_rates(ghl_m1_neutrino_nue, 0.08, 0.10, 0.10, 1.0, 2.0, 0.0, 0.0, &thin_rates);
    ghl_m1_neutrino_source_options thin_options = branched_options();
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
    require_error(
          ghl_m1_solve_neutrino_source_update(
                &thin_options, m1_params, &bounded_nu, &metric, &prims, &thin_rates,
                &state, &state, 0.1, 3.0, &state_out, &exchange, &diagnostics,
                &neutrino_diagnostics),
          ghl_error_m1_invalid_state, "dispatcher thin endpoint bounds", 1283);
    require_condition(
          diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
                && neutrino_diagnostics.source_failures == 1,
          "dispatcher thin bounds failure was not published", 1283);
    require_zero_exchange(&exchange, "dispatcher thin bounds exchange", 1283);
  }

  /* A superluminal frozen coordinate velocity fails the Eulerian-velocity
   * walk inside the thick/scattering branch publication. */
  {
    ghl_primitive_quantities superluminal_prims = prims;
    superluminal_prims.vU[0] = 1.5;
    ghl_m1_neutrino_rates thick_rates;
    make_rates(ghl_m1_neutrino_nue, 0.08, 40.0, 0.10, 1.0, 2.0, 0.0, 0.0, &thick_rates);
    ghl_m1_neutrino_source_options thick_options = branched_options();
    thick_options.thick_equilibrium_threshold = 0.0;
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
    require_error(
          ghl_m1_solve_neutrino_source_update(
                &thick_options, m1_params, &nu_params, &metric, &superluminal_prims,
                &thick_rates, &state, &state, 1.0, 3.0, &state_out, &exchange,
                &diagnostics, &neutrino_diagnostics),
          ghl_error_u0_singular, "dispatcher thick eulerian velocity", 1283);
    require_condition(
          diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
                && neutrino_diagnostics.source_failures == 1,
          "dispatcher thick velocity failure was not published", 1283);
    require_zero_exchange(&exchange, "dispatcher thick velocity exchange", 1283);
  }

  /* The dispatcher validates the immutable metric at its own boundary after
   * the pointer checks, so a zero volume must fail there with transactional
   * publication before any state validation. */
  ghl_metric_quantities bad_dispatcher_metric = metric;
  bad_dispatcher_metric.sqrt_detgamma = 0.0;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              NULL, m1_params, &nu_params, &bad_dispatcher_metric, &prims, &rates,
              &state, &state, 0.1, 3.0, &state_out, &exchange, &diagnostics,
              &neutrino_diagnostics),
        ghl_error_m1_invalid_metric, "invalid dispatcher metric", 1284);
  require_condition(
        diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && neutrino_diagnostics.source_failures == 1,
        "dispatcher metric failure was not published", 1284);
  require_zero_exchange(&exchange, "dispatcher metric exchange", 1284);
  const double invalid_lapse_volume[] = { NAN, 0.0, -1.0, INFINITY };
  for(int variant = 0; variant < 4; ++variant) {
    bad_dispatcher_metric = metric;
    bad_dispatcher_metric.lapse = invalid_lapse_volume[variant];
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
    require_error(
          ghl_m1_solve_neutrino_source_update(
                NULL, m1_params, &nu_params, &bad_dispatcher_metric, &prims, &rates,
                &state, &state, 0.1, 3.0, &state_out, &exchange, &diagnostics,
                &neutrino_diagnostics),
          ghl_error_m1_invalid_metric, "invalid dispatcher lapse", 1284);
    require_condition(
          diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
                && neutrino_diagnostics.source_failures == 1,
          "dispatcher lapse failure was not published", 1284);
    require_zero_exchange(&exchange, "dispatcher lapse exchange", 1284);
    bad_dispatcher_metric = metric;
    bad_dispatcher_metric.sqrt_detgamma = invalid_lapse_volume[variant];
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
    require_error(
          ghl_m1_solve_neutrino_source_update(
                NULL, m1_params, &nu_params, &bad_dispatcher_metric, &prims, &rates,
                &state, &state, 0.1, 3.0, &state_out, &exchange, &diagnostics,
                &neutrino_diagnostics),
          ghl_error_m1_invalid_metric, "invalid dispatcher volume", 1284);
    require_condition(
          diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
                && neutrino_diagnostics.source_failures == 1,
          "dispatcher volume failure was not published", 1284);
    require_zero_exchange(&exchange, "dispatcher volume exchange", 1284);
  }
}

static void test_source_compatibility_and_selector_failures(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 0.08, 0.10, 0.10, 1.0, 2.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_state state;
  make_state(rng, &metric, &state);
  ghl_m1_neutrino_state state_out;
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics neutrino_diagnostics;

  /* Compatibility mode may reject closure fallback, so exercise the
   * successful pre-closure check with fallback explicitly disabled. The
   * chosen opacity is thin and the ordinary converged closure is the source
   * of the endpoint. */
  ghl_m1_neutrino_source_options strict_options = branched_options();
  strict_options.allow_closure_fallback = false;
  ghl_m1_neutrino_state strict_out;
  ghl_m1_neutrino_exchange strict_exchange;
  ghl_m1_neutrino_source_diagnostics strict_diagnostics;
  ghl_m1_neutrino_diagnostics strict_nd;
  ghl_m1_neutrino_diagnostics_initialize(&strict_nd);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &strict_options, m1_params, &nu_params, &metric, &prims, &rates, &state,
              &state, 0.1, 3.0, &strict_out, &strict_exchange, &strict_diagnostics,
              &strict_nd),
        ghl_success, "strict compatibility source update", 1280);
  require_condition(
        strict_diagnostics.path == ghl_m1_neutrino_source_path_thin_explicit
              && !strict_diagnostics.closure_fallback_used
              && !strict_diagnostics.terminal_no_update,
        "strict compatibility path was not accepted", 1280);
  require_state_finite_and_admissible(m1_params, &nu_params, &metric, &strict_out, 1280);
  require_exchange_contract(&state, &strict_out, &strict_exchange, &metric, 3.0, 1280);

  /* Compatibility mode performs a strict pre-closure check when fallback is
   * disabled. Exercise both the closure error and the documented fallback
   * denial before branch selection. */
  ghl_m1_neutrino_source_options strict_failure_options = branched_options();
  strict_failure_options.allow_closure_fallback = false;
  ghl_primitive_quantities bad_preclosure_prims = prims;
  bad_preclosure_prims.vU[0] = NAN;
  ghl_m1_neutrino_diagnostics preclosure_nd;
  ghl_m1_neutrino_diagnostics_initialize(&preclosure_nd);
  ghl_m1_neutrino_state preclosure_out = state;
  ghl_m1_neutrino_exchange preclosure_exchange;
  ghl_m1_neutrino_source_diagnostics preclosure_diagnostics;
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &strict_failure_options, m1_params, &nu_params, &metric,
              &bad_preclosure_prims, &rates, &state, &state, 0.1, 3.0, &preclosure_out,
              &preclosure_exchange, &preclosure_diagnostics, &preclosure_nd),
        ghl_error_m1_invalid_state, "compatibility pre-closure error", 1283);
  require_condition(
        preclosure_diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && preclosure_nd.source_failures == 1,
        "compatibility pre-closure error was not published", 1283);
  require_zero_exchange(
        &preclosure_exchange, "compatibility pre-closure error exchange", 1283);

  ghl_primitive_quantities moving_prims = prims;
  moving_prims.vU[0] = 0.4;
  const ghl_m1_neutrino_state zero_flux_state
        = { .N = 1.0, .E = 1.0, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_neutrino_diagnostics_initialize(&preclosure_nd);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &strict_failure_options, m1_params, &nu_params, &metric, &moving_prims,
              &rates, &zero_flux_state, &zero_flux_state, 0.1, 3.0, &preclosure_out,
              &preclosure_exchange, &preclosure_diagnostics, &preclosure_nd),
        ghl_error_m1_implicit_solve_failure, "compatibility pre-closure fallback denial",
        1284);
  require_condition(
        preclosure_diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && preclosure_diagnostics.closure_fallback_used
              && preclosure_nd.source_failures == 1,
        "compatibility fallback denial was not published", 1284);
  require_condition(
        memcmp(&preclosure_out, &zero_flux_state, sizeof(preclosure_out)) == 0,
        "compatibility fallback denial changed state", 1284);
  require_zero_exchange(
        &preclosure_exchange, "compatibility fallback denial exchange", 1284);

  /* Choose a finite opacity that is neither thin nor selected by either
   * shortcut threshold, so the compatibility policy reaches general implicit.
   */
  ghl_m1_neutrino_rates general_rates;
  make_rates(ghl_m1_neutrino_nue, 0.10, 2.0, 0.0, 1.0, 2.0, 0.0, 0.0, &general_rates);
  ghl_m1_neutrino_source_options general_options = branched_options();
  general_options.thick_equilibrium_threshold = DBL_MAX;
  general_options.scattering_threshold = DBL_MAX;
  ghl_m1_parameters terminal_params = *m1_params;
  terminal_params.newton_max_iterations = 1;
  terminal_params.newton_tolerance = DBL_MIN;
  terminal_params.newton_absolute_tolerance = DBL_MIN;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &general_options, &terminal_params, &nu_params, &metric, &prims,
              &general_rates, &state, &state, 1.0, 3.0, &state_out, &exchange,
              &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_implicit_terminal_fallback, "compatibility general terminal route",
        1285);
  require_condition(
        diagnostics.path == ghl_m1_neutrino_source_path_terminal_no_update
              && diagnostics.terminal_no_update,
        "compatibility general terminal status", 1285);
  require_condition(
        memcmp(&state_out, &state, sizeof(state_out)) == 0,
        "compatibility general terminal changed state", 1285);
  require_zero_exchange(&exchange, "compatibility general terminal exchange", 1285);

  /* The same general route must convert a non-retryable implicit hard failure
   * into its transactional dispatcher publication. */
  ghl_primitive_quantities bad_general_prims = prims;
  bad_general_prims.vU[0] = NAN;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &general_options, m1_params, &nu_params, &metric, &bad_general_prims,
              &general_rates, &state, &state, 1.0, 3.0, &state_out, &exchange,
              &diagnostics, &neutrino_diagnostics),
        ghl_error_m1_invalid_state, "compatibility general hard failure", 1286);
  require_condition(
        diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && !diagnostics.terminal_no_update,
        "compatibility general hard-failure status", 1286);
  require_condition(
        memcmp(&state_out, &state, sizeof(state_out)) == 0,
        "compatibility general hard failure changed state", 1286);
  require_zero_exchange(&exchange, "compatibility general hard-failure exchange", 1286);
}

static void test_pair_source_conservation(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  for(int case_index = 0; case_index < PAIR_RANDOM_CASES; ++case_index) {
    ghl_metric_quantities metric;
    make_metric(rng, m1_test_rng_between(rng, 0.75, 1.25), &metric);
    ghl_primitive_quantities prims;
    make_primitives(&prims);
    ghl_m1_neutrino_parameters nu_params[2];
    make_neutrino_parameters(&nu_params[0]);
    make_neutrino_parameters(&nu_params[1]);
    ghl_m1_neutrino_rates active_rates[2], base_rates[2];
    make_rates(
          ghl_m1_neutrino_nue, 0.08, 0.12, 0.10, 1.0, 2.0, 0.02, 0.03, &active_rates[0]);
    make_rates(
          ghl_m1_neutrino_anue, 0.08, 0.12, 0.10, 1.0, 2.1, 0.02, 0.04,
          &active_rates[1]);
    base_rates[0] = active_rates[0];
    base_rates[1] = active_rates[1];
    for(int species = 0; species < 2; ++species) {
      for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
        base_rates[species].eta_N_pair[process] = 0.0;
        base_rates[species].eta_E_pair[process] = 0.0;
      }
    }

    ghl_m1_neutrino_state state_input[2], state_transport[2];
    for(int species = 0; species < 2; ++species) {
      make_state(rng, &metric, &state_transport[species]);
      state_input[species] = state_transport[species];
      state_input[species].E *= 0.8;
      state_input[species].N *= 1.1;
      for(int i = 0; i < 3; ++i) {
        state_input[species].F[i] *= 0.8;
      }
    }
    ghl_m1_neutrino_state active_out[2], base_out[2];
    ghl_m1_neutrino_exchange active_exchange[2], base_exchange[2];
    ghl_m1_neutrino_source_diagnostics active_diag[2], base_diag[2];
    ghl_m1_neutrino_diagnostics active_nd[2], base_nd[2];
    for(int species = 0; species < 2; ++species) {
      ghl_m1_neutrino_diagnostics_initialize(&active_nd[species]);
      ghl_m1_neutrino_diagnostics_initialize(&base_nd[species]);
    }
    ghl_error_codes_t error = ghl_m1_solve_neutrino_pair_source_update(
          m1_params, nu_params, &metric, &prims, active_rates, state_input,
          state_transport, 0.25, 3.0, active_out, active_exchange, active_diag,
          active_nd);
    require_error(error, ghl_success, "active pair source update", case_index);
    error = ghl_m1_solve_neutrino_pair_source_update(
          m1_params, nu_params, &metric, &prims, base_rates, state_input,
          state_transport, 0.25, 3.0, base_out, base_exchange, base_diag, base_nd);
    require_error(error, ghl_success, "base pair source update", case_index);

    const double pair_delta_nue = active_out[0].N - base_out[0].N;
    const double pair_delta_anue = active_out[1].N - base_out[1].N;
    require_close(
          pair_delta_nue, pair_delta_anue, 2.0e-9, 2.0e-12,
          "shared pair number increment", case_index);
    for(int species = 0; species < 2; ++species) {
      require_state_finite_and_admissible(
            m1_params, &nu_params[species], &metric, &active_out[species],
            2000 + case_index);
      require_exchange_contract(
            &state_transport[species], &active_out[species], &active_exchange[species],
            &metric, 3.0, 2000 + case_index);
      require_close(
            active_exchange[species].dYe_matter, base_exchange[species].dYe_matter,
            2.0e-9, 2.0e-12, "pair-independent electron-fraction exchange",
            2000 + case_index);
      require_close(
            active_exchange[species].dL_rad_cc, base_exchange[species].dL_rad_cc, 2.0e-9,
            2.0e-12, "pair-independent charged-current exchange", 2000 + case_index);
      require_condition(
            active_nd[species].source_converged == base_nd[species].source_converged + 1,
            "active pair did not add one E/F convergence", 2000 + case_index);
    }
  }
}

static void
test_pair_effective_opacity_overflow(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params[2];
  ghl_m1_neutrino_rates rates[2];
  const ghl_m1_neutrino_state state[2]
        = { { .N = 1.0, .E = 1.0 }, { .N = 1.0, .E = 1.0 } };
  for(int species = 0; species < 2; ++species) {
    make_neutrino_parameters(&nu_params[species]);
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0,
          DBL_MIN, 1.0, 0.0, 1.0, &rates[species]);
    require_error(
          ghl_m1_validate_neutrino_rates(&rates[species], NULL), ghl_success,
          "pair effective-opacity overflow input rates", 2101);
  }
  ghl_m1_neutrino_state output[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics nd[2] = { { 0 }, { 0 } };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state, state, 1.0, 1.0,
              output, exchange, diagnostics, nd),
        ghl_error_m1_invalid_state, "unrepresentable pair effective opacity", 2101);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&output[species], &state[species], sizeof(output[species])) == 0
                && diagnostics[species].path == ghl_m1_neutrino_source_path_hard_failure
                && nd[species].source_failures == 1,
          "pair effective-opacity failure was not transactional", 2101);
    require_zero_exchange(
          &exchange[species], "pair effective-opacity overflow exchange", 2101);
  }

  /* The direct effective opacity can also overflow before the scaled
   * fallback, which then succeeds with a representable quotient. */
  for(int species = 0; species < 2; ++species) {
    make_neutrino_parameters(&nu_params[species]);
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 1.0,
          1.0, 0.0, DBL_MAX, &rates[species]);
    require_error(
          ghl_m1_validate_neutrino_rates(&rates[species], NULL), ghl_success,
          "pair effective-opacity scaled-fallback input rates", 2102);
  }
  for(int species = 0; species < 2; ++species) {
    nd[species] = (ghl_m1_neutrino_diagnostics){ 0 };
  }
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state, state, 1.0, 1.0,
              output, exchange, diagnostics, nd),
        ghl_error_m1_invalid_state, "overflowing pair effective opacity fallback", 2102);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&output[species], &state[species], sizeof(output[species])) == 0
                && nd[species].source_failures == 1,
          "pair scaled-opacity failure was not transactional", 2102);
    require_zero_exchange(&exchange[species], "pair scaled-opacity exchange", 2102);
  }

  /* A subnormal-tiny energy emissivity over an enormous equilibrium energy
   * sends the direct quotient to zero and the scaled quotient below the
   * representable range: the pair operation must fail transactionally. */
  for(int species = 0; species < 2; ++species) {
    make_neutrino_parameters(&nu_params[species]);
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0,
          DBL_MAX, 1.0, 0.0, nextafter(0.0, 1.0), &rates[species]);
    require_error(
          ghl_m1_validate_neutrino_rates(&rates[species], NULL), ghl_success,
          "pair effective-opacity underflow input rates", 2103);
  }
  for(int species = 0; species < 2; ++species) {
    nd[species] = (ghl_m1_neutrino_diagnostics){ 0 };
  }
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state, state, 1.0, 1.0,
              output, exchange, diagnostics, nd),
        ghl_error_m1_invalid_state, "underflowed pair effective opacity", 2103);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&output[species], &state[species], sizeof(output[species])) == 0
                && nd[species].source_failures == 1,
          "pair underflow-opacity failure was not transactional", 2103);
    require_zero_exchange(&exchange[species], "pair underflow-opacity exchange", 2103);
  }
}

/* Independent oracle for the pair-only case below.  With a flat metric and
 * zero fluid velocity, Gamma_N=1 and the model in PAIR_SOURCE_MODEL reduces
 * to
 *
 *   d = h q [1 - (N_e+d)(N_a+d)/(n_eq_e n_eq_a)]
 *
 * for the shared number increment.  The E/F equations are then scalar
 * backward-Euler updates with
 *   kappa_s = eta_E,s/J_eq,s * (N_partner+d)/n_eq,partner.
 * This deliberately does not call any GRHayL source or closure helper. */
static double pair_oracle_number_root(
      const double h,
      const double q,
      const double N_e,
      const double N_a,
      const double n_eq_e,
      const double n_eq_a) {
  const long double H = (long double)h * (long double)q;
  const long double D = (long double)n_eq_e * (long double)n_eq_a;
  const long double B = D / H + (long double)N_e + (long double)N_a;
  /* Preserve an exact zero at equilibrium even when multiply-add is fused. */
  const long double delta_e = (long double)N_e - (long double)n_eq_e;
  const long double delta_a = (long double)N_a - (long double)n_eq_a;
  const long double C = delta_e * (long double)N_a + (long double)n_eq_e * delta_a;
  const long double discriminant = B * B - 4.0L * C;
  require_condition(
        isfinite((double)discriminant) && discriminant >= 0.0L,
        "pair oracle quadratic is invalid", 2050);
  /* This form retains the root continuous with d=0 as h approaches zero. */
  const long double d = -2.0L * C / (B + sqrtl(discriminant));
  require_condition(isfinite((double)d), "pair oracle extent is nonfinite", 2050);
  return (double)d;
}

static void make_pair_oracle_rates(
      ghl_m1_neutrino_rates rates[2],
      const double n_eq[2],
      const double mean_energy[2],
      const double eta_N_pair[3],
      const double eta_E_pair[2][3]) {
  make_rates(
        ghl_m1_neutrino_nue, 0.0, 0.0, 0.0, n_eq[0], mean_energy[0], 0.0, 0.0,
        &rates[0]);
  make_rates(
        ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, n_eq[1], mean_energy[1], 0.0, 0.0,
        &rates[1]);
  for(int process = 0; process < 3; ++process) {
    rates[0].eta_N_pair[process] = eta_N_pair[process];
    rates[1].eta_N_pair[process] = eta_N_pair[process];
    rates[0].eta_E_pair[process] = eta_E_pair[0][process];
    rates[1].eta_E_pair[process] = eta_E_pair[1][process];
  }
}

static void
test_pair_source_independent_oracle(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric = { 0 };
  metric.lapse = 1.0;
  metric.lapseinv = 1.0;
  metric.lapseinv2 = 1.0;
  metric.detgamma = 1.0;
  metric.sqrt_detgamma = 1.0;
  for(int direction = 0; direction < 3; ++direction) {
    metric.gammaDD[direction][direction] = 1.0;
    metric.gammaUU[direction][direction] = 1.0;
  }
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params[2];
  make_neutrino_parameters(&nu_params[0]);
  make_neutrino_parameters(&nu_params[1]);

  const double n_eq[2] = { 0.8, 1.4 };
  const double mean_energy[2] = { 2.0, 3.0 };
  const double eta_N_pair[3] = { 0.07, 0.02, 0.01 };
  const double eta_E_pair[2][3] = { { 0.20, 0.08, 0.03 }, { 0.35, 0.12, 0.05 } };
  ghl_m1_neutrino_rates rates[2];
  make_pair_oracle_rates(rates, n_eq, mean_energy, eta_N_pair, eta_E_pair);

  const double dt = 0.25;
  const double state_number[2] = { 0.35, 1.20 };
  ghl_m1_neutrino_state state_transport[2]
        = { { .N = state_number[0], .E = 1.7, .F = { 0.30, -0.20, 0.10 } },
            { .N = state_number[1], .E = 1.3, .F = { -0.25, 0.15, 0.05 } } };
  const ghl_m1_neutrino_state state_input[2]
        = { state_transport[0], state_transport[1] };
  ghl_m1_neutrino_state state_out[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics neutrino_diagnostics[2];
  for(int species = 0; species < 2; ++species) {
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics[species]);
  }

  const ghl_error_codes_t error = ghl_m1_solve_neutrino_pair_source_update(
        m1_params, nu_params, &metric, &prims, rates, state_input, state_transport, dt,
        3.0, state_out, exchange, diagnostics, neutrino_diagnostics);
  require_error(error, ghl_success, "pair independent oracle update", 2050);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          diagnostics[species].implicit.fallback_substeps == 1
                && !diagnostics[species].implicit.used_fallback_substepping
                && (diagnostics[species].implicit.solution_path_flags
                    & ghl_m1_solution_path_primary_convergence)
                         != 0u
                && (diagnostics[species].implicit.solution_path_flags
                    & ghl_m1_solution_path_substepping)
                         == 0u,
          "pair emission oracle did not use one substep", 2050);
  }

  double q = 0.0;
  for(int process = 0; process < 3; ++process) {
    q += eta_N_pair[process];
  }
  const double d = pair_oracle_number_root(
        dt, q, state_number[0], state_number[1], n_eq[0], n_eq[1]);
  const double expected_number[2] = { state_number[0] + d, state_number[1] + d };
  for(int species = 0; species < 2; ++species) {
    require_close(
          state_out[species].N, expected_number[species], 2.0e-9, 2.0e-12,
          "pair oracle number", 2051 + species);
    const int partner = 1 - species;
    double eta_E = 0.0;
    for(int process = 0; process < 3; ++process) {
      eta_E += eta_E_pair[species][process];
    }
    const double partner_ratio = expected_number[partner] / n_eq[partner];
    const double kappa_E = eta_E / rates[species].J_eq * partner_ratio;
    const double denominator = 1.0 + dt * kappa_E;
    const double expected_E = (state_transport[species].E + dt * eta_E) / denominator;
    require_close(
          state_out[species].E, expected_E, 2.0e-9, 2.0e-12, "pair oracle energy",
          2053 + species);
    for(int direction = 0; direction < 3; ++direction) {
      const double expected_F = state_transport[species].F[direction] / denominator;
      require_close(
            state_out[species].F[direction], expected_F, 2.0e-9, 2.0e-12,
            "pair oracle flux", 2055 + species * 3 + direction);
    }
    require_exchange_contract(
          &state_transport[species], &state_out[species], &exchange[species], &metric,
          3.0, 2060 + species);
  }

  /* Pair energy emission can be active while the shared number reaction is
   * exactly zero. This reaches the q==0 source branch and proves that the
   * partner-weighted E/F damping still uses the unchanged partner number. */
  ghl_m1_neutrino_rates zero_number_rates[2] = { rates[0], rates[1] };
  for(int species = 0; species < 2; ++species) {
    for(int process = 0; process < 3; ++process) {
      zero_number_rates[species].eta_N_pair[process] = 0.0;
    }
  }
  ghl_m1_neutrino_state zero_number_out[2];
  ghl_m1_neutrino_exchange zero_number_exchange[2];
  ghl_m1_neutrino_source_diagnostics zero_number_diagnostics[2];
  ghl_m1_neutrino_diagnostics zero_number_nd[2];
  for(int species = 0; species < 2; ++species) {
    ghl_m1_neutrino_diagnostics_initialize(&zero_number_nd[species]);
  }
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, zero_number_rates, state_input,
              state_transport, dt, 3.0, zero_number_out, zero_number_exchange,
              zero_number_diagnostics, zero_number_nd),
        ghl_success, "pair zero-number source update", 2065);
  for(int species = 0; species < 2; ++species) {
    const int partner = 1 - species;
    double eta_E = 0.0;
    for(int process = 0; process < 3; ++process) {
      eta_E += eta_E_pair[species][process];
    }
    const double partner_ratio = state_number[partner] / n_eq[partner];
    const double denominator = 1.0 + dt * eta_E / rates[species].J_eq * partner_ratio;
    require_close(
          zero_number_out[species].N, state_number[species], 0.0, 0.0,
          "pair zero-number extent", 2066 + species);
    require_close(
          zero_number_out[species].E,
          (state_transport[species].E + dt * eta_E) / denominator, 2.0e-9, 2.0e-12,
          "pair zero-number energy", 2068 + species);
    for(int direction = 0; direction < 3; ++direction) {
      require_close(
            zero_number_out[species].F[direction],
            state_transport[species].F[direction] / denominator, 2.0e-9, 2.0e-12,
            "pair zero-number flux", 2070 + species * 3 + direction);
    }
  }

  ghl_m1_neutrino_state zero_step_out[2];
  ghl_m1_neutrino_exchange zero_step_exchange[2];
  ghl_m1_neutrino_source_diagnostics zero_step_diagnostics[2];
  ghl_m1_neutrino_diagnostics zero_step_nd[2];
  for(int species = 0; species < 2; ++species) {
    ghl_m1_neutrino_diagnostics_initialize(&zero_step_nd[species]);
  }
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, state_transport,
              0.0, 3.0, zero_step_out, zero_step_exchange, zero_step_diagnostics,
              zero_step_nd),
        ghl_success, "pair zero-timestep source update", 2078);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(
                &zero_step_out[species], &state_transport[species],
                sizeof(zero_step_out[species]))
                == 0,
          "pair zero-timestep changed state", 2078);
    require_zero_exchange(
          &zero_step_exchange[species], "pair zero-timestep exchange", 2078);
  }

  /* The same independent pair model is also checked on net absorption and
   * exact equilibrium inputs.  These cases retain nonzero flux for the
   * absorption endpoint, and require the production schedule to accept the
   * documented first (single-substep) attempt. */
  const ghl_m1_neutrino_state additional_transport[2][2]
        = { { { .N = 1.0, .E = 1.9, .F = { 0.18, -0.11, 0.07 } },
              { .N = 1.3, .E = 1.8, .F = { -0.14, 0.09, 0.05 } } },
            { { .N = n_eq[0], .E = rates[0].J_eq, .F = { 0.0, 0.0, 0.0 } },
              { .N = n_eq[1], .E = rates[1].J_eq, .F = { 0.0, 0.0, 0.0 } } } };
  const char *const additional_labels[2]
        = { "pair oracle net absorption", "pair oracle exact equilibrium" };
  for(int additional_case = 0; additional_case < 2; ++additional_case) {
    ghl_m1_neutrino_state additional_input[2]
          = { additional_transport[additional_case][0],
              additional_transport[additional_case][1] };
    ghl_m1_neutrino_state additional_out[2];
    ghl_m1_neutrino_exchange additional_exchange[2];
    ghl_m1_neutrino_source_diagnostics additional_diagnostics[2];
    ghl_m1_neutrino_diagnostics additional_nd[2];
    for(int species = 0; species < 2; ++species) {
      ghl_m1_neutrino_diagnostics_initialize(&additional_nd[species]);
    }
    require_error(
          ghl_m1_solve_neutrino_pair_source_update(
                m1_params, nu_params, &metric, &prims, rates, additional_input,
                additional_transport[additional_case], dt, 3.0, additional_out,
                additional_exchange, additional_diagnostics, additional_nd),
          ghl_success, additional_labels[additional_case], 2080 + additional_case);
    for(int species = 0; species < 2; ++species) {
      require_condition(
            additional_diagnostics[species].implicit.fallback_substeps == 1
                  && !additional_diagnostics[species].implicit.used_fallback_substepping
                  && (additional_diagnostics[species].implicit.solution_path_flags
                      & ghl_m1_solution_path_primary_convergence)
                           != 0u
                  && (additional_diagnostics[species].implicit.solution_path_flags
                      & ghl_m1_solution_path_substepping)
                           == 0u,
            "pair oracle did not use one substep", 2082 + additional_case);
      const int partner = 1 - species;
      double q_case = 0.0;
      double eta_E = 0.0;
      for(int process = 0; process < 3; ++process) {
        q_case += eta_N_pair[process];
        eta_E += eta_E_pair[species][process];
      }
      const double d_case = pair_oracle_number_root(
            dt, q_case, additional_transport[additional_case][0].N,
            additional_transport[additional_case][1].N, n_eq[0], n_eq[1]);
      if(additional_case == 0) {
        require_condition(d_case < 0.0, "pair absorption oracle did not absorb", 2083);
      }
      else {
        require_close(d_case, 0.0, 0.0, 0.0, "pair equilibrium oracle extent", 2083);
      }
      const double expected_number
            = additional_transport[additional_case][species].N + d_case;
      const double partner_number
            = additional_transport[additional_case][partner].N + d_case;
      const double denominator
            = 1.0 + dt * eta_E / rates[species].J_eq * partner_number / n_eq[partner];
      const double expected_energy
            = (additional_transport[additional_case][species].E + dt * eta_E)
              / denominator;
      require_close(
            additional_out[species].N, expected_number, 2.0e-9, 2.0e-12,
            "pair additional oracle number", 2084 + additional_case * 20 + species);
      require_close(
            additional_out[species].E, expected_energy, 2.0e-9, 2.0e-12,
            "pair additional oracle energy", 2086 + additional_case * 20 + species);
      for(int direction = 0; direction < 3; ++direction) {
        require_close(
              additional_out[species].F[direction],
              additional_transport[additional_case][species].F[direction] / denominator,
              2.0e-9, 2.0e-12, "pair additional oracle flux",
              2088 + additional_case * 20 + species * 3 + direction);
      }
      require_exchange_contract(
            &additional_transport[additional_case][species], &additional_out[species],
            &additional_exchange[species], &metric, 3.0,
            2092 + additional_case * 20 + species);
      if(additional_case == 1) {
        require_close(
              additional_out[species].N, n_eq[species], 2.0e-9, 2.0e-12,
              "pair equilibrium number", 2094 + species);
        require_close(
              additional_out[species].E, rates[species].J_eq, 2.0e-9, 2.0e-12,
              "pair equilibrium energy", 2096 + species);
        for(int direction = 0; direction < 3; ++direction) {
          require_close(
                additional_out[species].F[direction], 0.0, 0.0, 0.0,
                "pair equilibrium flux", 2098 + species * 3 + direction);
        }
      }
    }
  }
}

static void test_pair_inactive_delegation(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params[2];
  make_neutrino_parameters(&nu_params[0]);
  make_neutrino_parameters(&nu_params[1]);
  ghl_m1_neutrino_rates rates[2];
  make_rates(ghl_m1_neutrino_nue, 0.08, 0.12, 0.10, 1.0, 2.0, 0.0, 0.0, &rates[0]);
  make_rates(ghl_m1_neutrino_anue, 0.08, 0.12, 0.10, 1.1, 2.1, 0.0, 0.0, &rates[1]);

  ghl_m1_neutrino_state state_input[2], state_transport[2];
  for(int species = 0; species < 2; ++species) {
    make_state(rng, &metric, &state_transport[species]);
    state_input[species] = state_transport[species];
    state_input[species].N *= 1.15;
    state_input[species].E *= 0.83;
    for(int direction = 0; direction < 3; ++direction) {
      state_input[species].F[direction] *= 0.83;
    }
  }

  ghl_m1_neutrino_state pair_out[2], single_out[2];
  ghl_m1_neutrino_exchange pair_exchange[2], single_exchange[2];
  ghl_m1_neutrino_source_diagnostics pair_diagnostics[2], single_diagnostics[2];
  ghl_m1_neutrino_diagnostics pair_nd[2], single_nd[2];
  for(int species = 0; species < 2; ++species) {
    ghl_m1_neutrino_diagnostics_initialize(&pair_nd[species]);
    ghl_m1_neutrino_diagnostics_initialize(&single_nd[species]);
  }

  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, state_transport,
              0.25, 3.0, pair_out, pair_exchange, pair_diagnostics, pair_nd),
        ghl_success, "inactive pair source update", 2120);
  for(int species = 0; species < 2; ++species) {
    require_error(
          ghl_m1_solve_neutrino_source_update(
                NULL, m1_params, &nu_params[species], &metric, &prims, &rates[species],
                &state_input[species], &state_transport[species], 0.25, 3.0,
                &single_out[species], &single_exchange[species],
                &single_diagnostics[species], &single_nd[species]),
          ghl_success, "inactive pair single-species reference", 2121 + species);
    require_condition(
          pair_diagnostics[species].path == ghl_m1_neutrino_source_path_general_implicit
                && !pair_diagnostics[species].terminal_no_update,
          "inactive pair did not publish general path", 2123 + species);
    require_condition(
          memcmp(&pair_out[species], &single_out[species], sizeof(pair_out[species]))
                == 0,
          "inactive pair changed the single-species state", 2125 + species);
    require_condition(
          memcmp(
                &pair_exchange[species], &single_exchange[species],
                sizeof(pair_exchange[species]))
                == 0,
          "inactive pair changed the single-species exchange", 2127 + species);
    require_condition(
          pair_nd[species].source_converged == single_nd[species].source_converged
                && pair_nd[species].source_failures
                         == single_nd[species].source_failures,
          "inactive pair changed source diagnostics", 2129 + species);
  }
}

static void
test_pair_final_mean_energy_bounds(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric = { 0 };
  metric.lapse = metric.lapseinv = metric.lapseinv2 = 1.0;
  metric.detgamma = metric.sqrt_detgamma = 1.0;
  for(int direction = 0; direction < 3; ++direction) {
    metric.gammaDD[direction][direction] = 1.0;
    metric.gammaUU[direction][direction] = 1.0;
  }
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  const ghl_m1_neutrino_state state[2]
        = { { .N = 0.1, .E = 0.18 }, { .N = 0.1, .E = 0.18 } };
  const double n_eq[2] = { 1.0, 1.0 };
  const double mean_energy[2] = { 2.0, 2.0 };
  const double eta_N_pair[3] = { 100.0, 0.0, 0.0 };
  const double eta_E_pair[2][3] = { { 10.0, 0.0, 0.0 }, { 10.0, 0.0, 0.0 } };
  ghl_m1_neutrino_rates rates[2];
  make_pair_oracle_rates(rates, n_eq, mean_energy, eta_N_pair, eta_E_pair);
  for(int species = 0; species < 2; ++species) {
    rates[species].kappa_a_E = 10.0;
    rates[species].kappa_tr = 10.0;
    rates[species].eta_E = 20.0;
  }

  /* At rest, the independent stage leaves N fixed while E obeys scalar
   * backward Euler. Its mean exceeds the bound, but the pair stage brings
   * the final mean back inside. Both initial means already satisfy it. */
  const double independent_E
        = (state[0].E + rates[0].eta_E) / (1.0 + rates[0].kappa_a_E);
  const double expected_N
        = state[0].N
          + pair_oracle_number_root(
                1.0, eta_N_pair[0], state[0].N, state[1].N, n_eq[0], n_eq[1]);
  const double expected_E
        = (independent_E + eta_E_pair[0][0])
          / (1.0 + eta_E_pair[0][0] / rates[0].J_eq * expected_N / n_eq[1]);
  require_condition(
        state[0].E / state[0].N > 1.5 && state[0].E / state[0].N < 1.9
              && independent_E / state[0].N > 2.1 && expected_E / expected_N > 1.9
              && expected_E / expected_N < 2.1,
        "pair endpoint bound fixture", 2130);

  /* Accept the valid final endpoint; reject a tighter final bound; preserve
   * endpoint checks when the pair stage is inactive. */
  const double upper_bounds[] = { 2.1, 1.9, 2.1 };
  ghl_m1_neutrino_source_diagnostics expected_hard_diagnostics = { 0 };
  ghl_m1_initialize_implicit_solve_diagnostics(&expected_hard_diagnostics.implicit);
  expected_hard_diagnostics.path = ghl_m1_neutrino_source_path_hard_failure;
  for(int test_case = 0; test_case < 3; ++test_case) {
    ghl_m1_neutrino_parameters nu_params[2];
    for(int species = 0; species < 2; ++species) {
      make_neutrino_parameters(&nu_params[species]);
      nu_params[species].enforce_mean_energy_bounds = 1;
      nu_params[species].mean_energy_min = 1.5;
      nu_params[species].mean_energy_max = upper_bounds[test_case];
      if(test_case == 2) {
        rates[species].eta_N_pair[0] = 0.0;
        rates[species].eta_E_pair[0] = 0.0;
      }
    }
    const ghl_m1_neutrino_parameters saved_params[2] = { nu_params[0], nu_params[1] };
    ghl_m1_neutrino_state output[2];
    ghl_m1_neutrino_exchange exchange[2];
    ghl_m1_neutrino_source_diagnostics diagnostics[2];
    ghl_m1_neutrino_diagnostics neutrino_diagnostics[2] = { { 0 }, { 0 } };
    for(int species = 0; species < 2; ++species) {
      diagnostics[species] = (ghl_m1_neutrino_source_diagnostics){
        .path = ghl_m1_neutrino_source_path_thin_explicit,
        .closure_fallback_used = true,
        .terminal_no_update = true,
        .implicit = { .newton_iterations = 13 + species,
                      .line_search_backtracks = 14 + species,
                      .fallback_substeps = 15 + species,
                      .used_fallback_substepping = true,
                      .residual_max_norm = 16.0 + species,
                      .residual_scaled_norm = 17.0 + species,
                      .solution_path_flags = ghl_m1_solution_path_projection }
      };
    }
    const ghl_m1_neutrino_diagnostics neutrino_diagnostics_before[2]
          = { neutrino_diagnostics[0], neutrino_diagnostics[1] };
    require_error(
          ghl_m1_solve_neutrino_pair_source_update(
                m1_params, nu_params, &metric, &prims, rates, state, state, 1.0, 3.0,
                output, exchange, diagnostics, neutrino_diagnostics),
          test_case == 0 ? ghl_success : ghl_error_m1_invalid_state,
          "pair final mean-energy bounds", 2131 + test_case);
    require_condition(
          memcmp(nu_params, saved_params, sizeof(nu_params)) == 0,
          "pair solve changed caller parameters", 2134 + test_case);
    for(int species = 0; species < 2; ++species) {
      if(test_case == 0) {
        require_close(
              output[species].N, expected_N, 2.0e-9, 2.0e-12,
              "bounded pair endpoint number", 2137 + species);
        require_close(
              output[species].E, expected_E, 2.0e-9, 2.0e-12,
              "bounded pair endpoint energy", 2139 + species);
      }
      else {
        require_condition(
              memcmp(&output[species], &state[species], sizeof(output[species])) == 0,
              "pair bound rejection changed state", 2141 + species);
        require_zero_exchange(
              &exchange[species], "pair bound rejection exchange", 2143 + species);
        require_condition(
              memcmp(
                    &diagnostics[species], &expected_hard_diagnostics,
                    sizeof(diagnostics[species]))
                    == 0,
              "pair bound rejection source diagnostics were not reset", 2145 + species);
        ghl_m1_neutrino_diagnostics expected_neutrino_diagnostics
              = neutrino_diagnostics_before[species];
        expected_neutrino_diagnostics.source_failures++;
        require_condition(
              memcmp(
                    &neutrino_diagnostics[species], &expected_neutrino_diagnostics,
                    sizeof(neutrino_diagnostics[species]))
                    == 0,
              "pair bound rejection diagnostics changed unexpectedly", 2147 + species);
      }
    }
  }
}

/* The pair stage derives each species' current from the state it is about to
 * advance (start of every substep) and from the published endpoint. J_floor is
 * checked only by that derivation, so absorption that drives J below it must
 * reject the update transactionally at either site. */
static void
test_pair_state_derivation_failures(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric = { 0 };
  metric.lapse = metric.lapseinv = metric.lapseinv2 = 1.0;
  metric.detgamma = metric.sqrt_detgamma = 1.0;
  for(int direction = 0; direction < 3; ++direction) {
    metric.gammaDD[direction][direction] = 1.0;
    metric.gammaUU[direction][direction] = 1.0;
  }
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  const ghl_m1_neutrino_state state[2]
        = { { .N = 1.0, .E = 1.0 }, { .N = 1.0, .E = 1.0 } };

  /* Case 0: no number source, so one substep succeeds and only the endpoint
   * derivation sees E = 0.13 < J_floor. Case 1: the one-substep number update
   * (0.0116) falls below N_floor and is retried; after the first of two
   * substeps (N = 0.0129, E = 0.52) the next substep's start derivation sees
   * J < J_floor. */
  const double pair_number[2] = { 0.0, 1.0 };
  const double N_floor[2] = { 1.0e-12, 0.012 };
  const double J_floor[2] = { 0.5, 0.7 };
  for(int test_case = 0; test_case < 2; ++test_case) {
    ghl_m1_neutrino_parameters nu_params[2];
    ghl_m1_neutrino_rates rates[2];
    for(int species = 0; species < 2; ++species) {
      make_neutrino_parameters(&nu_params[species]);
      nu_params[species].N_floor = N_floor[test_case];
      nu_params[species].J_floor = J_floor[test_case];
      make_rates(
            species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0,
            0.01, 2.0, pair_number[test_case], 0.01, &rates[species]);
    }
    ghl_m1_neutrino_state output[2];
    ghl_m1_neutrino_exchange exchange[2];
    ghl_m1_neutrino_source_diagnostics diagnostics[2];
    ghl_m1_neutrino_diagnostics neutrino_diagnostics[2] = { { 0 }, { 0 } };
    require_error(
          ghl_m1_solve_neutrino_pair_source_update(
                m1_params, nu_params, &metric, &prims, rates, state, state, 1.0, 1.0,
                output, exchange, diagnostics, neutrino_diagnostics),
          ghl_error_m1_invalid_state, "pair current derivation failure",
          2150 + test_case);
    for(int species = 0; species < 2; ++species) {
      require_condition(
            memcmp(&output[species], &state[species], sizeof(output[species])) == 0,
            "pair derivation failure changed state", 2152 + test_case);
      require_zero_exchange(
            &exchange[species], "pair derivation failure exchange", 2154 + test_case);
      require_condition(
            diagnostics[species].path == ghl_m1_neutrino_source_path_hard_failure
                  && !diagnostics[species].terminal_no_update
                  && neutrino_diagnostics[species].source_failures == 1
                  && neutrino_diagnostics[species].source_converged == 0,
            "pair derivation failure diagnostics", 2156 + test_case);
    }
  }
}

/* The stage diagnostic is false when only the pair E/F solve leaves the
 * converged closure branch, so the published flag must come from the pair
 * solve's own solution-path flags. */
static void
test_pair_newton_only_closure_fallback(const ghl_m1_parameters *restrict m1_params) {
  ghl_m1_parameters parameters = *m1_params;
  parameters.closure_root_max_iterations = 8;
  parameters.closure_root_residual_tolerance = 1.0;
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 2.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params[2];
  ghl_m1_neutrino_rates active_rates[2], inactive_rates[2];
  for(int species = 0; species < 2; ++species) {
    const ghl_m1_neutrino_species_t kind
          = species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue;
    make_neutrino_parameters(&nu_params[species]);
    make_rates(kind, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.1, &active_rates[species]);
    make_rates(kind, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, &inactive_rates[species]);
  }
  const ghl_m1_neutrino_state input[2]
        = { { .N = 1.0, .E = 1.0, .F = { -0.8, 0.0, 0.0 } },
            { .N = 1.0, .E = 1.0, .F = { -0.8, 0.0, 0.0 } } };
  ghl_m1_neutrino_state output[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics neutrino_diagnostics[2] = { { 0 }, { 0 } };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              &parameters, nu_params, &metric, &prims, inactive_rates, input, input, 1.0,
              1.0, output, exchange, diagnostics, neutrino_diagnostics),
        ghl_success, "pair-inactive closure precheck", 2160);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          !diagnostics[species].closure_fallback_used,
          "closure fallback fixture failed its independent-stage precheck", 2161);
    neutrino_diagnostics[species] = (ghl_m1_neutrino_diagnostics){ 0 };
  }
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              &parameters, nu_params, &metric, &prims, active_rates, input, input, 1.0,
              1.0, output, exchange, diagnostics, neutrino_diagnostics),
        ghl_success, "pair-only closure fallback", 2162);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          diagnostics[species].closure_fallback_used
                && (diagnostics[species].implicit.solution_path_flags
                    & ghl_m1_solution_path_closure_fallback)
                         != 0,
          "pair-only closure fallback was not reported", 2163 + species);
  }
}

static void test_pair_public_transaction_publication(
      const ghl_m1_parameters *restrict initialized_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params[2];
  make_neutrino_parameters(&nu_params[0]);
  make_neutrino_parameters(&nu_params[1]);

  ghl_m1_neutrino_rates failure_rates[2];
  make_rates(
        ghl_m1_neutrino_nue, 0.08, 0.12, 0.10, 1.0, 2.0, 0.02, 0.03, &failure_rates[0]);
  make_rates(
        ghl_m1_neutrino_anue, 0.08, 0.12, 0.10, 1.0, 2.1, 0.02, 0.04, &failure_rates[1]);
  ghl_m1_neutrino_state state_transport[2];
  ghl_m1_neutrino_state state_input[2];
  for(int species = 0; species < 2; ++species) {
    make_state(rng, &metric, &state_transport[species]);
    state_input[species] = state_transport[species];
    state_input[species].N *= 1.1;
  }

  /* An invalid timestep is rejected before either independent source stage.
   * The public pair API must still publish its complete hard-failure packet. */
  ghl_m1_neutrino_state failure_out[2]
        = { { .N = -1.0, .E = -2.0, .F = { -3.0, -4.0, -5.0 } },
            { .N = -6.0, .E = -7.0, .F = { -8.0, -9.0, -10.0 } } };
  ghl_m1_neutrino_exchange failure_exchange[2]
        = { { .dN_rad_total = -1.0,
              .dL_rad_cc = -2.0,
              .dE_rad = -3.0,
              .dF_rad = { -4.0, -5.0, -6.0 },
              .dTau_matter = -7.0,
              .dS_matter = { -8.0, -9.0, -10.0 },
              .dYe_matter = -11.0 },
            { .dN_rad_total = -12.0,
              .dL_rad_cc = -13.0,
              .dE_rad = -14.0,
              .dF_rad = { -15.0, -16.0, -17.0 },
              .dTau_matter = -18.0,
              .dS_matter = { -19.0, -20.0, -21.0 },
              .dYe_matter = -22.0 } };
  ghl_m1_neutrino_source_diagnostics failure_diagnostics[2] = { { 0 }, { 0 } };
  ghl_m1_neutrino_diagnostics failure_nd[2];
  ghl_m1_neutrino_diagnostics failure_nd_before[2];
  for(int species = 0; species < 2; ++species) {
    failure_diagnostics[species].path = ghl_m1_neutrino_source_path_thin_explicit;
    failure_diagnostics[species].closure_fallback_used = true;
    failure_diagnostics[species].terminal_no_update = true;
    failure_diagnostics[species].implicit = (ghl_m1_implicit_solve_diagnostics){
      .newton_iterations = 2 + species,
      .line_search_backtracks = 3 + species,
      .fallback_substeps = 4 + species,
      .used_fallback_substepping = true,
      .residual_max_norm = 5.0 + species,
      .residual_scaled_norm = 6.0 + species,
      .solution_path_flags = ghl_m1_solution_path_closure_fallback
    };
    ghl_m1_neutrino_diagnostics_initialize(&failure_nd[species]);
    failure_nd[species].provider_validation_failures = 2 + species;
    failure_nd[species].source_converged = 7 + species;
    failure_nd[species].source_terminal_fallbacks = 11 + species;
    failure_nd[species].source_failures = 13 + species;
    failure_nd[species].N_floor_repairs = 17 + species;
    failure_nd[species].EF_repairs = 19 + species;
    failure_nd[species].limiter_reductions = 23 + species;
    failure_nd[species].mean_energy_diag = 0.25 + species;
    failure_nd[species].mean_energy_consistent = 29 + species;
    failure_nd[species].Jeq_over_neq_consistent = 31 + species;
    failure_nd[species].mean_energy_diag_invalid = 37 + species;
    failure_nd[species].repair_dN = 0.25 + species;
    failure_nd[species].repair_dE = 0.5 + species;
    failure_nd[species].repair_dF[0] = 0.75 + species;
    failure_nd[species].repair_dF[1] = 1.0 + species;
    failure_nd[species].repair_dF[2] = 1.25 + species;
    failure_nd[species].repair_dL_e = 1.5 + species;
    failure_nd[species].repair_stage = species == 0
                                             ? ghl_m1_neutrino_repair_stage_post_source
                                             : ghl_m1_neutrino_repair_stage_face_state;
    failure_nd[species].repair_lepton_weight = 0.5 + species;
    failure_nd[species].rate_product_underflows = 41 + species;
    failure_nd_before[species] = failure_nd[species];
  }
  ghl_m1_neutrino_source_diagnostics expected_failure_source_diagnostics = { 0 };
  ghl_m1_initialize_implicit_solve_diagnostics(
        &expected_failure_source_diagnostics.implicit);
  expected_failure_source_diagnostics.path = ghl_m1_neutrino_source_path_hard_failure;
  ghl_error_codes_t error = ghl_m1_solve_neutrino_pair_source_update(
        initialized_params, nu_params, &metric, &prims, failure_rates, state_input,
        state_transport, -0.25, 3.0, failure_out, failure_exchange, failure_diagnostics,
        failure_nd);
  require_error(
        error, ghl_error_m1_invalid_state, "pair invalid-timestep publication", 3300);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(
                &failure_out[species], &state_transport[species],
                sizeof(failure_out[species]))
                == 0,
          "pair hard failure changed transport state", 3300);
    require_zero_exchange(
          &failure_exchange[species], "pair hard-failure exchange", 3300);
    require_condition(
          memcmp(
                &failure_diagnostics[species], &expected_failure_source_diagnostics,
                sizeof(failure_diagnostics[species]))
                == 0,
          "pair hard-failure source diagnostics were not reset", 3300);
    ghl_m1_neutrino_diagnostics expected_failure_neutrino_diagnostics
          = failure_nd_before[species];
    expected_failure_neutrino_diagnostics.source_failures++;
    require_condition(
          memcmp(
                &failure_nd[species], &expected_failure_neutrino_diagnostics,
                sizeof(failure_nd[species]))
                == 0,
          "pair hard failure diagnostics were not isolated", 3300);
  }

  /* Keep the independent stage trivial and valid, then make only the active
   * pair E/F solve exhaust every retry schedule. This reaches the public
   * terminal publication path rather than the earlier hard-failure path. */
  ghl_m1_neutrino_rates terminal_rates[2];
  make_rates(
        ghl_m1_neutrino_nue, 0.0, 0.0, 0.0, 1.0, 2.0, 0.02, 0.03, &terminal_rates[0]);
  make_rates(
        ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 1.0, 2.1, 0.02, 0.03, &terminal_rates[1]);
  ghl_m1_parameters terminal_params = *initialized_params;
  terminal_params.newton_max_iterations = 1;
  terminal_params.newton_tolerance = DBL_MIN;
  terminal_params.newton_absolute_tolerance = DBL_MIN;
  ghl_m1_neutrino_state terminal_out[2]
        = { { .N = -23.0, .E = -24.0, .F = { -25.0, -26.0, -27.0 } },
            { .N = -28.0, .E = -29.0, .F = { -30.0, -31.0, -32.0 } } };
  ghl_m1_neutrino_exchange terminal_exchange[2]
        = { { .dN_rad_total = -23.0,
              .dL_rad_cc = -24.0,
              .dE_rad = -25.0,
              .dF_rad = { -26.0, -27.0, -28.0 },
              .dTau_matter = -29.0,
              .dS_matter = { -30.0, -31.0, -32.0 },
              .dYe_matter = -33.0 },
            { .dN_rad_total = -34.0,
              .dL_rad_cc = -35.0,
              .dE_rad = -36.0,
              .dF_rad = { -37.0, -38.0, -39.0 },
              .dTau_matter = -40.0,
              .dS_matter = { -41.0, -42.0, -43.0 },
              .dYe_matter = -44.0 } };
  ghl_m1_neutrino_source_diagnostics terminal_diagnostics[2] = { { 0 }, { 0 } };
  ghl_m1_neutrino_diagnostics terminal_nd[2];
  ghl_m1_neutrino_diagnostics terminal_nd_before[2];
  for(int species = 0; species < 2; ++species) {
    terminal_diagnostics[species].path = ghl_m1_neutrino_source_path_thick_equilibrium;
    terminal_diagnostics[species].closure_fallback_used = true;
    terminal_diagnostics[species].terminal_no_update = false;
    terminal_diagnostics[species].implicit = (ghl_m1_implicit_solve_diagnostics){
      .newton_iterations = 8 + species,
      .line_search_backtracks = 9 + species,
      .fallback_substeps = 10 + species,
      .used_fallback_substepping = true,
      .residual_max_norm = 11.0 + species,
      .residual_scaled_norm = 12.0 + species,
      .solution_path_flags = ghl_m1_solution_path_terminal_failure
    };
    ghl_m1_neutrino_diagnostics_initialize(&terminal_nd[species]);
    terminal_nd[species].provider_validation_failures = 3 + species;
    terminal_nd[species].source_converged = 11 + species;
    terminal_nd[species].source_terminal_fallbacks = 17 + species;
    terminal_nd[species].source_failures = 19 + species;
    terminal_nd[species].N_floor_repairs = 23 + species;
    terminal_nd[species].EF_repairs = 29 + species;
    terminal_nd[species].limiter_reductions = 31 + species;
    terminal_nd[species].mean_energy_diag = 0.75 + species;
    terminal_nd[species].mean_energy_consistent = 37 + species;
    terminal_nd[species].Jeq_over_neq_consistent = 41 + species;
    terminal_nd[species].mean_energy_diag_invalid = 43 + species;
    terminal_nd[species].repair_dE = 0.5 + species;
    terminal_nd[species].repair_dN = 0.75 + species;
    terminal_nd[species].repair_dF[0] = 1.0 + species;
    terminal_nd[species].repair_dF[1] = 1.25 + species;
    terminal_nd[species].repair_dF[2] = 1.5 + species;
    terminal_nd[species].repair_dL_e = 1.75 + species;
    terminal_nd[species].repair_stage
          = species == 0 ? ghl_m1_neutrino_repair_stage_post_transport
                         : ghl_m1_neutrino_repair_stage_post_source;
    terminal_nd[species].repair_lepton_weight = 0.75 + species;
    terminal_nd[species].rate_product_underflows = 47 + species;
    terminal_nd_before[species] = terminal_nd[species];
  }
  ghl_m1_neutrino_source_diagnostics expected_terminal_source_diagnostics = { 0 };
  ghl_m1_initialize_implicit_solve_diagnostics(
        &expected_terminal_source_diagnostics.implicit);
  expected_terminal_source_diagnostics.path
        = ghl_m1_neutrino_source_path_terminal_no_update;
  expected_terminal_source_diagnostics.terminal_no_update = true;
  error = ghl_m1_solve_neutrino_pair_source_update(
        &terminal_params, nu_params, &metric, &prims, terminal_rates, state_input,
        state_transport, 0.5, 3.0, terminal_out, terminal_exchange, terminal_diagnostics,
        terminal_nd);
  require_error(
        error, ghl_error_m1_implicit_terminal_fallback, "pair terminal publication",
        3301);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(
                &terminal_out[species], &state_transport[species],
                sizeof(terminal_out[species]))
                == 0,
          "pair terminal changed transport state", 3301);
    require_zero_exchange(&terminal_exchange[species], "pair terminal exchange", 3301);
    require_condition(
          memcmp(
                &terminal_diagnostics[species], &expected_terminal_source_diagnostics,
                sizeof(terminal_diagnostics[species]))
                == 0,
          "pair terminal source diagnostics were not reset", 3301);
    ghl_m1_neutrino_diagnostics expected_terminal_neutrino_diagnostics
          = terminal_nd_before[species];
    expected_terminal_neutrino_diagnostics.source_terminal_fallbacks++;
    require_condition(
          memcmp(
                &terminal_nd[species], &expected_terminal_neutrino_diagnostics,
                sizeof(terminal_nd[species]))
                == 0,
          "pair terminal diagnostics were not isolated", 3301);
  }
}

static void test_repairs_and_transactional_failures(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_parameters repair_nu_params = nu_params;
  repair_nu_params.N_floor = 1.0;
  ghl_m1_neutrino_state repair_state = { .N = 0.25,
                                         .E = 0.5 * m1_params->E_floor,
                                         .F = { 2.0 * m1_params->E_floor, 0.0, 0.0 } };
  ghl_m1_neutrino_diagnostics repair_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&repair_diagnostics);
  repair_diagnostics.repair_stage = ghl_m1_neutrino_repair_stage_post_source;
  repair_diagnostics.repair_lepton_weight = 1.0;
  const ghl_error_codes_t repair_error = ghl_m1_repair_neutrino_state(
        m1_params, &repair_nu_params, &metric, &repair_state, &repair_diagnostics);
  require_error(repair_error, ghl_success, "neutrino state repair", 3000);
  require_close(
        repair_state.N, repair_nu_params.N_floor, 0.0, 0.0, "number floor repair", 3000);
  require_close(
        repair_state.E, m1_params->E_floor, 0.0, 0.0, "energy floor repair", 3000);
  require_condition(
        repair_diagnostics.N_floor_repairs == 1 && repair_diagnostics.EF_repairs == 1,
        "repair diagnostics did not record both repairs", 3000);
  require_condition(
        repair_diagnostics.repair_dN > 0.0 && repair_diagnostics.repair_dE > 0.0,
        "repair magnitude budget was not updated", 3000);

  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 0.10, 0.10, 0.10, 1.0, 2.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_state state_transport;
  make_state(rng, &metric, &state_transport);
  ghl_m1_neutrino_state state_input = state_transport;
  ghl_m1_neutrino_state state_out = { .N = -9.0, .E = -8.0, .F = { -7.0, -6.0, -5.0 } };
  ghl_m1_neutrino_exchange exchange = { .dN_rad_total = -1.0,
                                        .dL_rad_cc = -2.0,
                                        .dE_rad = -3.0,
                                        .dF_rad = { -4.0, -5.0, -6.0 },
                                        .dTau_matter = -7.0,
                                        .dS_matter = { -8.0, -9.0, -10.0 },
                                        .dYe_matter = -11.0 };
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics neutrino_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
  rates.eta_N = -1.0;
  const ghl_error_codes_t rate_error = ghl_m1_solve_neutrino_source_update(
        NULL, m1_params, &nu_params, &metric, &prims, &rates, &state_input,
        &state_transport, 0.1, 3.0, &state_out, &exchange, &diagnostics,
        &neutrino_diagnostics);
  require_error(
        rate_error, ghl_error_m1_microphysics_failure, "invalid-rate source update",
        3001);
  require_condition(
        memcmp(&state_out, &state_transport, sizeof(state_out)) == 0,
        "invalid-rate failure changed output state", 3001);
  require_zero_exchange(&exchange, "invalid-rate transactional exchange", 3001);
  require_condition(
        diagnostics.path == ghl_m1_neutrino_source_path_hard_failure,
        "invalid-rate failure path was not reported", 3001);

  make_rates(ghl_m1_neutrino_nue, 0.10, 0.10, 0.10, 1.0, 2.0, 0.0, 0.0, &rates);
  state_input.E = NAN;
  state_out = (ghl_m1_neutrino_state){ .N = -1.0, .E = -2.0, .F = { -3.0, -4.0, -5.0 } };
  exchange.dE_rad = -13.0;
  const ghl_error_codes_t state_error = ghl_m1_solve_neutrino_source_update(
        NULL, m1_params, &nu_params, &metric, &prims, &rates, &state_input,
        &state_transport, 0.1, 3.0, &state_out, &exchange, &diagnostics,
        &neutrino_diagnostics);
  require_error(
        state_error, ghl_error_m1_invalid_state, "invalid-state source update", 3002);
  require_condition(
        memcmp(&state_out, &state_transport, sizeof(state_out)) == 0,
        "invalid-state failure changed output state", 3002);
  require_zero_exchange(&exchange, "invalid-state transactional exchange", 3002);

  state_input = state_transport;
  ghl_m1_neutrino_source_options double_application = branched_options();
  double_application.interaction_sources_already_applied = true;
  state_out = (ghl_m1_neutrino_state){ .N = -1.0, .E = -2.0, .F = { -3.0, -4.0, -5.0 } };
  const ghl_error_codes_t double_error = ghl_m1_solve_neutrino_source_update(
        &double_application, m1_params, &nu_params, &metric, &prims, &rates,
        &state_input, &state_transport, 0.1, 3.0, &state_out, &exchange, &diagnostics,
        &neutrino_diagnostics);
  require_error(
        double_error, ghl_error_m1_source_double_application,
        "double-application source update", 3003);
  require_condition(
        memcmp(&state_out, &state_transport, sizeof(state_out)) == 0,
        "double-application failure changed output state", 3003);
  require_zero_exchange(&exchange, "double-application exchange", 3003);
}

static void
test_stiff_charged_current_conservation(const ghl_m1_parameters *restrict m1_params) {
  m1_test_rng rng = { .state = UINT64_C(0x43435f5354494646) };
  ghl_metric_quantities metric;
  make_metric(&rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  const ghl_m1_neutrino_source_options branched = branched_options();
  const ghl_m1_neutrino_state input = { .N = 0.5, .E = 1.0, .F = { 0.0, 0.0, 0.0 } };
  for(int species = ghl_m1_neutrino_nue; species <= ghl_m1_neutrino_anue; ++species) {
    for(int route = 0; route < 4; ++route) {
      ghl_m1_neutrino_rates rates;
      make_rates(
            species, 1.0e16, route == 2 ? 1.0 : 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, &rates);
      ghl_m1_neutrino_state output;
      ghl_m1_neutrino_exchange exchange;
      ghl_m1_neutrino_diagnostics nd;
      ghl_m1_neutrino_diagnostics_initialize(&nd);
      ghl_error_codes_t error;
      if(route == 3) {
        ghl_m1_implicit_solve_diagnostics diagnostics;
        error = ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &metric, &prims, &rates, 1.0, 2.0, &input, &output,
              &exchange, &diagnostics, &nd);
      }
      else {
        ghl_m1_neutrino_source_diagnostics diagnostics;
        error = ghl_m1_solve_neutrino_source_update(
              route == 0 ? NULL : &branched, m1_params, &nu_params, &metric, &prims,
              &rates, &input, &input, 1.0, 2.0, &output, &exchange, &diagnostics, &nd);
        const ghl_m1_neutrino_source_path_t expected[]
              = { ghl_m1_neutrino_source_path_general_implicit,
                  ghl_m1_neutrino_source_path_thin_explicit,
                  ghl_m1_neutrino_source_path_thick_equilibrium };
        require_condition(
              diagnostics.path == expected[route], "stiff charged-current route", route);
      }
      require_error(error, ghl_success, "stiff charged-current update", route);
      require_close(
            output.N, 1.0, 1.0e-12, 1.0e-13, "stiff charged-current number endpoint",
            route);
      require_close(
            exchange.dL_rad_cc, rates.lepton_weight * (output.N - input.N), 1.0e-12,
            1.0e-13, "stiff charged-current conservation", route);
      require_close(
            2.0 * exchange.dYe_matter + exchange.dL_rad_cc, 0.0, 0.0, 1.0e-13,
            "stiff matter-lepton compensation", route);
    }
  }
}

static void test_number_floor_charged_current_accounting(
      const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  m1_test_rng rng = { .state = UINT64_C(0x4d315f464c4f4f52) };
  make_metric(&rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);

  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  nu_params.N_floor = 0.9;

  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 1.0, 0.0, 0.0, 0.1, 1.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_rates thick_rates;
  make_rates(ghl_m1_neutrino_nue, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, &thick_rates);

  const ghl_m1_neutrino_state state_transport
        = { .N = 1.0, .E = 1.0, .F = { 0.0, 0.0, 0.0 } };
  const ghl_m1_neutrino_state state_input = state_transport;
  const double expected_number_endpoint
        = (state_transport.N + rates.eta_N) / (1.0 + rates.kappa_a_N);
  const double expected_dL_rad_cc
        = rates.lepton_weight
          * (rates.eta_N_cc - rates.kappa_a_N_cc * expected_number_endpoint);
  const double expected_repair = nu_params.N_floor - expected_number_endpoint;
  const ghl_m1_neutrino_source_options branched = branched_options();
  const ghl_m1_neutrino_source_options *const option_cases[3]
        = { NULL, &branched, &branched };
  const ghl_m1_neutrino_rates *const rate_cases[3] = { &rates, &rates, &thick_rates };
  const double expected_endpoints[3]
        = { expected_number_endpoint, expected_number_endpoint, 0.5 };
  const double expected_dL_rad_cc_cases[3]
        = { expected_dL_rad_cc, expected_dL_rad_cc, -0.5 };
  const double expected_repairs[3] = { expected_repair, expected_repair,
                                       nu_params.N_floor - expected_endpoints[2] };
  const ghl_m1_neutrino_source_path_t expected_paths[3]
        = { ghl_m1_neutrino_source_path_general_implicit,
            ghl_m1_neutrino_source_path_thin_explicit,
            ghl_m1_neutrino_source_path_thick_equilibrium };

  for(int case_index = 0; case_index < 3; ++case_index) {
    ghl_m1_neutrino_state state_out;
    ghl_m1_neutrino_exchange exchange;
    ghl_m1_neutrino_source_diagnostics diagnostics;
    ghl_m1_neutrino_diagnostics neutrino_diagnostics;
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics);
    neutrino_diagnostics.repair_lepton_weight = rate_cases[case_index]->lepton_weight;

    const ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
          option_cases[case_index], m1_params, &nu_params, &metric, &prims,
          rate_cases[case_index], &state_input, &state_transport, 1.0, 1.0, &state_out,
          &exchange, &diagnostics, &neutrino_diagnostics);
    require_error(error, ghl_success, "number-floor source update", 3050 + case_index);
    require_condition(
          diagnostics.path == expected_paths[case_index], "number-floor source path",
          3050 + case_index);
    require_close(
          state_out.N, nu_params.N_floor, 1.0e-12, 1.0e-13, "number-floor endpoint",
          3050 + case_index);
    require_close(
          exchange.dL_rad_cc, expected_dL_rad_cc_cases[case_index], 1.0e-12, 1.0e-13,
          "charged-current exchange before repair", 3050 + case_index);
    require_close(
          exchange.dYe_matter, -expected_dL_rad_cc_cases[case_index], 1.0e-12, 1.0e-13,
          "matter lepton exchange", 3050 + case_index);
    require_close(
          neutrino_diagnostics.repair_dN, expected_repairs[case_index], 1.0e-12, 1.0e-13,
          "number-floor repair budget", 3050 + case_index);
    require_close(
          neutrino_diagnostics.repair_dL_e, expected_repairs[case_index], 1.0e-12,
          1.0e-13, "lepton repair budget", 3050 + case_index);
    require_close(
          exchange.dN_rad_total + exchange.dYe_matter, expected_repairs[case_index],
          1.0e-12, 1.0e-13, "radiation-matter lepton accounting", 3050 + case_index);
  }
}

static void test_lepton_policies(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_state state;
  make_state(rng, &metric, &state);
  ghl_m1_neutrino_rates rates;
  ghl_m1_neutrino_state output;
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics nd;

  for(int species = ghl_m1_neutrino_nue; species <= ghl_m1_neutrino_anue; ++species) {
    make_rates(
          (ghl_m1_neutrino_species_t)species, 0.12, 0.20, 0.10, 1.0, 2.0, 0.0, 0.0,
          &rates);
    ghl_m1_neutrino_diagnostics_initialize(&nd);
    ghl_m1_neutrino_source_options options = branched_options();
    options.thick_equilibrium_threshold = 0.0;
    options.scattering_threshold = 0.0;
    ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
          &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state, 0.1,
          4.0, &output, &exchange, &diagnostics, &nd);
    require_error(error, ghl_success, "lepton-policy source update", 3100 + species);
    require_exchange_contract(&state, &output, &exchange, &metric, 4.0, 3100 + species);
    require_condition(
          (species == ghl_m1_neutrino_nue && rates.lepton_weight == 1.0)
                || (species == ghl_m1_neutrino_anue && rates.lepton_weight == -1.0),
          "wrong electron-flavor lepton weight", 3100 + species);

    options.ye_policy = ghl_m1_neutrino_ye_from_signed_total_number;
    ghl_m1_neutrino_diagnostics_initialize(&nd);
    ghl_m1_neutrino_exchange signed_exchange;
    error = ghl_m1_solve_neutrino_source_update(
          &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state, 0.1,
          4.0, &output, &signed_exchange, &diagnostics, &nd);
    require_error(error, ghl_success, "signed-total lepton policy", 3110 + species);
    require_close(
          signed_exchange.dYe_matter, exchange.dYe_matter, 2.0e-9, 2.0e-12,
          "charged-current versus signed-total Ye", 3110 + species);
  }

  make_rates(ghl_m1_neutrino_nux, 0.12, 0.20, 0.10, 1.0, 2.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_diagnostics_initialize(&nd);
  ghl_m1_neutrino_source_options options = branched_options();
  options.ye_policy = ghl_m1_neutrino_ye_from_signed_total_number;
  options.thick_equilibrium_threshold = 0.0;
  options.scattering_threshold = 0.0;
  const ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
        &options, m1_params, &nu_params, &metric, &prims, &rates, &state, &state, 0.1,
        4.0, &output, &exchange, &diagnostics, &nd);
  require_error(error, ghl_success, "heavy-flavor lepton policy", 3120);
  require_close(
        exchange.dL_rad_cc, 0.0, 0.0, 0.0, "heavy-flavor charged-current exchange",
        3120);
  require_close(
        exchange.dYe_matter, 0.0, 0.0, 0.0, "heavy-flavor electron-fraction exchange",
        3120);
}

static void test_retry_and_terminal_no_update(
      const ghl_m1_parameters *restrict initialized_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 40.0, 80.0, 40.0, 1.0, 2.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_state state;
  make_state(rng, &metric, &state);
  ghl_m1_neutrino_state output;
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics nd;

  /* A deliberately short Newton budget exercises the bounded implicit solve
   * while retaining a valid successful endpoint. The terminal case below
   * covers the exhausted retry schedule. */
  ghl_m1_parameters retry_params = *initialized_params;
  retry_params.newton_max_iterations = 2;
  retry_params.newton_tolerance = 1.0e-12;
  retry_params.newton_absolute_tolerance = 1.0e-14;
  ghl_m1_neutrino_diagnostics_initialize(&nd);
  ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
        NULL, &retry_params, &nu_params, &metric, &prims, &rates, &state, &state, 0.5,
        3.0, &output, &exchange, &diagnostics, &nd);
  require_error(error, ghl_success, "retrying implicit source update", 3200);
  require_condition(
        diagnostics.implicit.fallback_substeps >= 1
              && diagnostics.implicit.newton_iterations > 0
              && isfinite(diagnostics.implicit.residual_scaled_norm),
        "bounded implicit solve diagnostics are invalid", 3200);
  require_state_finite_and_admissible(
        initialized_params, &nu_params, &metric, &output, 3200);

  /* Make convergence effectively impossible while retaining a valid M1
   * configuration. Every schedule must then return the contract's terminal
   * no-update status and preserve the transport state with zero exchange. */
  ghl_m1_parameters terminal_params = *initialized_params;
  terminal_params.newton_max_iterations = 1;
  terminal_params.newton_tolerance = DBL_MIN;
  terminal_params.newton_absolute_tolerance = DBL_MIN;
  ghl_m1_neutrino_diagnostics_initialize(&nd);
  output = (ghl_m1_neutrino_state){ .N = -1.0, .E = -2.0, .F = { -3.0, -4.0, -5.0 } };
  exchange.dE_rad = -3.0;
  error = ghl_m1_solve_neutrino_source_update(
        NULL, &terminal_params, &nu_params, &metric, &prims, &rates, &state, &state, 0.5,
        3.0, &output, &exchange, &diagnostics, &nd);
  require_error(
        error, ghl_error_m1_implicit_terminal_fallback,
        "terminal implicit source update", 3201);
  require_condition(
        diagnostics.path == ghl_m1_neutrino_source_path_terminal_no_update
              && diagnostics.terminal_no_update,
        "terminal source path was not reported", 3201);
  require_condition(
        memcmp(&output, &state, sizeof(output)) == 0,
        "terminal fallback changed the transport state", 3201);
  require_zero_exchange(&exchange, "terminal fallback exchange", 3201);
  require_condition(
        nd.source_terminal_fallbacks > 0, "terminal fallback was not counted", 3201);

  /* The paired operation must propagate an independent-stage terminal
   * failure atomically across both species, before its pair stage runs. */
  ghl_m1_neutrino_parameters pair_params[2] = { nu_params, nu_params };
  ghl_m1_neutrino_rates pair_rates[2] = { rates, rates };
  pair_rates[1].species = ghl_m1_neutrino_anue;
  pair_rates[1].lepton_weight = -1.0;
  const ghl_m1_neutrino_state pair_input[2] = { state, state };
  ghl_m1_neutrino_state pair_output[2];
  ghl_m1_neutrino_exchange pair_exchange[2];
  ghl_m1_neutrino_source_diagnostics pair_diagnostics[2];
  ghl_m1_neutrino_diagnostics pair_nd[2] = { { 0 }, { 0 } };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              &terminal_params, pair_params, &metric, &prims, pair_rates, pair_input,
              pair_input, 0.5, 3.0, pair_output, pair_exchange, pair_diagnostics,
              pair_nd),
        ghl_error_m1_implicit_terminal_fallback,
        "pair independent-stage terminal fallback", 3202);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&pair_output[species], &pair_input[species], sizeof(state)) == 0
                && pair_diagnostics[species].terminal_no_update
                && pair_diagnostics[species].path
                         == ghl_m1_neutrino_source_path_terminal_no_update
                && pair_nd[species].source_terminal_fallbacks == 1,
          "pair independent-stage terminal failure was not atomic", 3202);
    require_zero_exchange(&pair_exchange[species], "pair terminal exchange", 3202);
  }
}

static void
test_implicit_predictor_recovery(const ghl_m1_parameters *restrict initialized_params) {
  ghl_metric_quantities metric = { 0 };
  metric.lapse = metric.lapseinv = metric.lapseinv2 = 1.0;
  metric.detgamma = metric.sqrt_detgamma = 1.0;
  for(int direction = 0; direction < 3; ++direction) {
    metric.gammaDD[direction][direction] = 1.0;
    metric.gammaUU[direction][direction] = 1.0;
  }
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);

  /* With zero flux, E=4 and J_eq=2, the explicit predictor is E*=0 while
   * the implicit endpoint is 8/3. This is a supported homogeneous-solver
   * input and selects the predictor reset before Newton starts. */
  ghl_m1_neutrino_rates reset_rates;
  make_rates(ghl_m1_neutrino_nue, 0.0, 2.0, 0.0, 1.0, 2.0, 0.0, 0.0, &reset_rates);
  const ghl_m1_neutrino_state reset_input
        = { .N = 1.0, .E = 4.0, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_neutrino_state reset_output
        = { .N = -1.0, .E = -2.0, .F = { -3.0, -4.0, -5.0 } };
  ghl_m1_neutrino_exchange reset_exchange = { .dE_rad = -6.0 };
  ghl_m1_implicit_solve_diagnostics reset_solve_diagnostics;
  ghl_m1_neutrino_diagnostics reset_neutrino_diagnostics;
  ghl_m1_initialize_implicit_solve_diagnostics(&reset_solve_diagnostics);
  ghl_m1_neutrino_diagnostics_initialize(&reset_neutrino_diagnostics);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              initialized_params, &nu_params, &metric, &prims, &reset_rates, 1.0, 1.0,
              &reset_input, &reset_output, &reset_exchange, &reset_solve_diagnostics,
              &reset_neutrino_diagnostics),
        ghl_success, "implicit predictor reset recovery", 3302);
  require_close(
        reset_output.E, 8.0 / 3.0, 2.0e-9, 2.0e-12,
        "implicit predictor reset endpoint energy", 3302);
  require_close(
        reset_output.N, reset_input.N, 0.0, 0.0,
        "implicit predictor reset endpoint number", 3302);
  require_condition(
        reset_output.F[0] == 0.0 && reset_output.F[1] == 0.0 && reset_output.F[2] == 0.0
              && reset_solve_diagnostics.newton_iterations > 0
              && reset_neutrino_diagnostics.source_converged == 1,
        "implicit predictor reset did not publish a converged endpoint", 3302);

  /* This finite SPD metric is rejected by the public current-state/closure
   * validation before Newton is entered: lapse-squared and densitization
   * overflow in that validated boundary. The public transaction must publish
   * the hard failure without modifying its output state or exchange. */
  const double spatial_scale = ldexp(1.0, 500);
  ghl_metric_quantities overflow_metric = { 0 };
  overflow_metric.lapse = ldexp(1.0, 600);
  overflow_metric.lapseinv = ldexp(1.0, -600);
  overflow_metric.lapseinv2 = ldexp(1.0, -1200);
  overflow_metric.detgamma = ldexp(1.0, 1000);
  overflow_metric.sqrt_detgamma = spatial_scale;
  overflow_metric.gammaDD[0][0] = spatial_scale;
  overflow_metric.gammaDD[1][1] = spatial_scale;
  overflow_metric.gammaDD[2][2] = 1.0;
  overflow_metric.gammaUU[0][0] = ldexp(1.0, -500);
  overflow_metric.gammaUU[1][1] = ldexp(1.0, -500);
  overflow_metric.gammaUU[2][2] = 1.0;
  const ghl_m1_neutrino_state overflow_input
        = { .N = 1.0, .E = 1.0, .F = { 0.0, 0.0, 0.0 } };
  ghl_m1_neutrino_rates transparent_rates;
  make_rates(ghl_m1_neutrino_nue, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, &transparent_rates);
  ghl_m1_neutrino_state overflow_output
        = { .N = -7.0, .E = -8.0, .F = { -9.0, -10.0, -11.0 } };
  ghl_m1_neutrino_exchange overflow_exchange = { .dE_rad = -12.0 };
  ghl_m1_implicit_solve_diagnostics overflow_solve_diagnostics;
  ghl_m1_neutrino_diagnostics overflow_neutrino_diagnostics;
  ghl_m1_initialize_implicit_solve_diagnostics(&overflow_solve_diagnostics);
  ghl_m1_neutrino_diagnostics_initialize(&overflow_neutrino_diagnostics);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              initialized_params, &nu_params, &overflow_metric, &prims,
              &transparent_rates, 1.0, 1.0, &overflow_input, &overflow_output,
              &overflow_exchange, &overflow_solve_diagnostics,
              &overflow_neutrino_diagnostics),
        ghl_error_m1_invalid_state, "implicit predictor error reset", 3303);
  require_condition(
        memcmp(&overflow_output, &overflow_input, sizeof(overflow_output)) == 0
              && memcmp(
                       &overflow_exchange, &(ghl_m1_neutrino_exchange){ 0 },
                       sizeof(overflow_exchange))
                       == 0
              && overflow_neutrino_diagnostics.source_failures == 1,
        "implicit predictor error was not transactional", 3303);
}

static void require_invalid_rates(
      const ghl_m1_neutrino_rates *restrict rates,
      const int case_index) {
  ghl_m1_neutrino_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&diagnostics);
  require_error(
        ghl_m1_validate_neutrino_rates(rates, &diagnostics),
        ghl_error_m1_microphysics_failure, "direct invalid-rate validation", case_index);
  require_condition(
        diagnostics.provider_validation_failures == 1,
        "invalid rate was not counted exactly once", case_index);
}

static void test_direct_neutrino_validation_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 0.10, 0.20, 0.10, 1.0, 2.0, 0.0, 0.0, &rates);

  const ghl_m1_rad_state valid_rad_state = { .E = 1.0, .F = { 0.1, 0.0, 0.0 } };
  double flux_factor_sq = -1.0;
  require_error(
        ghl_m1_validate_realizability_state(
              m1_params, &metric, &valid_rad_state, 128.0, &flux_factor_sq),
        ghl_success, "realizability flux-factor output", 3999);
  require_close(
        flux_factor_sq,
        metric.gammaUU[0][0] * valid_rad_state.F[0] * valid_rad_state.F[0]
              / (valid_rad_state.E * valid_rad_state.E),
        1.0e-14, 1.0e-15, "realizability flux-factor square", 3999);

  /* Exercise every independently mutable field in the provider validation
   * boundary. Each case starts from the same valid Kirchhoff bundle, so the
   * assertion identifies the selected predicate rather than a later failure. */
  ghl_m1_neutrino_rates candidate = rates;
  candidate.species = (ghl_m1_neutrino_species_t)99;
  require_invalid_rates(&candidate, 4000);
  candidate = rates;
  candidate.lepton_weight = 0.0;
  require_invalid_rates(&candidate, 4001);
  candidate = rates;
  candidate.species = ghl_m1_neutrino_anue;
  candidate.lepton_weight = 0.0;
  require_invalid_rates(&candidate, 4002);
  candidate = rates;
  candidate.species = ghl_m1_neutrino_nux;
  candidate.lepton_weight = 1.0;
  require_invalid_rates(&candidate, 4003);

  for(int field = 0; field < 11; ++field) {
    candidate = rates;
    switch(field) {
      case 0:
        candidate.eta_N = NAN;
        break;
      case 1:
        candidate.eta_E = NAN;
        break;
      case 2:
        candidate.kappa_a_N = NAN;
        break;
      case 3:
        candidate.kappa_a_E = NAN;
        break;
      case 4:
        candidate.kappa_s = NAN;
        break;
      case 5:
        candidate.kappa_tr = NAN;
        break;
      case 6:
        candidate.n_eq = NAN;
        break;
      case 7:
        candidate.J_eq = NAN;
        break;
      case 8:
        candidate.mean_energy = NAN;
        break;
      case 9:
        candidate.lepton_weight = NAN;
        break;
      default:
        candidate.eta_N_cc = NAN;
        break;
    }
    require_invalid_rates(&candidate, 4010 + field);
  }

  candidate = rates;
  candidate.eta_N = -DBL_MIN;
  require_invalid_rates(&candidate, 4020);
  candidate = rates;
  candidate.eta_E = -DBL_MIN;
  require_invalid_rates(&candidate, 4021);
  candidate = rates;
  candidate.n_eq = -DBL_MIN;
  require_invalid_rates(&candidate, 4022);
  candidate = rates;
  candidate.J_eq = -DBL_MIN;
  require_invalid_rates(&candidate, 4023);

  candidate = rates;
  candidate.kappa_a_N = -DBL_MIN;
  require_invalid_rates(&candidate, 4030);
  candidate = rates;
  candidate.kappa_a_E = -DBL_MIN;
  require_invalid_rates(&candidate, 4031);
  candidate = rates;
  candidate.kappa_s = -DBL_MIN;
  require_invalid_rates(&candidate, 4032);
  candidate = rates;
  candidate.kappa_tr = -DBL_MIN;
  require_invalid_rates(&candidate, 4033);
  candidate = rates;
  candidate.eta_N_cc = -DBL_MIN;
  require_invalid_rates(&candidate, 4034);
  candidate = rates;
  candidate.kappa_a_N_cc = -DBL_MIN;
  require_invalid_rates(&candidate, 4035);
  candidate = rates;
  candidate.mean_energy = 0.0;
  require_invalid_rates(&candidate, 4036);

  candidate = rates;
  candidate.eta_N_pair[0] = NAN;
  require_invalid_rates(&candidate, 4040);
  candidate = rates;
  candidate.eta_E_pair[0] = -DBL_MIN;
  require_invalid_rates(&candidate, 4041);
  candidate = rates;
  candidate.eta_N_pair[0] = 0.1;
  candidate.n_eq = 0.0;
  candidate.J_eq = 0.0;
  candidate.eta_N = 0.0;
  candidate.eta_E = 0.0;
  candidate.kappa_a_N = 0.0;
  candidate.kappa_a_E = 0.0;
  candidate.kappa_s = 0.0;
  candidate.kappa_tr = 0.0;
  candidate.eta_N_cc = 0.0;
  candidate.kappa_a_N_cc = 0.0;
  require_invalid_rates(&candidate, 4042);
  candidate = rates;
  candidate.species = ghl_m1_neutrino_nux;
  candidate.lepton_weight = 0.0;
  candidate.eta_N_cc = 0.0;
  candidate.kappa_a_N_cc = 0.0;
  require_error(
        ghl_m1_validate_neutrino_rates(&candidate, NULL), ghl_success,
        "valid heavy-flavor rate control", 4043);
  candidate.eta_N_pair[0] = 0.1;
  candidate.eta_E_pair[0] = 0.1;
  require_invalid_rates(&candidate, 4043);

  candidate = rates;
  candidate.eta_N_cc = rates.eta_N + 0.01;
  require_invalid_rates(&candidate, 4050);
  candidate = rates;
  candidate.kappa_a_N_cc = rates.kappa_a_N + 0.01;
  require_invalid_rates(&candidate, 4051);
  candidate = rates;
  candidate.kappa_a_E = DBL_MAX;
  candidate.kappa_s = DBL_MAX;
  candidate.kappa_tr = DBL_MAX;
  candidate.eta_E = DBL_MAX;
  require_invalid_rates(&candidate, 4052);

  /* The rounding-mode guard is an actual process configuration, not a
   * fabricated object state. Restore the caller's mode immediately. */
  const int saved_rounding_mode = fegetround();
  require_condition(
        saved_rounding_mode != -1, "could not read floating-point rounding mode", 4053);
  require_condition(
        fesetround(FE_DOWNWARD) == 0, "could not select downward rounding mode", 4053);
  require_invalid_rates(&rates, 4053);
  require_condition(
        fesetround(saved_rounding_mode) == 0,
        "could not restore floating-point rounding mode", 4053);

  candidate = rates;
  candidate.kappa_a_N = DBL_MIN;
  candidate.n_eq = 0x1p-53;
  candidate.J_eq = 0x1p-52;
  candidate.eta_N = 0.0;
  candidate.eta_E = candidate.kappa_a_E * candidate.J_eq;
  candidate.kappa_a_N_cc = DBL_MIN;
  candidate.eta_N_cc = 0.0;
  ghl_m1_neutrino_diagnostics underflow_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&underflow_diagnostics);
  require_error(
        ghl_m1_validate_neutrino_rates(&candidate, &underflow_diagnostics), ghl_success,
        "accepted representational rate underflow", 4054);
  require_condition(
        underflow_diagnostics.rate_product_underflows >= 2,
        "accepted underflow products were not diagnosed", 4054);

  candidate = rates;
  candidate.kappa_a_E = 0.0;
  candidate.kappa_s = 0.0;
  candidate.kappa_tr = 0.0;
  candidate.eta_E = 0.0;
  require_error(
        ghl_m1_validate_neutrino_rates(&candidate, NULL), ghl_success,
        "zero transport-opacity sum", 4055);

  candidate = rates;
  candidate.kappa_a_N = DBL_MAX;
  candidate.n_eq = DBL_MAX;
  candidate.mean_energy = DBL_MIN;
  candidate.J_eq = DBL_MAX * DBL_MIN;
  candidate.eta_N = DBL_MAX;
  candidate.kappa_a_N_cc = DBL_MAX;
  candidate.eta_N_cc = DBL_MAX;
  candidate.kappa_a_E = 0.0;
  candidate.eta_E = 0.0;
  candidate.kappa_s = 0.0;
  candidate.kappa_tr = 0.0;
  require_invalid_rates(&candidate, 4056);

  /* Single-species validation rejects pair or non-CC electron number data,
   * while the aggregate validator accepts a physically complete bundle. */
  candidate = rates;
  candidate.eta_N_pair[0] = 0.01;
  candidate.eta_E_pair[0] = 0.02;
  ghl_m1_neutrino_diagnostics single_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&single_diagnostics);
  require_error(
        ghl_m1_neutrino_validate_single_species_rates(&candidate, &single_diagnostics),
        ghl_error_m1_microphysics_failure, "single-species pair-rate validation", 4057);
  require_condition(
        single_diagnostics.provider_validation_failures == 1,
        "single-species partner requirement was not counted", 4057);
  candidate = rates;
  require_error(
        ghl_m1_neutrino_validate_single_species_rates(&candidate, NULL), ghl_success,
        "single-species electron rates", 4058);
  candidate.species = ghl_m1_neutrino_nux;
  candidate.lepton_weight = 0.0;
  candidate.eta_N_cc = 0.0;
  candidate.kappa_a_N_cc = 0.0;
  require_error(
        ghl_m1_neutrino_validate_single_species_rates(&candidate, NULL), ghl_success,
        "single-species heavy-flavor rates", 4059);

  double dYe = 17.0;
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(NULL, 1.0, 1.0, &dYe),
        ghl_error_m1_null_pointer, "null lepton rates", 4060);
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&rates, 1.0, 1.0, NULL),
        ghl_error_m1_null_pointer, "null lepton output", 4061);
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&rates, NAN, 1.0, &dYe),
        ghl_error_m1_invalid_state, "nonfinite lepton increment", 4062);
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&rates, 1.0, 0.0, &dYe),
        ghl_error_m1_invalid_state, "zero baryon normalization", 4063);
  candidate = rates;
  candidate.lepton_weight = NAN;
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&candidate, 1.0, 1.0, &dYe),
        ghl_error_m1_microphysics_failure, "invalid lepton weight", 4064);
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&rates, DBL_MAX, DBL_MIN, &dYe),
        ghl_error_m1_invalid_state, "overflowed lepton increment", 4065);
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&rates, -2.0, 4.0, &dYe), ghl_success,
        "signed lepton increment", 4066);
  require_close(dYe, 0.5, 0.0, 0.0, "signed lepton increment value", 4066);

  candidate = rates;
  candidate.lepton_weight = 0.0;
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&candidate, 1.0, 1.0, &dYe),
        ghl_error_m1_microphysics_failure, "wrong electron-neutrino lepton weight",
        4067);
  candidate = rates;
  candidate.species = ghl_m1_neutrino_anue;
  candidate.lepton_weight = 0.0;
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&candidate, 1.0, 1.0, &dYe),
        ghl_error_m1_microphysics_failure, "wrong electron-antineutrino lepton weight",
        4068);
  candidate = rates;
  candidate.species = ghl_m1_neutrino_nux;
  candidate.lepton_weight = 1.0;
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&candidate, 1.0, 1.0, &dYe),
        ghl_error_m1_microphysics_failure, "wrong heavy-flavor lepton weight", 4069);
  candidate = rates;
  candidate.species = ghl_m1_neutrino_nux;
  candidate.lepton_weight = 0.0;
  candidate.eta_N_cc = 0.0;
  candidate.kappa_a_N_cc = 0.0;
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&candidate, 1.0, 1.0, &dYe),
        ghl_error_m1_microphysics_failure,
        "heavy-flavor nonzero charged-current exchange", 4070);
  require_error(
        ghl_m1_compute_neutrino_lepton_increment(&candidate, -0.0, 1.0, &dYe),
        ghl_success, "heavy-flavor signed-zero charged-current exchange", 4071);
  require_close(dYe, 0.0, 0.0, 0.0, "heavy-flavor zero composition source", 4071);

  double N_out = -1.0;
  bool floor_applied = false;
  require_error(
        ghl_m1_apply_neutrino_number_floor(NULL, 0.0, &N_out, &floor_applied),
        ghl_error_m1_null_pointer, "null number-floor parameters", 4070);
  require_error(
        ghl_m1_apply_neutrino_number_floor(&nu_params, 0.0, NULL, &floor_applied),
        ghl_error_m1_null_pointer, "null number-floor output", 4071);
  ghl_m1_neutrino_parameters bad_nu_params = nu_params;
  bad_nu_params.N_floor = NAN;
  require_error(
        ghl_m1_apply_neutrino_number_floor(&bad_nu_params, 0.0, &N_out, &floor_applied),
        ghl_error_m1_invalid_state, "invalid number floor", 4072);
  require_error(
        ghl_m1_apply_neutrino_number_floor(&nu_params, NAN, &N_out, &floor_applied),
        ghl_error_m1_invalid_state, "nonfinite number input", 4073);
  require_error(
        ghl_m1_apply_neutrino_number_floor(&nu_params, 0.0, &N_out, NULL), ghl_success,
        "number floor without optional flag", 4074);
  require_close(N_out, nu_params.N_floor, 0.0, 0.0, "applied number floor", 4074);
  require_error(
        ghl_m1_apply_neutrino_number_floor(&nu_params, 1.0, &N_out, &floor_applied),
        ghl_success, "unapplied number floor", 4075);
  require_condition(
        !floor_applied && N_out == 1.0,
        "number floor changed an already-admissible value", 4075);
  ghl_m1_neutrino_diagnostics_initialize(NULL);

  ghl_m1_neutrino_state state_in, state_out;
  make_state(rng, &metric, &state_in);
  state_out = state_in;
  state_out.E += 0.1;
  state_out.F[0] += 0.02;
  ghl_m1_neutrino_exchange exchange = { 0 };

  /* The explicit-thin helper is also a public caller boundary. Exercise each
   * required pointer independently; the helper must reject these before it
   * initializes the transactional outputs or reads a later argument. */
  int thin_selected = 1;
  ghl_m1_neutrino_state thin_state = state_in;
  ghl_m1_neutrino_exchange thin_exchange = { .dE_rad = -1.0 };
  ghl_m1_neutrino_diagnostics thin_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&thin_diagnostics);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              NULL, &nu_params, &metric, &prims, &rates, 0.1, 2.0, &state_in,
              &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_null_pointer, "null thin M1 parameters", 4500);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, NULL, &metric, &prims, &rates, 0.1, 2.0, &state_in,
              &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_null_pointer, "null thin neutrino parameters", 4501);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, NULL, &prims, &rates, 0.1, 2.0, &state_in,
              &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_null_pointer, "null thin metric", 4502);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, NULL, &rates, 0.1, 2.0, &state_in,
              &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_null_pointer, "null thin primitives", 4503);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, NULL, 0.1, 2.0, &state_in,
              &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_null_pointer, "null thin rates", 4504);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &rates, 0.1, 2.0, NULL,
              &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_null_pointer, "null thin input state", 4505);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &rates, 0.1, 2.0, &state_in, NULL,
              &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_null_pointer, "null thin selection flag", 4506);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &rates, 0.1, 2.0, &state_in,
              &thin_selected, NULL, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_null_pointer, "null thin state output", 4507);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &rates, 0.1, 2.0, &state_in,
              &thin_selected, &thin_state, NULL, &thin_diagnostics),
        ghl_error_m1_null_pointer, "null thin exchange output", 4508);

  /* A thick scattering opacity reaches the second inequality after the
   * absorption side passes. */
  ghl_m1_neutrino_rates thick_scattering = rates;
  thick_scattering.kappa_a_E = 0.10;
  thick_scattering.kappa_s = 20.0;
  thick_scattering.kappa_tr = 20.10;
  thick_scattering.eta_E = thick_scattering.kappa_a_E * thick_scattering.J_eq;
  thin_state = state_in;
  thin_exchange = (ghl_m1_neutrino_exchange){ .dE_rad = -3.0 };
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &thick_scattering, 0.1, 2.0,
              &state_in, &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_success, "thin scattering non-selection", 4509);
  require_condition(
        thin_selected == 0 && memcmp(&thin_state, &state_in, sizeof(thin_state)) == 0,
        "scattering non-selection changed state", 4509);
  require_zero_exchange(&thin_exchange, "scattering non-selection exchange", 4509);

  /* An invalid provider rate and a non-finite lapse are rejected after the
   * output transaction is initialized, preserving the input state and zero
   * exchange. */
  ghl_m1_neutrino_rates bad_thin_rates = rates;
  bad_thin_rates.eta_E = -DBL_MIN;
  thin_state = state_in;
  thin_exchange = (ghl_m1_neutrino_exchange){ .dE_rad = -4.0 };
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &bad_thin_rates, 0.1, 2.0,
              &state_in, &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_microphysics_failure, "invalid thin rates", 4510);
  require_condition(
        memcmp(&thin_state, &state_in, sizeof(thin_state)) == 0,
        "invalid thin rates changed state", 4510);
  require_zero_exchange(&thin_exchange, "invalid thin rates exchange", 4510);
  ghl_metric_quantities bad_thin_metric = metric;
  bad_thin_metric.lapse = NAN;
  thin_state = state_in;
  thin_exchange = (ghl_m1_neutrino_exchange){ .dE_rad = -5.0 };
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &bad_thin_metric, &prims, &rates, 0.1, 2.0,
              &state_in, &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_invalid_state, "invalid thin lapse", 4511);
  require_condition(
        memcmp(&thin_state, &state_in, sizeof(thin_state)) == 0,
        "invalid thin lapse changed state", 4511);
  require_zero_exchange(&thin_exchange, "invalid thin lapse exchange", 4511);

  require_error(
        ghl_m1_neutrino_assemble_exchange(
              NULL, &state_out, &rates, 0.1, 1.0, 2.0, &exchange),
        ghl_error_m1_null_pointer, "null exchange input state", 4080);
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, NULL, &rates, 0.1, 1.0, 2.0, &exchange),
        ghl_error_m1_null_pointer, "null exchange output state", 4081);
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &state_out, NULL, 0.1, 1.0, 2.0, &exchange),
        ghl_error_m1_null_pointer, "null exchange rates", 4082);
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &state_out, &rates, 0.1, 1.0, 2.0, NULL),
        ghl_error_m1_null_pointer, "null exchange output", 4083);
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &state_out, &rates, 0.1, 0.0, 2.0, &exchange),
        ghl_error_m1_invalid_state, "nonpositive exchange determinant", 4084);
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &state_out, &rates, NAN, 1.0, 2.0, &exchange),
        ghl_error_m1_invalid_state, "nonfinite lepton exchange", 4085);
  ghl_m1_neutrino_state bad_state = state_out;
  bad_state.N = NAN;
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &bad_state, &rates, 0.1, 1.0, 2.0, &exchange),
        ghl_error_m1_invalid_state, "nonfinite number exchange", 4086);
  bad_state = state_out;
  bad_state.E = NAN;
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &bad_state, &rates, 0.1, 1.0, 2.0, &exchange),
        ghl_error_m1_invalid_state, "nonfinite energy exchange", 4087);
  bad_state = state_out;
  bad_state.F[0] = NAN;
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &bad_state, &rates, 0.1, 1.0, 2.0, &exchange),
        ghl_error_m1_invalid_state, "nonfinite flux exchange", 4088);
  bad_state = state_out;
  bad_state.E = DBL_MAX;
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &bad_state, &rates, 0.1, DBL_MAX, 2.0, &exchange),
        ghl_error_m1_invalid_state, "overflowed matter energy exchange", 4089);
  bad_state = state_out;
  bad_state.F[0] = DBL_MAX;
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &bad_state, &rates, 0.1, DBL_MAX, 2.0, &exchange),
        ghl_error_m1_invalid_state, "overflowed matter momentum exchange", 4090);
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &state_out, &rates, DBL_MAX, 1.0, DBL_MIN, &exchange),
        ghl_error_m1_invalid_state, "overflowed electron exchange", 4091);
  candidate = rates;
  candidate.lepton_weight = 0.0;
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &state_out, &candidate, 0.1, 1.0, 2.0, &exchange),
        ghl_error_m1_microphysics_failure, "exchange lepton policy", 4092);
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &state_out, &rates, 0.1, 1.0, 0.0, &exchange),
        ghl_error_m1_invalid_state, "exchange baryon normalization", 4093);
  require_error(
        ghl_m1_neutrino_assemble_exchange(
              &state_in, &state_out, &rates, 0.1, 1.0, 2.0, &exchange),
        ghl_success, "valid exchange assembly", 4094);
  require_exchange_contract(&state_in, &state_out, &exchange, &metric, 2.0, 4094);

  double N_update = -1.0;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              NULL, &rates, 0.1, 1.0, 1.0, &N_update),
        ghl_error_m1_null_pointer, "null BE parameters", 4100);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, NULL, 0.1, 1.0, 1.0, &N_update),
        ghl_error_m1_null_pointer, "null BE rates", 4101);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &rates, 0.1, 1.0, 1.0, NULL),
        ghl_error_m1_null_pointer, "null BE output", 4102);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &rates, NAN, 1.0, 1.0, &N_update),
        ghl_error_m1_invalid_state, "nonfinite BE timestep", 4103);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &rates, -0.1, 1.0, 1.0, &N_update),
        ghl_error_m1_invalid_state, "negative BE timestep", 4104);
  ghl_m1_neutrino_parameters zero_gamma_floor = nu_params;
  zero_gamma_floor.Gamma_N_floor = 0.0;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &zero_gamma_floor, &rates, 0.1, 1.0, 1.0, &N_update),
        ghl_success, "default BE Gamma floor", 4105);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &rates, 0.1, nu_params.Gamma_N_floor, 1.0, &N_update),
        ghl_error_m1_invalid_state, "BE Gamma floor boundary", 4106);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &rates, 0.1, NAN, 1.0, &N_update),
        ghl_error_m1_invalid_state, "nonfinite BE Gamma", 4107);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &rates, 0.1, 1.0, NAN, &N_update),
        ghl_error_m1_invalid_state, "nonfinite BE number", 4108);
  candidate = rates;
  candidate.eta_N = -1.0;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &candidate, 0.1, 1.0, 1.0, &N_update),
        ghl_error_m1_microphysics_failure, "invalid BE rates", 4109);
  candidate = rates;
  candidate.kappa_a_N = DBL_MAX;
  candidate.n_eq = 1.0;
  candidate.eta_N = DBL_MAX;
  candidate.kappa_a_N_cc = DBL_MAX;
  candidate.eta_N_cc = DBL_MAX;
  candidate.kappa_a_E = 0.0;
  candidate.kappa_s = 0.0;
  candidate.kappa_tr = 0.0;
  candidate.eta_E = 0.0;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &candidate, DBL_MAX, 1.0, 1.0, &N_update),
        ghl_success, "scaled BE denominator", 4110);
  require_close(N_update, 1.0, 2.0e-15, 0.0, "scaled BE denominator endpoint", 4110);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &candidate, 1.0, 1.0, DBL_MAX, &N_update),
        ghl_success, "scaled BE numerator", 4111);
  require_close(N_update, 2.0, 2.0e-15, 0.0, "scaled BE numerator endpoint", 4111);

  /* dt_alpha*kappa overflows before division by Gamma_N although the
   * endpoint is finite and close to the equilibrium value eta/kappa = 0.01. */
  make_rates(ghl_m1_neutrino_nue, 1.0e150, 0.0, 0.0, 0.01, 1.0, 0.0, 0.0, &candidate);
  require_error(
        ghl_m1_validate_neutrino_rates(&candidate, NULL), ghl_success,
        "large finite BE rate bundle", 4113);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &candidate, 1.0e159, 100.0, 1.0, &N_update),
        ghl_success, "BE division after overflowing product", 4113);
  require_close(N_update, 1.0, 2.0e-15, 0.0, "scaled finite BE endpoint", 4113);

  /* Scaled evaluation must still reject an endpoint whose mathematical value
   * exceeds binary64. Here the stiff limit is Gamma_N*n_eq = 2*DBL_MAX. */
  make_rates(ghl_m1_neutrino_nue, 1.0, 0.0, 0.0, DBL_MAX, 1.0, 0.0, 0.0, &candidate);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &candidate, DBL_MAX, 2.0, 0.0, &N_update),
        ghl_error_m1_invalid_state, "nonrepresentable BE endpoint", 4114);
  candidate = rates;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &candidate, 0.0, 1.0, 1.0, &N_update),
        ghl_success, "zero-step BE update", 4112);
  require_close(N_update, 1.0, 0.0, 0.0, "zero-step BE number", 4112);

  /* Repair and EN-bound APIs expose independent validation order and optional
   * diagnostics behavior. */
  ghl_m1_neutrino_diagnostics repair_nd;
  ghl_m1_neutrino_diagnostics_initialize(&repair_nd);
  ghl_m1_neutrino_state repair_state = state_in;
  require_error(
        ghl_m1_repair_neutrino_state(
              NULL, &nu_params, &metric, &repair_state, &repair_nd),
        ghl_error_m1_null_pointer, "null repair M1 parameters", 4120);
  require_error(
        ghl_m1_repair_neutrino_state(
              m1_params, NULL, &metric, &repair_state, &repair_nd),
        ghl_error_m1_null_pointer, "null repair neutrino parameters", 4121);
  require_error(
        ghl_m1_repair_neutrino_state(
              m1_params, &nu_params, NULL, &repair_state, &repair_nd),
        ghl_error_m1_null_pointer, "null repair metric", 4122);
  require_error(
        ghl_m1_repair_neutrino_state(m1_params, &nu_params, &metric, NULL, &repair_nd),
        ghl_error_m1_null_pointer, "null repair state", 4123);
  bad_nu_params = nu_params;
  bad_nu_params.N_floor = NAN;
  require_error(
        ghl_m1_repair_neutrino_state(
              m1_params, &bad_nu_params, &metric, &repair_state, &repair_nd),
        ghl_error_m1_invalid_state, "invalid repair floor", 4125);
  repair_state = state_in;
  repair_state.N = NAN;
  require_error(
        ghl_m1_repair_neutrino_state(
              m1_params, &nu_params, &metric, &repair_state, &repair_nd),
        ghl_error_m1_invalid_state, "nonfinite repair number", 4126);
  repair_state = state_in;
  repair_state.F[1] = NAN;
  require_error(
        ghl_m1_repair_neutrino_state(
              m1_params, &nu_params, &metric, &repair_state, &repair_nd),
        ghl_error_m1_invalid_state, "nonfinite repair flux", 4127);
  repair_state = state_in;
  const ghl_m1_neutrino_state repair_before = repair_state;
  require_error(
        ghl_m1_repair_neutrino_state(
              m1_params, &nu_params, &metric, &repair_state, NULL),
        ghl_success, "repair without diagnostics", 4128);
  require_condition(
        memcmp(&repair_state, &repair_before, sizeof(repair_state)) == 0,
        "no-op repair changed an admissible state", 4128);

  ghl_m1_neutrino_current current = { 0 };
  current.J = 2.0;
  current.Gamma_N = 2.0;
  ghl_m1_neutrino_state bound_state = state_in;
  bound_state.N = 2.0;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(NULL, &nu_params, &current),
        ghl_error_m1_null_pointer, "null EN state", 4130);
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, NULL, &current),
        ghl_error_m1_null_pointer, "null EN parameters", 4131);
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &nu_params, NULL),
        ghl_error_m1_null_pointer, "null EN current", 4132);
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &nu_params, &current), ghl_success,
        "disabled EN bounds", 4133);
  ghl_m1_neutrino_parameters bounded_nu_params = nu_params;
  bounded_nu_params.enforce_mean_energy_bounds = 1;
  bounded_nu_params.mean_energy_min = NAN;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_error_m1_invalid_state, "nonfinite EN lower bound", 4134);
  bounded_nu_params.mean_energy_min = 3.0;
  bounded_nu_params.mean_energy_max = 2.0;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_error_m1_invalid_state, "reversed EN bounds", 4135);
  bounded_nu_params.mean_energy_min = 0.0;
  bounded_nu_params.mean_energy_max = 0.0;
  bound_state.N = NAN;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_error_m1_invalid_state, "nonfinite EN number", 4136);
  bound_state.N = 0.5 * nu_params.N_floor;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_error_m1_invalid_state, "below-floor EN number", 4137);
  bound_state.N = 2.0;
  current.J = NAN;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_error_m1_invalid_state, "nonfinite EN current", 4138);
  current.J = DBL_MAX;
  current.Gamma_N = DBL_MAX;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_error_m1_invalid_state, "overflowed EN ratio", 4139);
  current.J = 2.0;
  current.Gamma_N = 2.0;
  bounded_nu_params.mean_energy_min = 3.0;
  bounded_nu_params.mean_energy_max = 4.0;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_error_m1_invalid_state, "low EN ratio", 4140);
  bounded_nu_params.mean_energy_min = 0.0;
  bounded_nu_params.mean_energy_max = 1.0;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_error_m1_invalid_state, "high EN ratio", 4141);
  bounded_nu_params.mean_energy_max = 4.0;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_success, "valid EN ratio", 4142);
  bounded_nu_params.N_floor = 0.0;
  bound_state.N = 0.0;
  require_error(
        ghl_m1_neutrino_check_EN_bounds(&bound_state, &bounded_nu_params, &current),
        ghl_success, "zero-number EN skip", 4143);

  /* Public source boundaries are exercised with a known valid closure input;
   * the tests below focus on their own pointer/state/configuration guards. */
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(&state_in);
  ghl_m1_sources sources = { 0 };
  double N_source = 0.0;
  bool fallback = false;
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              NULL, &metric, &prims, &rad_state, &rates, &fallback, &sources),
        ghl_error_m1_null_pointer, "null diagnostic source parameters", 4150);
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              m1_params, NULL, &prims, &rad_state, &rates, &fallback, &sources),
        ghl_error_m1_null_pointer, "null diagnostic source metric", 4151);
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              m1_params, &metric, NULL, &rad_state, &rates, &fallback, &sources),
        ghl_error_m1_null_pointer, "null diagnostic source primitives", 4152);
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              m1_params, &metric, &prims, NULL, &rates, &fallback, &sources),
        ghl_error_m1_null_pointer, "null diagnostic source radiation", 4153);
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              m1_params, &metric, &prims, &rad_state, NULL, &fallback, &sources),
        ghl_error_m1_null_pointer, "null diagnostic source rates", 4154);
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              m1_params, &metric, &prims, &rad_state, &rates, &fallback, NULL),
        ghl_error_m1_null_pointer, "null diagnostic source output", 4155);
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              m1_params, &metric, &prims, &rad_state, &rates, NULL, &sources),
        ghl_success, "diagnostic source without fallback flag", 4156);

  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, &prims, &state_in, &rates, &sources,
              &N_source),
        ghl_success, "validated public interaction source", 4157);
  ghl_m1_neutrino_parameters invalid_source_nu = nu_params;
  invalid_source_nu.N_floor = NAN;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &invalid_source_nu, &metric, &prims, &state_in, &rates,
              &sources, &N_source),
        ghl_error_m1_invalid_state, "invalid public source number floor", 4158);
  ghl_m1_neutrino_rates invalid_source_rates = rates;
  invalid_source_rates.eta_E = NAN;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, &prims, &state_in, &invalid_source_rates,
              &sources, &N_source),
        ghl_error_m1_microphysics_failure, "invalid public source rates", 4158);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, &prims, &state_in, NULL, &rates, &sources,
              &N_source),
        ghl_error_m1_null_pointer, "null source closure", 4159);
  ghl_m1_closure closure;
  require_error(
        ghl_m1_compute_neutrino_closure(m1_params, &metric, &prims, &state_in, &closure),
        ghl_success, "source closure construction", 4160);
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, &prims, &state_in, &closure,
              &invalid_source_rates, &sources, &N_source),
        ghl_error_m1_microphysics_failure, "invalid closure-source rates", 4161);
  ghl_m1_closure invalid_source_closure = closure;
  invalid_source_closure.P[0][0] = NAN;
  sources = (ghl_m1_sources){ .S_E = -76.0, .S = { -77.0, -78.0, -79.0 } };
  N_source = -80.0;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources_from_closure(
              m1_params, &nu_params, &metric, &prims, &state_in, &invalid_source_closure,
              &rates, &sources, &N_source),
        ghl_error_m1_invalid_state, "invalid closure-source tensor", 4162);
  if(sources.S_E != -76.0 || sources.S[0] != -77.0 || N_source != -80.0) {
    ghl_error("M1 source-update case 4163: invalid closure-source published outputs\n");
  }

  thin_selected = 0;
  thin_exchange = (ghl_m1_neutrino_exchange){ 0 };
  thin_state = state_in;
  ghl_m1_neutrino_diagnostics_initialize(&thin_diagnostics);
  ghl_m1_neutrino_rates thick_rates = rates;
  thick_rates.kappa_a_E = 20.0;
  thick_rates.kappa_s = 20.0;
  thick_rates.kappa_tr = 40.0;
  thick_rates.eta_E = thick_rates.kappa_a_E * thick_rates.J_eq;
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &thick_rates, 0.1, 2.0, &state_in,
              &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_success, "thin non-selection", 4160);
  require_condition(
        thin_selected == 0 && memcmp(&thin_state, &state_in, sizeof(thin_state)) == 0,
        "thin non-selection mutated state", 4160);
  require_zero_exchange(&thin_exchange, "thin non-selection exchange", 4160);
  require_error(
        ghl_m1_try_neutrino_explicit_thin_update_with_diagnostics(
              m1_params, &nu_params, &metric, &prims, &rates, NAN, 2.0, &state_in,
              &thin_selected, &thin_state, &thin_exchange, &thin_diagnostics),
        ghl_error_m1_invalid_state, "thin nonfinite timestep", 4161);
}

static void test_implicit_kernel_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  make_metric(rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 0.10, 0.20, 0.10, 1.0, 2.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_state state;
  make_state(rng, &metric, &state);
  const double U_base[4] = { state.E, state.F[0], state.F[1], state.F[2] };
  const double U[4] = { state.E, state.F[0], state.F[1], state.F[2] };
  double residual[4] = { 0.0, 0.0, 0.0, 0.0 };

  ghl_m1_rad_state trial = { 0 };
  require_error(
        ghl_m1_neutrino_build_trial_state(NULL, U, &trial), ghl_error_m1_null_pointer,
        "null trial metric", 4200);
  require_error(
        ghl_m1_neutrino_build_trial_state(&metric, NULL, &trial),
        ghl_error_m1_null_pointer, "null trial conservative state", 4201);
  require_error(
        ghl_m1_neutrino_build_trial_state(&metric, U, NULL), ghl_error_m1_null_pointer,
        "null trial output", 4202);
  ghl_metric_quantities bad_metric = metric;
  bad_metric.sqrt_detgamma = 0.0;
  require_error(
        ghl_m1_neutrino_build_trial_state(&bad_metric, U, &trial),
        ghl_error_m1_invalid_metric, "invalid trial metric", 4203);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = NAN;
  require_error(
        ghl_m1_neutrino_build_trial_state(&bad_metric, U, &trial),
        ghl_error_m1_invalid_metric, "nonfinite trial metric", 4203);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = -1.0;
  require_error(
        ghl_m1_neutrino_build_trial_state(&bad_metric, U, &trial),
        ghl_error_m1_invalid_metric, "negative trial metric", 4203);
  double bad_U[4] = { U[0], U[1], U[2], U[3] };
  bad_U[0] = NAN;
  require_error(
        ghl_m1_neutrino_build_trial_state(&metric, bad_U, &trial),
        ghl_error_m1_implicit_admissibility, "nonfinite trial energy", 4204);
  bad_U[0] = U[0];
  bad_U[2] = NAN;
  require_error(
        ghl_m1_neutrino_build_trial_state(&metric, bad_U, &trial),
        ghl_error_m1_implicit_admissibility, "nonfinite trial flux", 4205);
  require_error(
        ghl_m1_neutrino_build_trial_state(&metric, U, &trial), ghl_success,
        "valid trial state", 4206);
  require_close(trial.E, state.E, 0.0, 0.0, "trial energy", 4206);

  require_error(
        ghl_m1_neutrino_check_trial_admissibility(NULL, &metric, &trial),
        ghl_error_m1_null_pointer, "null admissibility parameters", 4210);
  require_error(
        ghl_m1_neutrino_check_trial_admissibility(m1_params, NULL, &trial),
        ghl_error_m1_null_pointer, "null admissibility metric", 4211);
  require_error(
        ghl_m1_neutrino_check_trial_admissibility(m1_params, &metric, NULL),
        ghl_error_m1_null_pointer, "null admissibility state", 4212);
  ghl_m1_rad_state inadmissible_trial = trial;
  inadmissible_trial.E = NAN;
  require_error(
        ghl_m1_neutrino_check_trial_admissibility(
              m1_params, &metric, &inadmissible_trial),
        ghl_error_m1_implicit_admissibility, "nonfinite trial admissibility state",
        4213);
  require_error(
        ghl_m1_neutrino_check_trial_admissibility(m1_params, &metric, &trial),
        ghl_success, "valid trial admissibility", 4214);

  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              NULL, &metric, &prims, &rates, 0.1, U_base, U, residual),
        ghl_error_m1_null_pointer, "null residual parameters", 4220);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, NULL, &prims, &rates, 0.1, U_base, U, residual),
        ghl_error_m1_null_pointer, "null residual metric", 4221);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, NULL, &rates, 0.1, U_base, U, residual),
        ghl_error_m1_null_pointer, "null residual primitives", 4222);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, NULL, 0.1, U_base, U, residual),
        ghl_error_m1_null_pointer, "null residual rates", 4223);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, &rates, 0.1, NULL, U, residual),
        ghl_error_m1_null_pointer, "null residual base", 4224);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, &rates, 0.1, U_base, NULL, residual),
        ghl_error_m1_null_pointer, "null residual trial", 4225);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, &rates, 0.1, U_base, U, NULL),
        ghl_error_m1_null_pointer, "null residual output", 4226);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, &rates, NAN, U_base, U, residual),
        ghl_error_m1_invalid_state, "nonfinite residual timestep", 4227);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &bad_metric, &prims, &rates, 0.1, U_base, U, residual),
        ghl_error_m1_invalid_metric, "invalid residual metric", 4228);
  bad_metric = metric;
  bad_metric.lapse = NAN;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &bad_metric, &prims, &rates, 0.1, U_base, U, residual),
        ghl_error_m1_invalid_metric, "nonfinite residual lapse", 4228);
  bad_metric = metric;
  bad_metric.lapse = 0.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &bad_metric, &prims, &rates, 0.1, U_base, U, residual),
        ghl_error_m1_invalid_metric, "nonpositive residual lapse", 4228);
  bad_metric = metric;
  bad_metric.lapse = -1.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &bad_metric, &prims, &rates, 0.1, U_base, U, residual),
        ghl_error_m1_invalid_metric, "negative residual lapse", 4228);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = NAN;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &bad_metric, &prims, &rates, 0.1, U_base, U, residual),
        ghl_error_m1_invalid_metric, "nonfinite residual volume", 4228);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = -1.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &bad_metric, &prims, &rates, 0.1, U_base, U, residual),
        ghl_error_m1_invalid_metric, "negative residual volume", 4228);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = 0.0;
  double bad_base[4] = { U_base[0], U_base[1], U_base[2], U_base[3] };
  bad_base[0] = NAN;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, &rates, 0.1, bad_base, U, residual),
        ghl_error_m1_invalid_state, "nonfinite residual base", 4229);
  bad_U[0] = NAN;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, &rates, 0.1, U_base, bad_U, residual),
        ghl_error_m1_implicit_admissibility, "nonfinite residual trial", 4230);
  bad_U[0] = U[0];
  ghl_m1_neutrino_rates bad_rates = rates;
  bad_rates.eta_N = -1.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, &bad_rates, 0.1, U_base, U, residual),
        ghl_error_m1_microphysics_failure, "invalid residual rates", 4232);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, &rates, 0.0, U_base, U, residual),
        ghl_success, "valid implicit residual", 4233);
  for(int i = 0; i < 4; ++i) {
    require_close(residual[i], 0.0, 0.0, 0.0, "zero-step implicit residual", 4233);
  }

  require_error(
        ghl_m1_neutrino_compute_implicit_residual(
              m1_params, &nu_params, NULL, &prims, &rates, &state, 0.1, U, residual),
        ghl_error_m1_null_pointer, "null public residual metric", 4240);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual(
              m1_params, &nu_params, &metric, &prims, &rates, NULL, 0.1, U, residual),
        ghl_error_m1_null_pointer, "null public residual input state", 4241);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = 0.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual(
              m1_params, &nu_params, &bad_metric, &prims, &rates, &state, 0.1, U,
              residual),
        ghl_error_m1_invalid_metric, "invalid public residual determinant", 4242);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = NAN;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual(
              m1_params, &nu_params, &bad_metric, &prims, &rates, &state, 0.1, U,
              residual),
        ghl_error_m1_invalid_metric, "nonfinite public residual determinant", 4242);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = -1.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_residual(
              m1_params, &nu_params, &bad_metric, &prims, &rates, &state, 0.1, U,
              residual),
        ghl_error_m1_invalid_metric, "negative public residual determinant", 4242);
  require_error(
        ghl_m1_neutrino_compute_implicit_residual(
              m1_params, NULL, &metric, &prims, &rates, &state, 0.1, U, residual),
        ghl_error_m1_null_pointer, "null public residual neutrino parameters", 4243);

  double jacobian[4][4] = { { 0.0 } };
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              NULL, &metric, &prims, &rates, 0.1, U_base, U, residual, jacobian),
        ghl_error_m1_null_pointer, "null Jacobian parameters", 4250);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, NULL, &prims, &rates, 0.1, U_base, U, residual, jacobian),
        ghl_error_m1_null_pointer, "null Jacobian metric", 4251);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, NULL, &rates, 0.1, U_base, U, residual, jacobian),
        ghl_error_m1_null_pointer, "null Jacobian primitives", 4252);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, NULL, 0.1, U_base, U, residual, jacobian),
        ghl_error_m1_null_pointer, "null Jacobian rates", 4253);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, 0.1, NULL, U, residual, jacobian),
        ghl_error_m1_null_pointer, "null Jacobian base", 4254);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, 0.1, U_base, NULL, residual, jacobian),
        ghl_error_m1_null_pointer, "null Jacobian trial", 4255);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, 0.1, U_base, U, NULL, jacobian),
        ghl_error_m1_null_pointer, "null Jacobian residual", 4256);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, 0.1, U_base, U, residual, NULL),
        ghl_error_m1_null_pointer, "null Jacobian output", 4257);
  bad_base[0] = NAN;
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, 0.1, bad_base, U, residual, jacobian),
        ghl_error_m1_invalid_state, "nonfinite Jacobian base", 4258);
  bad_base[0] = U_base[0];
  bad_U[0] = NAN;
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, 0.1, U_base, bad_U, residual,
              jacobian),
        ghl_error_m1_invalid_state, "nonfinite Jacobian trial", 4259);
  bad_U[0] = U[0];
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, 0.0, U_base, U, residual, jacobian),
        ghl_success, "valid implicit Jacobian", 4260);
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      require_close(
            jacobian[i][j], i == j ? 1.0 : 0.0, 2.0e-8, 2.0e-10,
            "zero-step implicit Jacobian identity", 4260);
    }
  }

  /* Put one flux component just inside the realizability cone. The forward FD
   * perturbation crosses the cone, so the documented backward one-sided
   * fallback is selected and still produces a valid Jacobian. */
  double near_cone_base[4] = { 1.0, 0.0, 0.0, 0.0 };
  const double permitted_flux = sqrt(1.0 - m1_params->epsilon_c);
  const double near_cone_delta = ghl_m1_compute_fd_delta(
        m1_params, &metric, near_cone_base[1], permitted_flux, near_cone_base[0],
        near_cone_base[0]);
  double near_cone_U[4] = { 1.0, permitted_flux - 0.5 * near_cone_delta, 0.0, 0.0 };
  double forward_cone_U[4] = { near_cone_U[0], near_cone_U[1] + near_cone_delta,
                               near_cone_U[2], near_cone_U[3] };
  double backward_cone_U[4] = { near_cone_U[0], near_cone_U[1] - near_cone_delta,
                                near_cone_U[2], near_cone_U[3] };
  ghl_m1_rad_state forward_cone_state, backward_cone_state;
  require_error(
        ghl_m1_neutrino_build_trial_state(&metric, forward_cone_U, &forward_cone_state),
        ghl_success, "finite forward cone trial", 4261);
  require_error(
        ghl_m1_neutrino_check_trial_admissibility(
              m1_params, &metric, &forward_cone_state),
        ghl_error_m1_implicit_admissibility,
        "configured FD forward perturbation crossed cone", 4261);
  require_error(
        ghl_m1_neutrino_build_trial_state(
              &metric, backward_cone_U, &backward_cone_state),
        ghl_success, "finite backward cone trial", 4261);
  require_error(
        ghl_m1_neutrino_check_trial_admissibility(
              m1_params, &metric, &backward_cone_state),
        ghl_success, "configured FD backward perturbation stayed admissible", 4261);
  double near_cone_residual[4] = { 0.0, 0.0, 0.0, 0.0 };
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, &metric, &prims, &rates, 0.0, near_cone_base, near_cone_U,
              near_cone_residual),
        ghl_success, "near-cone residual", 4261);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, 0.0, near_cone_base, near_cone_U,
              near_cone_residual, jacobian),
        ghl_success, "one-sided implicit Jacobian", 4261);
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      require_close(
            jacobian[i][j], i == j ? 1.0 : 0.0, 2.0e-8, 2.0e-10,
            "one-sided zero-step Jacobian identity", 4261);
    }
  }

  /* If the forward FD trial is outside the cone while the backward trial is
   * admissible, an ordinary residual failure in the backward trial is
   * propagated instead of being relabeled as a Jacobian failure. */
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, DBL_MAX, near_cone_base, near_cone_U,
              near_cone_residual, jacobian),
        ghl_error_m1_invalid_state, "backward implicit Jacobian residual failure", 4262);

  /* Set both FD trials outside the cone. The forward admissibility error is
   * followed by the backward admissibility error, which is the only route
   * that publishes an invalid implicit Jacobian. */
  double both_cone_U[4] = { 1.0, permitted_flux + 1.5 * near_cone_delta, 0.0, 0.0 };
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              m1_params, &metric, &prims, &rates, 0.0, near_cone_base, both_cone_U,
              near_cone_residual, jacobian),
        ghl_error_m1_invalid_implicit_jacobian,
        "two-sided inadmissible implicit Jacobian", 4263);

  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian(
              m1_params, &nu_params, NULL, &prims, &rates, &state, 0.1, U, residual,
              jacobian),
        ghl_error_m1_null_pointer, "null public Jacobian metric", 4270);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian(
              m1_params, &nu_params, &metric, &prims, &rates, NULL, 0.1, U, residual,
              jacobian),
        ghl_error_m1_null_pointer, "null public Jacobian input state", 4271);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian(
              m1_params, NULL, &metric, &prims, &rates, &state, 0.1, U, residual,
              jacobian),
        ghl_error_m1_null_pointer, "null public Jacobian neutrino parameters", 4272);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian(
              m1_params, &nu_params, &bad_metric, &prims, &rates, &state, 0.1, U,
              residual, jacobian),
        ghl_error_m1_invalid_metric, "invalid public Jacobian metric", 4273);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = NAN;
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian(
              m1_params, &nu_params, &bad_metric, &prims, &rates, &state, 0.1, U,
              residual, jacobian),
        ghl_error_m1_invalid_metric, "nonfinite public Jacobian metric", 4273);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = -1.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian(
              m1_params, &nu_params, &bad_metric, &prims, &rates, &state, 0.1, U,
              residual, jacobian),
        ghl_error_m1_invalid_metric, "negative public Jacobian metric", 4273);
  /* Productive public entries: each validated public helper must also accept a
   * valid metric so its guard's continue edge stays executable. */
  ghl_m1_rad_state valid_trial;
  require_error(
        ghl_m1_neutrino_build_trial_state(&metric, U, &valid_trial), ghl_success,
        "productive trial state", 4273);
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian(
              m1_params, &nu_params, &metric, &prims, &rates, &state, 0.0, U, residual,
              jacobian),
        ghl_success, "productive public Jacobian", 4273);

  /* The full homogeneous solver owns a separate transactional validation
   * boundary from the residual/Jacobian helpers above. Exercise each
   * post-output-initialization rejection with otherwise valid operands. */
  ghl_m1_neutrino_state solver_output = state;
  ghl_m1_neutrino_exchange solver_exchange;
  ghl_m1_implicit_solve_diagnostics solver_diagnostics;
  ghl_m1_neutrino_diagnostics solver_nd;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &metric, &prims, &rates, -0.1, 3.0, &state,
              &solver_output, &solver_exchange, &solver_diagnostics, &solver_nd),
        ghl_error_m1_invalid_state, "negative homogeneous timestep", 4274);
  require_condition(
        memcmp(&solver_output, &state, sizeof(state)) == 0,
        "negative timestep changed homogeneous state", 4274);
  require_zero_exchange(&solver_exchange, "negative homogeneous timestep", 4274);
  require_condition(
        solver_nd.source_failures == 1, "negative timestep was not diagnosed", 4274);

  bad_metric = metric;
  bad_metric.sqrt_detgamma = 0.0;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &bad_metric, &prims, &rates, 0.1, 3.0, &state,
              &solver_output, &solver_exchange, &solver_diagnostics, &solver_nd),
        ghl_error_m1_invalid_metric, "invalid homogeneous metric", 4276);
  require_condition(
        solver_nd.source_failures == 1, "invalid homogeneous metric was not diagnosed",
        4276);
  bad_metric = metric;
  bad_metric.lapse = NAN;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &bad_metric, &prims, &rates, 0.1, 3.0, &state,
              &solver_output, &solver_exchange, &solver_diagnostics, &solver_nd),
        ghl_error_m1_invalid_metric, "nonfinite homogeneous lapse", 4276);
  require_condition(
        solver_nd.source_failures == 1, "nonfinite homogeneous lapse was not diagnosed",
        4276);
  bad_metric = metric;
  bad_metric.lapse = 0.0;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &bad_metric, &prims, &rates, 0.1, 3.0, &state,
              &solver_output, &solver_exchange, &solver_diagnostics, &solver_nd),
        ghl_error_m1_invalid_metric, "nonpositive homogeneous lapse", 4276);
  require_condition(
        solver_nd.source_failures == 1,
        "nonpositive homogeneous lapse was not diagnosed", 4276);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = NAN;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &bad_metric, &prims, &rates, 0.1, 3.0, &state,
              &solver_output, &solver_exchange, &solver_diagnostics, &solver_nd),
        ghl_error_m1_invalid_metric, "nonfinite homogeneous volume", 4276);
  require_condition(
        solver_nd.source_failures == 1, "nonfinite homogeneous volume was not diagnosed",
        4276);
  bad_metric = metric;
  bad_metric.sqrt_detgamma = -1.0;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &bad_metric, &prims, &rates, 0.1, 3.0, &state,
              &solver_output, &solver_exchange, &solver_diagnostics, &solver_nd),
        ghl_error_m1_invalid_metric, "negative homogeneous volume", 4276);
  require_condition(
        solver_nd.source_failures == 1, "negative homogeneous volume was not diagnosed",
        4276);

  /* solve_diagnostics is a required output: a hard failure must leave a
   * defined "no implicit solve" record rather than the caller's prior bytes.
   * Poison the struct first so an untouched buffer is detectable. */
  {
    ghl_m1_implicit_solve_diagnostics poisoned = { .newton_iterations = 1234,
                                                   .line_search_backtracks = 5678,
                                                   .fallback_substeps = 91,
                                                   .used_fallback_substepping = true,
                                                   .residual_max_norm = -3.5,
                                                   .residual_scaled_norm = -7.25,
                                                   .solution_path_flags = 0xfeedu };
    ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
    require_error(
          ghl_m1_solve_neutrino_implicit_homogeneous_update(
                m1_params, &nu_params, &bad_metric, &prims, &rates, 0.1, 3.0, &state,
                &solver_output, &solver_exchange, &poisoned, &solver_nd),
          ghl_error_m1_invalid_metric, "poisoned solve diagnostics hard failure", 4277);
    require_condition(
          poisoned.newton_iterations == 0 && poisoned.line_search_backtracks == 0
                && poisoned.fallback_substeps == 1 && !poisoned.used_fallback_substepping
                && poisoned.residual_max_norm == INFINITY
                && poisoned.residual_scaled_norm == INFINITY
                && poisoned.solution_path_flags == 0u,
          "hard failure left solve diagnostics uninitialized", 4277);
  }

  bad_metric = metric;
  bad_metric.lapse = DBL_MAX;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &bad_metric, &prims, &rates, DBL_MAX, 3.0, &state,
              &solver_output, &solver_exchange, &solver_diagnostics, &solver_nd),
        ghl_error_m1_invalid_state, "overflowed homogeneous dt alpha", 4277);
  require_condition(
        solver_nd.source_failures == 1,
        "overflowed homogeneous dt alpha was not diagnosed", 4277);

  ghl_m1_neutrino_state bad_solver_state = state;
  bad_solver_state.N = NAN;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &metric, &prims, &rates, 0.1, 3.0,
              &bad_solver_state, &solver_output, &solver_exchange, &solver_diagnostics,
              &solver_nd),
        ghl_error_m1_invalid_state, "nonfinite homogeneous number", 4278);
  bad_solver_state = state;
  bad_solver_state.F[0] = NAN;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &metric, &prims, &rates, 0.1, 3.0,
              &bad_solver_state, &solver_output, &solver_exchange, &solver_diagnostics,
              &solver_nd),
        ghl_error_m1_invalid_state, "nonfinite homogeneous flux", 4279);

  ghl_m1_neutrino_rates bad_solver_rates = rates;
  bad_solver_rates.eta_N = -DBL_MIN;
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &metric, &prims, &bad_solver_rates, 0.1, 3.0,
              &state, &solver_output, &solver_exchange, &solver_diagnostics, &solver_nd),
        ghl_error_m1_microphysics_failure, "invalid homogeneous rates", 4289);

  /* A finite, validated rate bundle can still make the interaction source
   * nonrepresentable at the trial state. The schedule must publish that
   * ordinary non-retryable error immediately rather than treating it as an
   * admissible substep failure. */
  ghl_m1_neutrino_rates nonretryable_rates = rates;
  nonretryable_rates.n_eq = 1.0;
  nonretryable_rates.J_eq = 1.0;
  nonretryable_rates.mean_energy = 1.0;
  nonretryable_rates.eta_E = DBL_MAX;
  nonretryable_rates.kappa_a_E = DBL_MAX;
  nonretryable_rates.kappa_s = 0.0;
  nonretryable_rates.kappa_tr = DBL_MAX;
  const ghl_m1_neutrino_state nonretryable_state
        = { .N = 1.0, .E = 2.0, .F = { 0.5, 0.0, 0.0 } };
  ghl_m1_neutrino_diagnostics_initialize(&solver_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &metric, &prims, &nonretryable_rates, 0.1, 3.0,
              &nonretryable_state, &solver_output, &solver_exchange, &solver_diagnostics,
              &solver_nd),
        ghl_error_m1_invalid_state, "non-retryable homogeneous residual failure", 4290);
  require_condition(
        solver_nd.source_failures == 1,
        "non-retryable homogeneous failure was not diagnosed", 4290);

  ghl_m1_newton_diagnostics newton_diagnostics;
  bool closure_fallback = false;
  double U_out[4] = { 0.0, 0.0, 0.0, 0.0 };
  /* The predictor uses the frozen primitive-aware closure. If that closure
   * fails, the Newton bridge must reset its initial guess to U_in and let the
   * residual path report the same hard error. */
  ghl_primitive_quantities bad_prims = prims;
  bad_prims.vU[0] = NAN;
  bool predictor_fallback = false;
  const double predictor_U_in[4] = { 1.0, 0.0, 0.0, 0.0 };
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              m1_params, &metric, &bad_prims, &rates, 0.1, predictor_U_in, U_out,
              &newton_diagnostics, &predictor_fallback),
        ghl_error_m1_invalid_state, "Newton predictor error reset", 4273);
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              NULL, &metric, &prims, &rates, 0.1, U_base, U_out, &newton_diagnostics,
              &closure_fallback),
        ghl_error_m1_null_pointer, "null Newton parameters", 4280);
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              m1_params, NULL, &prims, &rates, 0.1, U_base, U_out, &newton_diagnostics,
              &closure_fallback),
        ghl_error_m1_null_pointer, "null Newton metric", 4281);
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              m1_params, &metric, NULL, &rates, 0.1, U_base, U_out, &newton_diagnostics,
              &closure_fallback),
        ghl_error_m1_null_pointer, "null Newton primitives", 4282);
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              m1_params, &metric, &prims, NULL, 0.1, U_base, U_out, &newton_diagnostics,
              &closure_fallback),
        ghl_error_m1_null_pointer, "null Newton rates", 4283);
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              m1_params, &metric, &prims, &rates, 0.1, NULL, U_out, &newton_diagnostics,
              &closure_fallback),
        ghl_error_m1_null_pointer, "null Newton input", 4284);
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              m1_params, &metric, &prims, &rates, 0.1, U_base, NULL, &newton_diagnostics,
              &closure_fallback),
        ghl_error_m1_null_pointer, "null Newton output", 4285);
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              m1_params, &metric, &prims, &rates, 0.1, U_base, U_out, NULL,
              &closure_fallback),
        ghl_error_m1_null_pointer, "null Newton diagnostics", 4286);
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              m1_params, &metric, &prims, &rates, 0.1, U_base, U_out,
              &newton_diagnostics, NULL),
        ghl_error_m1_null_pointer, "null Newton fallback flag", 4287);
  require_error(
        ghl_m1_neutrino_attempt_EF_newton_step(
              m1_params, &metric, &prims, &rates, 0.0, U_base, U_out,
              &newton_diagnostics, &closure_fallback),
        ghl_success, "valid Newton substep", 4288);
  for(int i = 0; i < 4; ++i) {
    require_close(U_out[i], U_base[i], 0.0, 0.0, "zero-step Newton identity", 4288);
  }

  ghl_m1_neutrino_diagnostics mean_diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&mean_diagnostics);
  ghl_m1_neutrino_state mean_state = state;
  mean_state.N = 2.0;
  ghl_m1_neutrino_current mean_current = { 0 };
  mean_current.J = 2.0;
  mean_current.Gamma_N = 2.0;
  ghl_m1_neutrino_populate_mean_energy_diagnostics(
        &mean_state, &mean_current, &nu_params, &rates, &mean_diagnostics);
  require_close(
        mean_diagnostics.mean_energy_diag, 2.0, 0.0, 0.0, "mean-energy diagnostic",
        4290);
  require_condition(
        mean_diagnostics.mean_energy_diag_invalid == 0
              && mean_diagnostics.mean_energy_consistent
              && mean_diagnostics.Jeq_over_neq_consistent,
        "consistent mean-energy diagnostic was not recorded", 4290);
  mean_diagnostics.mean_energy_diag = 9.0;
  ghl_m1_neutrino_populate_mean_energy_diagnostics(
        NULL, &mean_current, &nu_params, &rates, &mean_diagnostics);
  require_close(
        mean_diagnostics.mean_energy_diag, 9.0, 0.0, 0.0, "null mean-energy call", 4291);
  ghl_m1_neutrino_state invalid_mean_state = mean_state;
  invalid_mean_state.N = nu_params.N_floor;
  ghl_m1_neutrino_populate_mean_energy_diagnostics(
        &invalid_mean_state, &mean_current, &nu_params, &rates, &mean_diagnostics);
  require_condition(
        mean_diagnostics.mean_energy_diag_invalid == 1,
        "floor mean-energy diagnostic was not invalidated", 4292);
  invalid_mean_state = mean_state;
  mean_current.J = NAN;
  ghl_m1_neutrino_populate_mean_energy_diagnostics(
        &invalid_mean_state, &mean_current, &nu_params, &rates, &mean_diagnostics);
  require_condition(
        mean_diagnostics.mean_energy_diag_invalid == 1,
        "nonfinite current mean-energy diagnostic was accepted", 4293);
  mean_current.J = 2.0;
  rates.mean_energy = NAN;
  ghl_m1_neutrino_populate_mean_energy_diagnostics(
        &mean_state, &mean_current, &nu_params, &rates, &mean_diagnostics);
  require_condition(
        mean_diagnostics.mean_energy_diag_invalid == 0
              && !mean_diagnostics.mean_energy_consistent,
        "invalid provider mean-energy diagnostic was mishandled", 4294);
  rates.mean_energy = 2.0;
  rates.n_eq = 0.0;
  rates.J_eq = 0.0;
  ghl_m1_neutrino_populate_mean_energy_diagnostics(
        &mean_state, &mean_current, &nu_params, &rates, &mean_diagnostics);
  require_condition(
        !mean_diagnostics.Jeq_over_neq_consistent,
        "zero equilibrium target produced a ratio diagnostic", 4295);
  mean_current.J = DBL_MAX;
  mean_current.Gamma_N = DBL_MAX;
  ghl_m1_neutrino_populate_mean_energy_diagnostics(
        &mean_state, &mean_current, &nu_params, &rates, &mean_diagnostics);
  require_condition(
        mean_diagnostics.mean_energy_diag_invalid == 1,
        "overflowed mean-energy ratio was accepted", 4296);
}

static void
test_endpoint_failure_transactions(const ghl_m1_parameters *restrict m1_params) {
  m1_test_rng rng = { .state = UINT64_C(0x454e44605f494e54) };
  ghl_metric_quantities metric;
  make_metric(&rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  const ghl_m1_neutrino_state input = { .N = 1.0, .E = 1.0 };

  /* Transparent E/F isolates the finite number endpoint from a genuinely
   * overflowing exchange normalization. An absorbing E/F endpoint isolates
   * the current-floor check. */
  const double number_opacity[] = { DBL_MAX, 1.0, 0.0 };
  const double energy_opacity[] = { 0.0, 0.0, 1.0 };
  const double equilibrium_number[] = { 1.0, 16.0, 0.1 };
  const double timestep[] = { 2.0, 1.0, 1.0 };
  const double baryon_number[] = { 1.0, DBL_MIN, 1.0 };
  for(int scenario = 0; scenario < 3; ++scenario) {
    ghl_m1_neutrino_rates rates;
    make_rates(
          ghl_m1_neutrino_nue, number_opacity[scenario], energy_opacity[scenario], 0.0,
          equilibrium_number[scenario], 1.0, 0.0, 0.0, &rates);
    ghl_m1_neutrino_parameters endpoint_params = nu_params;
    /* The third endpoint has J=0.55, below the configured current floor,
     * although the input J=1 is valid and the E/F solve converges. */
    if(scenario == 2) {
      endpoint_params.J_floor = 0.9;
    }
    ghl_m1_neutrino_state output = { .N = -1.0, .E = -2.0 };
    ghl_m1_neutrino_exchange exchange = { .dN_rad_total = -3.0 };
    ghl_m1_implicit_solve_diagnostics solve_diagnostics;
    ghl_m1_initialize_implicit_solve_diagnostics(&solve_diagnostics);
    ghl_m1_neutrino_diagnostics nd;
    ghl_m1_neutrino_diagnostics_initialize(&nd);
    const ghl_error_codes_t endpoint_error
          = ghl_m1_solve_neutrino_implicit_homogeneous_update(
                m1_params, &endpoint_params, &metric, &prims, &rates, timestep[scenario],
                baryon_number[scenario], &input, &output, &exchange, &solve_diagnostics,
                &nd);
    if(scenario == 0) {
      require_error(
            endpoint_error, ghl_success, "finite implicit number endpoint", 4600);
      require_close(output.N, 1.0, 2.0e-15, 0.0, "finite implicit number", 4600);
      require_close(
            output.E, input.E, 2.0e-15, 0.0, "transparent implicit energy", 4600);
      require_zero_exchange(&exchange, "finite endpoint exchange", 4600);
      require_condition(
            nd.source_failures == 0 && nd.source_converged == 1,
            "finite endpoint diagnostics", 4600);
      continue;
    }
    require_error(
          endpoint_error, ghl_error_m1_invalid_state, "implicit endpoint rejection",
          4600 + scenario);
    require_condition(
          memcmp(&input, &output, sizeof(input)) == 0, "endpoint failure changed state",
          4600 + scenario);
    require_zero_exchange(&exchange, "endpoint failure exchange", 4600 + scenario);
    require_condition(
          nd.source_failures == 1 && nd.source_converged == 0,
          "endpoint failure diagnostics", 4600 + scenario);
  }

  /* E/F stays transparent, so Newton converges before the number endpoint
   * overflows. A late number failure must withdraw the converged E/F state. */
  ghl_m1_neutrino_rates overflowing_number_rates;
  make_rates(
        ghl_m1_neutrino_nux, DBL_MIN, 0.0, 0.0, DBL_MAX, DBL_MIN, 0.0, 0.0,
        &overflowing_number_rates);
  ghl_primitive_quantities moving_prims = prims;
  moving_prims.vU[0] = 0.9 / sqrt(metric.gammaDD[0][0]);
  moving_prims.u0 = 1.0 / sqrt(1.0 - 0.9 * 0.9);
  ghl_m1_neutrino_state moving_input = input;
  moving_input.F[0] = 0.945 / sqrt(metric.gammaDD[0][0]);
  ghl_m1_neutrino_state overflow_output = { .N = -1.0, .E = -2.0 };
  ghl_m1_neutrino_exchange overflow_exchange = { .dN_rad_total = -3.0 };
  ghl_m1_implicit_solve_diagnostics overflow_solve_diagnostics;
  ghl_m1_neutrino_diagnostics overflow_nd;
  ghl_m1_neutrino_diagnostics_initialize(&overflow_nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &metric, &moving_prims, &overflowing_number_rates,
              DBL_MAX, 1.0, &moving_input, &overflow_output, &overflow_exchange,
              &overflow_solve_diagnostics, &overflow_nd),
        ghl_error_m1_invalid_state, "late implicit number overflow", 4603);
  require_condition(
        memcmp(&moving_input, &overflow_output, sizeof(moving_input)) == 0
              && overflow_nd.source_failures == 1,
        "late number overflow published a partial state", 4603);
  require_zero_exchange(&overflow_exchange, "late number overflow exchange", 4603);

  /* Thermalized number projection can change N without a charged-current
   * source. Its optional signed-total Ye policy must still reject overflow
   * and withdraw the otherwise valid candidate. */
  ghl_m1_neutrino_source_options options = branched_options();
  options.thermalized_number_threshold = 0.0;
  options.ye_policy = ghl_m1_neutrino_ye_from_signed_total_number;
  ghl_m1_neutrino_rates projection_rates;
  make_rates(
        ghl_m1_neutrino_nue, 0.0, 0.0, 0.0, 1.0, 1.0 / 16.0, 0.0, 0.0,
        &projection_rates);
  ghl_m1_neutrino_state projected_output;
  ghl_m1_neutrino_exchange projected_exchange;
  ghl_m1_neutrino_source_diagnostics projected_diagnostics;
  ghl_m1_neutrino_diagnostics projected_nd;
  ghl_m1_neutrino_diagnostics_initialize(&projected_nd);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &projection_rates,
              &input, &input, 0.0, DBL_MIN, &projected_output, &projected_exchange,
              &projected_diagnostics, &projected_nd),
        ghl_error_m1_invalid_state, "signed projected-number overflow", 4603);
  require_condition(
        memcmp(&projected_output, &input, sizeof(input)) == 0,
        "signed-number rejection changed state", 4603);
  require_zero_exchange(&projected_exchange, "signed-number rejection exchange", 4603);
  require_condition(
        projected_nd.source_failures == 1, "signed-number rejection diagnostics", 4603);

  /* Pair updates cannot silently insert unpaired particles through the
   * independent stage's number-floor repair. Rejection is atomic for both
   * species, including the caller's diagnostic budgets. */
  ghl_m1_neutrino_parameters pair_params[2] = { nu_params, nu_params };
  ghl_m1_neutrino_rates pair_rates[2];
  const ghl_m1_neutrino_state pair_input[2] = { input, input };
  ghl_m1_neutrino_state pair_output[2];
  ghl_m1_neutrino_exchange pair_exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics nd[2];
  for(int species = 0; species < 2; ++species) {
    pair_params[species].N_floor = 0.9;
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 1.0, 0.0, 0.0, 0.1,
          1.0, 0.01, 0.0, &pair_rates[species]);
    ghl_m1_neutrino_diagnostics_initialize(&nd[species]);
  }
  const ghl_m1_neutrino_diagnostics before_nd[2] = { nd[0], nd[1] };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, pair_params, &metric, &prims, pair_rates, pair_input,
              pair_input, 1.0, 1.0, pair_output, pair_exchange, diagnostics, nd),
        ghl_error_m1_invalid_state, "pair independent number-floor rejection", 4604);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&pair_output[species], &input, sizeof(input)) == 0,
          "pair floor rejection changed state", 4604);
    require_zero_exchange(
          &pair_exchange[species], "pair floor rejection exchange", 4604);
    ghl_m1_neutrino_diagnostics expected_nd = before_nd[species];
    expected_nd.source_failures++;
    require_condition(
          memcmp(&nd[species], &expected_nd, sizeof(nd[species])) == 0,
          "pair floor rejection diagnostics were not isolated", 4604);
  }

  const ghl_m1_neutrino_diagnostics before_energy_overflow[2] = { nd[0], nd[1] };

  /* Each pair emissivity is finite and valid, but their energy sum is not
   * representable. This is a hard error, not retry exhaustion. */
  for(int species = 0; species < 2; ++species) {
    pair_params[species] = nu_params;
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 1.0,
          1.0, 0.01, DBL_MAX, &pair_rates[species]);
  }
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, pair_params, &metric, &prims, pair_rates, pair_input,
              pair_input, 1.0, 1.0, pair_output, pair_exchange, diagnostics, nd),
        ghl_error_m1_invalid_state, "pair energy-rate overflow", 4605);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&pair_output[species], &input, sizeof(input)) == 0,
          "pair rate overflow changed state", 4605);
    require_zero_exchange(&pair_exchange[species], "pair rate overflow exchange", 4605);
    ghl_m1_neutrino_diagnostics expected_nd = before_energy_overflow[species];
    expected_nd.source_failures++;
    require_condition(
          !diagnostics[species].terminal_no_update
                && memcmp(&nd[species], &expected_nd, sizeof(nd[species])) == 0,
          "pair rate overflow diagnostics were not isolated", 4605);
  }

  const ghl_m1_neutrino_diagnostics before_number_overflow[2] = { nd[0], nd[1] };

  /* A finite coefficient in every process can still overflow the shared
   * number-emission sum. Both species must retain their source bases. */
  for(int species = 0; species < 2; ++species) {
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 1.0,
          1.0, DBL_MAX, 0.0, &pair_rates[species]);
  }
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, pair_params, &metric, &prims, pair_rates, pair_input,
              pair_input, 1.0, 1.0, pair_output, pair_exchange, diagnostics, nd),
        ghl_error_m1_microphysics_failure, "pair number-rate overflow", 4606);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&pair_output[species], &input, sizeof(input)) == 0,
          "pair number overflow changed state", 4606);
    require_zero_exchange(
          &pair_exchange[species], "pair number overflow exchange", 4606);
    ghl_m1_neutrino_diagnostics expected_nd = before_number_overflow[species];
    expected_nd.source_failures++;
    require_condition(
          !diagnostics[species].terminal_no_update
                && memcmp(&nd[species], &expected_nd, sizeof(nd[species])) == 0,
          "pair number overflow diagnostics were not isolated", 4606);
  }
}

static void
test_pair_roundtrip_energy_floor(const ghl_m1_parameters *restrict m1_params) {
  ghl_m1_parameters parameters = *m1_params;
  parameters.E_floor = 1000.0;
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 3.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params[2];
  ghl_m1_neutrino_rates rates[2];
  const ghl_m1_neutrino_state input[2]
        = { { .N = 1.0, .E = 1000.0 }, { .N = 1.0, .E = 1000.0 } };
  for(int species = 0; species < 2; ++species) {
    make_neutrino_parameters(&nu_params[species]);
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 1.0,
          2.0, 0.02, 0.03, &rates[species]);
  }
  /* The pair stage divides its densitized result by sqrt(gamma). This
   * finite round trip lands one ulp below the configured floor. */
  volatile double densitized_energy = input[0].E * metric.sqrt_detgamma;
  const double roundtrip_energy = densitized_energy / metric.sqrt_detgamma;
  require_condition(
        roundtrip_energy < parameters.E_floor,
        "pair round-trip fixture did not cross energy floor", 4612);
  ghl_m1_neutrino_state output[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics nd[2] = { { 0 }, { 0 } };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              &parameters, nu_params, &metric, &prims, rates, input, input, 0.0, 1.0,
              output, exchange, diagnostics, nd),
        ghl_success, "pair round-trip floor repair", 4612);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&output[species], &input[species], sizeof(output[species])) == 0
                && nd[species].EF_repairs == 1,
          "pair round-trip repair lost state or accounting", 4612);
    require_close(
          nd[species].repair_dE, parameters.E_floor - roundtrip_energy, 0.0, 0.0,
          "pair round-trip energy repair budget", 4612);
    require_zero_exchange(&exchange[species], "pair round-trip exchange", 4612);
  }
}

static void
test_late_closure_fallback_transaction(const ghl_m1_parameters *restrict m1_params) {
  ghl_m1_parameters parameters = *m1_params;
  parameters.closure_root_max_iterations = 8;
  parameters.closure_root_residual_tolerance = 1.0;
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 2.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 0.01, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0, &rates);
  const ghl_m1_neutrino_state input = { .N = 1.0, .E = 1.0, .F = { -0.8, 0.0, 0.0 } };
  const ghl_m1_rad_state rad = ghl_m1_neutrino_project_rad_state(&input);
  ghl_m1_closure closure;
  require_error(
        ghl_m1_compute_closure_with_primitives(
              &parameters, &metric, &prims, &rad, &closure),
        ghl_success, "late fallback input closure", 4611);
  require_condition(
        closure.solve_status == ghl_m1_closure_solve_converged,
        "late fallback fixture failed its input precheck", 4611);
  ghl_m1_neutrino_source_options options = branched_options();
  ghl_m1_neutrino_state output;
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics nd = { 0 };
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, &parameters, &nu_params, &metric, &prims, &rates, &input, &input,
              1.0, 1.0, &output, &exchange, &diagnostics, &nd),
        ghl_success, "permitted late closure fallback", 4611);
  require_condition(
        diagnostics.closure_fallback_used, "endpoint did not exercise closure fallback",
        4611);
  options.allow_closure_fallback = false;
  nd = (ghl_m1_neutrino_diagnostics){ 0 };
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, &parameters, &nu_params, &metric, &prims, &rates, &input, &input,
              1.0, 1.0, &output, &exchange, &diagnostics, &nd),
        ghl_error_m1_implicit_solve_failure, "forbidden late closure fallback", 4611);
  require_condition(
        memcmp(&output, &input, sizeof(input)) == 0 && diagnostics.closure_fallback_used
              && nd.source_failures == 1 && nd.source_converged == 0,
        "late fallback rejection was not transactional", 4611);
  require_zero_exchange(&exchange, "late fallback exchange", 4611);
}

static void
test_pair_projected_backtracking(const ghl_m1_parameters *restrict m1_params) {
  /* A stiff moving-fluid pair reaction with a raised energy floor requires
   * projection and backtracking. The configured mixed tolerance is part of
   * this solve's acceptance contract, checked explicitly below. */
  ghl_m1_parameters parameters = *m1_params;
  parameters.E_floor = 0.3;
  require_error(
        ghl_m1_set_newton_tolerances(0.01, 0.001, &parameters), ghl_success,
        "pair projection tolerances", 4610);
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.5, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  prims.vU[0] = 0.95;
  ghl_m1_neutrino_parameters nu_params[2];
  ghl_m1_neutrino_rates rates[2];
  const ghl_m1_neutrino_state input[2]
        = { { .N = 1.0, .E = 1.0, .F = { 0.0, 0.0, 0.0 } },
            { .N = 1.0, .E = 1.0, .F = { 0.0, 0.0, 0.0 } } };
  for(int species = 0; species < 2; ++species) {
    make_neutrino_parameters(&nu_params[species]);
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 0.01,
          2.0, 0.02, 0.01, &rates[species]);
  }
  ghl_m1_neutrino_state output[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics nd[2] = { { 0 }, { 0 } };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              &parameters, nu_params, &metric, &prims, rates, input, input, 3.0, 1.0,
              output, exchange, diagnostics, nd),
        ghl_success, "projected pair reaction", 4610);
  for(int species = 0; species < 2; ++species) {
    const unsigned expected = ghl_m1_solution_path_projection
                              | ghl_m1_solution_path_line_search_backtracking;
    require_condition(
          (diagnostics[species].implicit.solution_path_flags & expected) == expected
                && diagnostics[species].implicit.line_search_backtracks > 0
                && diagnostics[species].implicit.residual_scaled_norm <= 1.0,
          "pair projection/backtracking diagnostics lost", 4610);
    require_condition(
          output[species].E >= parameters.E_floor,
          "projected pair endpoint violates energy floor", 4610);
    require_close(
          exchange[species].dTau_matter, -(output[species].E - input[species].E), 0.0,
          0.0, "projected pair energy exchange", 4610);
    require_close(
          exchange[species].dS_matter[0], -(output[species].F[0] - input[species].F[0]),
          0.0, 0.0, "projected pair momentum exchange", 4610);
    require_close(
          exchange[species].dYe_matter, 0.0, 0.0, 0.0, "projected pair lepton exchange",
          4610);
  }
  require_close(
        output[0].N - output[1].N, 0.0, 0.0, 0.0, "projected pair number difference",
        4610);
}

static void
test_pair_moving_closure_diagnostics(const ghl_m1_parameters *restrict m1_params) {
  m1_test_rng rng = { .state = UINT64_C(0x4d315f534f555243) };
  ghl_metric_quantities metric;
  make_metric(&rng, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  prims.vU[0] = 0.4;
  ghl_m1_neutrino_parameters nu_params[2];
  ghl_m1_neutrino_rates rates[2];
  const ghl_m1_neutrino_state input[2]
        = { { .N = 1.0, .E = 1.0 }, { .N = 1.0, .E = 1.0 } };
  ghl_m1_neutrino_state output[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics nd[2];
  for(int species = 0; species < 2; ++species) {
    make_neutrino_parameters(&nu_params[species]);
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 1.0,
          2.0, 0.02, 0.03, &rates[species]);
    ghl_m1_neutrino_diagnostics_initialize(&nd[species]);
  }
  /* The zero-flux moving-fluid closure uses the documented admissibility
   * fallback. The pair solve must preserve this diagnostic while publishing
   * a conserving, admissible endpoint for both species. */
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, input, input, 0.5, 1.0,
              output, exchange, diagnostics, nd),
        ghl_success, "moving pair source", 4607);
  for(int species = 0; species < 2; ++species) {
    require_state_finite_and_admissible(
          m1_params, &nu_params[species], &metric, &output[species], 4607);
    require_condition(
          diagnostics[species].closure_fallback_used
                && (diagnostics[species].implicit.solution_path_flags
                    & ghl_m1_solution_path_closure_fallback)
                         != 0
                && diagnostics[species].implicit.residual_scaled_norm <= 1.0,
          "moving pair lost converged fallback diagnostics", 4607);
    require_close(
          exchange[species].dL_rad_cc, 0.0, 0.0, 0.0,
          "pair-only charged-current exchange", 4607);
    require_close(
          exchange[species].dYe_matter, 0.0, 0.0, 0.0,
          "pair-only matter lepton exchange", 4607);
    require_close(
          exchange[species].dTau_matter,
          -metric.sqrt_detgamma * (output[species].E - input[species].E), 0.0, 0.0,
          "moving pair energy exchange", 4607);
    for(int direction = 0; direction < 3; ++direction) {
      require_close(
            exchange[species].dS_matter[direction],
            -metric.sqrt_detgamma
                  * (output[species].F[direction] - input[species].F[direction]),
            0.0, 0.0, "moving pair momentum exchange", 4607);
    }
  }
  require_close(
        output[0].N - output[1].N, input[0].N - input[1].N, 0.0, 0.0,
        "moving pair lepton-number conservation", 4607);
}

static void
test_number_endpoint_and_policy_overflow(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0, 1.0 / 16.0, 0.0, 0.0, 1.0 / 16.0, 0.0, 1.0 / 16.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(
        ghl_m1_neutrino_nue, 1.0 / 8.0, 0.0, 0.0, DBL_MAX, 1.0 / DBL_MAX, 0.0, 0.0,
        &rates);
  const ghl_m1_neutrino_state input = { .N = 1.0, .E = 1.0, .F = { 0.24, 0.0, 0.0 } };
  ghl_m1_neutrino_current current;
  require_error(
        ghl_m1_neutrino_derive_current(
              m1_params, &nu_params, &metric, &prims, &input, &current),
        ghl_success, "finite input number current", 4608);
  double number;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &rates, 4.0, current.Gamma_N, input.N, &number),
        ghl_success, "finite endpoint number", 4608);
  /* The physical speed is subluminal, but coordinate speed exceeds one.
   * N remains finite while its coordinate flux is no longer representable. */
  require_condition(
        isfinite(number) && current.number_transport_velocity[0] > DBL_MAX / number,
        "number endpoint did not isolate flux overflow", 4608);
  ghl_m1_neutrino_state output;
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_implicit_solve_diagnostics solve_diagnostics;
  ghl_m1_initialize_implicit_solve_diagnostics(&solve_diagnostics);
  ghl_m1_neutrino_diagnostics nd;
  ghl_m1_neutrino_diagnostics_initialize(&nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update(
              m1_params, &nu_params, &metric, &prims, &rates, 4.0, 1.0, &input, &output,
              &exchange, &solve_diagnostics, &nd),
        ghl_error_m1_invalid_state, "endpoint number-current overflow", 4608);
  require_condition(
        memcmp(&output, &input, sizeof(input)) == 0 && nd.source_failures == 1
              && nd.source_converged == 0,
        "number-current failure changed state or counters", 4608);
  require_zero_exchange(&exchange, "number-current failure exchange", 4608);

  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  const ghl_m1_neutrino_state equilibrium = { .N = 100.0, .E = 1.0 };
  make_rates(ghl_m1_neutrino_nue, 0.0, 2.0, 0.0, 1.0, 1.0, 0.0, 0.0, &rates);
  const double baryon_number = 4.0 * nextafter(0.0, 1.0);
  ghl_m1_neutrino_source_options options
        = { .policy = ghl_m1_neutrino_source_branched_compatibility,
            .thick_equilibrium_threshold = DBL_MAX,
            .scattering_threshold = DBL_MAX,
            .thermalized_number_threshold = 0.0,
            .ye_policy = ghl_m1_neutrino_ye_from_charged_current };
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics_initialize(&nd);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &equilibrium,
              &equilibrium, 1.0, baryon_number, &output, &exchange, &diagnostics, &nd),
        ghl_success, "charged-current normalization witness", 4609);
  /* Thermalized projection changes N even with zero charged-current rates.
   * A finite CC recommendation does not establish that the separately
   * selected total-number recommendation will be representable. */
  require_condition(
        isfinite(exchange.dYe_matter)
              && fabs(exchange.dN_rad_total) > DBL_MAX * baryon_number,
        "lepton policy fixture did not isolate total-number overflow", 4609);
  options.ye_policy = ghl_m1_neutrino_ye_from_signed_total_number;
  ghl_m1_neutrino_diagnostics_initialize(&nd);
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &nu_params, &metric, &prims, &rates, &equilibrium,
              &equilibrium, 1.0, baryon_number, &output, &exchange, &diagnostics, &nd),
        ghl_error_m1_invalid_state, "implicit signed-total overflow", 4609);
  require_condition(
        memcmp(&output, &equilibrium, sizeof(output)) == 0 && nd.source_failures == 1
              && diagnostics.path == ghl_m1_neutrino_source_path_hard_failure,
        "implicit lepton-policy failure publication", 4609);
  require_zero_exchange(&exchange, "implicit lepton-policy failure exchange", 4609);
}

static void
test_projected_charged_current_overflow(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, DBL_MAX, 0.0, 0.0, DBL_MIN, 1.0, 0.0, 0.0, &rates);
  require_error(
        ghl_m1_validate_neutrino_rates(&rates, NULL), ghl_success,
        "projected charged-current input rates", 4613);

  /* E/F remain unchanged, while the selected thermalized number N=J=2
   * makes kappa_a_N_cc*N overflow. The incomplete CC packet must never
   * escape the implicit solve. */
  const ghl_m1_neutrino_state input = { .N = 1.0, .E = 2.0 };
  ghl_m1_neutrino_state output = { .N = -1.0, .E = -2.0 };
  ghl_m1_neutrino_exchange exchange = { .dN_rad_total = -3.0 };
  ghl_m1_implicit_solve_diagnostics solve_diagnostics;
  ghl_m1_initialize_implicit_solve_diagnostics(&solve_diagnostics);
  ghl_m1_neutrino_diagnostics nd;
  ghl_m1_neutrino_diagnostics_initialize(&nd);
  require_error(
        ghl_m1_solve_neutrino_implicit_homogeneous_update_with_number_policy(
              m1_params, &nu_params, &metric, &prims, &rates, 1.0, 1.0, 0.0, &input,
              &output, &exchange, &solve_diagnostics, &nd),
        ghl_error_m1_invalid_state, "projected charged-current overflow", 4613);
  require_condition(
        memcmp(&output, &input, sizeof(output)) == 0 && nd.source_failures == 1
              && nd.source_converged == 0,
        "projected charged-current failure changed state or counters", 4613);
  require_zero_exchange(&exchange, "projected charged-current overflow exchange", 4613);

  /* The with-number-policy solver validates the immutable metric at its own
   * boundary; every lapse/volume reject arm is a separate condition. */
  const double invalid_solver_metric_values[] = { NAN, 0.0, -1.0, INFINITY };
  for(size_t variant = 0; variant < sizeof(invalid_solver_metric_values)
                                          / sizeof(invalid_solver_metric_values[0]);
      ++variant) {
    ghl_metric_quantities bad_solver_metric;
    ghl_initialize_metric(
          1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &bad_solver_metric);
    bad_solver_metric.lapse = invalid_solver_metric_values[variant];
    ghl_m1_neutrino_diagnostics_initialize(&nd);
    require_error(
          ghl_m1_solve_neutrino_implicit_homogeneous_update_with_number_policy(
                m1_params, &nu_params, &bad_solver_metric, &prims, &rates, 1.0, 1.0, 0.0,
                &input, &output, &exchange, &solve_diagnostics, &nd),
          ghl_error_m1_invalid_metric, "invalid number-policy lapse", 4613);
    require_condition(
          memcmp(&output, &input, sizeof(output)) == 0 && nd.source_failures == 1,
          "number-policy lapse failure was not transactional", 4613);
    ghl_initialize_metric(
          1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &bad_solver_metric);
    bad_solver_metric.sqrt_detgamma = invalid_solver_metric_values[variant];
    ghl_m1_neutrino_diagnostics_initialize(&nd);
    require_error(
          ghl_m1_solve_neutrino_implicit_homogeneous_update_with_number_policy(
                m1_params, &nu_params, &bad_solver_metric, &prims, &rates, 1.0, 1.0, 0.0,
                &input, &output, &exchange, &solve_diagnostics, &nd),
          ghl_error_m1_invalid_metric, "invalid number-policy volume", 4613);
    require_condition(
          memcmp(&output, &input, sizeof(output)) == 0 && nd.source_failures == 1,
          "number-policy volume failure was not transactional", 4613);
  }
}

static void
test_overflow_safe_number_source_terms(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  prims.vU[0] = 0.9;
  prims.u0 = 1.0 / sqrt(1.0 - prims.vU[0] * prims.vU[0]);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_state state = { .N = 1.0, .E = 1.0, .F = { 0.99, 0.0, 0.0 } };
  ghl_m1_neutrino_current current;
  require_error(
        ghl_m1_neutrino_derive_current(
              m1_params, &nu_params, &metric, &prims, &state, &current),
        ghl_success, "overflow-safe number source current", 4614);
  state.N = current.Gamma_N;

  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nue, 1.0e308, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, &rates);
  require_error(
        ghl_m1_validate_neutrino_rates(&rates, NULL), ghl_success,
        "overflow-safe number source rates", 4614);

  ghl_m1_sources sources = { 0 };
  double number_source = NAN;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, &prims, &state, &rates, &sources,
              &number_source),
        ghl_success, "overflow-safe interaction number source", 4614);
  require_close(
        number_source, 0.0, 0.0, 0.0, "cancelled interaction number source", 4614);

  /* The projected charged-current source evaluates the same kappa*N/Gamma
   * form. endpoint=Gamma=2 makes the exact source zero although the raw
   * kappa*endpoint intermediate overflows. */
  double dL_rad_cc = -17.0;
  require_error(
        ghl_m1_neutrino_charged_current_lepton_delta(
              &rates, 1.0, 0.0, true, 2.0, 2.0, &dL_rad_cc),
        ghl_success, "overflow-safe projected charged-current source", 4615);
  require_close(
        dL_rad_cc, 0.0, 0.0, 0.0, "cancelled projected charged-current source", 4615);

  /* Keep the legacy underflow-to-zero result for physically valid tiny rates:
   * both the source absorption product and the BE timestep products round to
   * zero in binary64, and the old direct expressions accepted those results. */
  const double true_min = nextafter(0.0, 1.0);
  nu_params.N_floor = 0.0;
  state.N = true_min;
  make_rates(ghl_m1_neutrino_nue, true_min, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, &rates);
  require_error(
        ghl_m1_validate_neutrino_rates(&rates, NULL), ghl_success,
        "tiny valid number source rates", 4616);
  number_source = NAN;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &nu_params, &metric, &prims, &state, &rates, &sources,
              &number_source),
        ghl_success, "tiny absorption source underflow", 4616);
  require_condition(
        number_source == true_min,
        "tiny absorption source lost its legacy underflow result", 4616);
  dL_rad_cc = -17.0;
  require_error(
        ghl_m1_neutrino_charged_current_lepton_delta(
              &rates, 1.0, 0.0, true, true_min, 1.0, &dL_rad_cc),
        ghl_success, "tiny projected absorption underflow", 4617);
  require_condition(
        dL_rad_cc == true_min, "tiny projected source lost its legacy underflow result",
        4617);
  double tiny_endpoint = -1.0;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &rates, true_min, 1.0, 0.0, &tiny_endpoint),
        ghl_success, "tiny BE products underflow", 4618);
  require_condition(
        tiny_endpoint == 0.0, "tiny BE endpoint changed its legacy rounded value", 4618);

  /* Emission overflow with a finite absorption term takes the scaled
   * fallback through the emission-nonderminateness arm. */
  ghl_m1_neutrino_rates emission_overflow_rates;
  make_rates(
        ghl_m1_neutrino_nue, ldexp(1.0, -50), 0.0, 0.0, ldexp(1.0, 600), 1.0, 0.0, 0.0,
        &emission_overflow_rates);
  require_error(
        ghl_m1_validate_neutrino_rates(&emission_overflow_rates, NULL), ghl_success,
        "emission-overflow BE rates", 4619);
  double emission_overflow_endpoint = -1.0;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &emission_overflow_rates, ldexp(1.0, 600), 1.0, 1.0,
              &emission_overflow_endpoint),
        ghl_success, "emission-overflow BE endpoint", 4619);
  require_close(
        emission_overflow_endpoint, ldexp(1.0, 600), 2.0e-15, 0.0,
        "emission-overflow BE equilibrium", 4619);

  /* An absorption term that overflows after a finite product (tiny Gamma)
   * takes the term-nonderminateness arm; the scaled fallback publishes the
   * tiny representable equilibrium. */
  ghl_m1_neutrino_rates term_overflow_rates;
  make_rates(
        ghl_m1_neutrino_nue, ldexp(1.0, 500), 0.0, 0.0, ldexp(1.0, -500), 1.0, 0.0, 0.0,
        &term_overflow_rates);
  make_rates(
        ghl_m1_neutrino_nue, ldexp(1.0, 504), 0.0, 0.0, ldexp(1.0, -44), ldexp(1.0, 44),
        0.0, 0.0, &term_overflow_rates);
  require_error(
        ghl_m1_validate_neutrino_rates(&term_overflow_rates, NULL), ghl_success,
        "term-overflow BE rates", 4620);
  double term_overflow_endpoint = -1.0;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &nu_params, &term_overflow_rates, ldexp(1.0, 482), 2.0e-12, 1.0,
              &term_overflow_endpoint),
        ghl_success, "term-overflow BE endpoint", 4620);
  require_condition(
        term_overflow_endpoint > 0.0 && isfinite(term_overflow_endpoint)
              && term_overflow_endpoint < 1.0,
        "term-overflow BE published a nonrepresentable endpoint", 4620);
}

static void initialize_pair_boundary_inputs(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng,
      ghl_metric_quantities *restrict metric,
      ghl_primitive_quantities *restrict prims,
      ghl_m1_neutrino_parameters nu_params[2],
      ghl_m1_neutrino_rates rates[2],
      ghl_m1_neutrino_state state_input[2],
      ghl_m1_neutrino_state state_transport[2],
      ghl_m1_neutrino_state state_out[2],
      ghl_m1_neutrino_exchange exchange[2],
      ghl_m1_neutrino_source_diagnostics diagnostics[2],
      ghl_m1_neutrino_diagnostics neutrino_diagnostics[2]) {
  make_metric(rng, 1.0, metric);
  make_primitives(prims);
  make_neutrino_parameters(&nu_params[0]);
  make_neutrino_parameters(&nu_params[1]);
  make_rates(ghl_m1_neutrino_nue, 0.08, 0.12, 0.10, 1.0, 2.0, 0.01, 0.02, &rates[0]);
  make_rates(ghl_m1_neutrino_anue, 0.08, 0.12, 0.10, 1.0, 2.1, 0.01, 0.03, &rates[1]);
  for(int species = 0; species < 2; ++species) {
    make_state(rng, metric, &state_transport[species]);
    state_input[species] = state_transport[species];
    state_input[species].N *= 1.1;
    state_out[species] = (ghl_m1_neutrino_state){ .N = -1.0 - species,
                                                  .E = -2.0 - species,
                                                  .F = { -3.0 - species, -4.0 - species,
                                                         -5.0 - species } };
    exchange[species] = (ghl_m1_neutrino_exchange){
      .dN_rad_total = -1.0 - species,
      .dL_rad_cc = -2.0 - species,
      .dE_rad = -3.0 - species,
      .dF_rad = { -4.0 - species, -5.0 - species, -6.0 - species },
      .dTau_matter = -7.0 - species,
      .dS_matter = { -8.0 - species, -9.0 - species, -10.0 - species },
      .dYe_matter = -11.0 - species
    };
    diagnostics[species].path = ghl_m1_neutrino_source_path_thin_explicit;
    diagnostics[species].terminal_no_update = true;
    ghl_m1_neutrino_diagnostics_initialize(&neutrino_diagnostics[species]);
  }
  (void)m1_params;
}

static void test_pair_dispatcher_boundaries(
      const ghl_m1_parameters *restrict m1_params,
      m1_test_rng *restrict rng) {
  ghl_metric_quantities metric;
  ghl_primitive_quantities prims;
  ghl_m1_neutrino_parameters nu_params[2];
  ghl_m1_neutrino_rates rates[2];
  ghl_m1_neutrino_state state_input[2], state_transport[2], state_out[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics neutrino_diagnostics[2];
  initialize_pair_boundary_inputs(
        m1_params, rng, &metric, &prims, nu_params, rates, state_input, state_transport,
        state_out, exchange, diagnostics, neutrino_diagnostics);

  /* The public pair API has a single required-input boundary. Vary one
   * pointer at a time so each operand is independently exercised. */
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              NULL, nu_params, &metric, &prims, rates, state_input, state_transport, 0.1,
              3.0, state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair M1 parameters", 4300);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, NULL, &metric, &prims, rates, state_input, state_transport, 0.1,
              3.0, state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair neutrino parameters", 4301);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, NULL, &prims, rates, state_input, state_transport,
              0.1, 3.0, state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair metric", 4302);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, NULL, rates, state_input, state_transport,
              0.1, 3.0, state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair primitives", 4303);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, NULL, state_input, state_transport,
              0.1, 3.0, state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair rates", 4304);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, NULL, state_transport, 0.1,
              3.0, state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair input state", 4305);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, NULL, 0.1, 3.0,
              state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair transport state", 4306);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, state_transport,
              0.1, 3.0, NULL, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair state output", 4307);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, state_transport,
              0.1, 3.0, state_out, NULL, diagnostics, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair exchange output", 4308);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, state_transport,
              0.1, 3.0, state_out, exchange, NULL, neutrino_diagnostics),
        ghl_error_m1_null_pointer, "null pair source diagnostics", 4309);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, state_transport,
              0.1, 3.0, state_out, exchange, diagnostics, NULL),
        ghl_error_m1_null_pointer, "null pair neutrino diagnostics", 4310);

  /* Validation failures must publish both hard-failure status and a complete
   * no-update packet, not leave the caller's sentinels visible. */
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, state_transport,
              NAN, 3.0, state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite pair timestep", 4311);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(
                &state_out[species], &state_transport[species],
                sizeof(state_out[species]))
                == 0,
          "pair invalid timestep changed state", 4311);
    require_zero_exchange(&exchange[species], "pair invalid timestep exchange", 4311);
    require_condition(
          diagnostics[species].path == ghl_m1_neutrino_source_path_hard_failure,
          "pair invalid timestep status", 4311);
  }

  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, state_transport,
              0.1, NAN, state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_invalid_state, "nonfinite pair baryon normalization", 4312);
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state_input, state_transport,
              0.1, 0.0, state_out, exchange, diagnostics, neutrino_diagnostics),
        ghl_error_m1_invalid_state, "zero pair baryon normalization", 4312);
  ghl_m1_neutrino_rates bad_rates[2] = { rates[0], rates[1] };
  bad_rates[0].species = ghl_m1_neutrino_anue;
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, bad_rates, state_input,
              state_transport, 0.1, 3.0, state_out, exchange, diagnostics,
              neutrino_diagnostics),
        ghl_error_m1_microphysics_failure, "wrong pair species order", 4313);
  bad_rates[0] = rates[0];
  bad_rates[1] = rates[1];
  bad_rates[1].species = ghl_m1_neutrino_nue;
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, bad_rates, state_input,
              state_transport, 0.1, 3.0, state_out, exchange, diagnostics,
              neutrino_diagnostics),
        ghl_error_m1_microphysics_failure, "wrong second pair species", 4313);
  bad_rates[0] = rates[0];
  bad_rates[1] = rates[1];
  bad_rates[0].eta_N = -1.0;
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, bad_rates, state_input,
              state_transport, 0.1, 3.0, state_out, exchange, diagnostics,
              neutrino_diagnostics),
        ghl_error_m1_microphysics_failure, "invalid pair rates", 4314);
  bad_rates[0] = rates[0];
  bad_rates[1] = rates[1];
  bad_rates[1].eta_N_pair[0] += 0.01;
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, bad_rates, state_input,
              state_transport, 0.1, 3.0, state_out, exchange, diagnostics,
              neutrino_diagnostics),
        ghl_error_m1_microphysics_failure, "mismatched pair number rates", 4315);

  /* A pair-active packet with no pair energy emissivity is physically valid:
   * it exercises the effective-rate zero-energy branch while retaining a
   * shared number reaction. */
  ghl_m1_neutrino_rates zero_pair_energy_rates[2] = { rates[0], rates[1] };
  for(int species = 0; species < 2; ++species) {
    for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
      zero_pair_energy_rates[species].eta_E_pair[process] = 0.0;
    }
  }
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, zero_pair_energy_rates, state_input,
              state_transport, 0.05, 3.0, state_out, exchange, diagnostics,
              neutrino_diagnostics),
        ghl_success, "pair number-only source update", 4317);
  for(int species = 0; species < 2; ++species) {
    require_state_finite_and_admissible(
          m1_params, &nu_params[species], &metric, &state_out[species], 4317);
    require_exchange_contract(
          &state_transport[species], &state_out[species], &exchange[species], &metric,
          3.0, 4317);
  }
}

/* Arithmetic helpers are shared across several source policies. Test their
 * failure contracts directly, including preservation of caller output. */
static void test_arithmetic_failure_contracts(void) {
  volatile double values[] = { NAN, INFINITY, -1.0, 0.0, DBL_MIN, DBL_MAX };
  for(int i = 0; i < 8; ++i) {
    double kappa = 1.0, number = 1.0, gamma = 1.0, absorption = 7.0;
    if(i == 1) {
      kappa = values[0];
    }
    if(i == 2) {
      kappa = values[2];
    }
    if(i == 3) {
      number = values[0];
    }
    if(i == 4) {
      number = values[2];
    }
    if(i == 5) {
      gamma = values[0];
    }
    if(i == 6) {
      gamma = values[2];
    }
    if(i == 7) {
      gamma = values[3];
    }
    require_condition(
          !ghl_m1_neutrino_absorption_number(
                kappa, number, gamma, i == 0 ? NULL : &absorption)
                && absorption == 7.0,
          "absorption accepted invalid input or changed output", 4650 + i);
  }
  double absorption;
  require_condition(
        ghl_m1_neutrino_absorption_number(values[3], 1.0, 1.0, &absorption)
              && absorption == 0.0,
        "zero opacity must give zero absorption", 4660);
  require_condition(
        ghl_m1_neutrino_absorption_number(1.0, values[3], 1.0, &absorption)
              && absorption == 0.0,
        "zero occupancy must give zero absorption", 4661);
  require_condition(
        ghl_m1_neutrino_absorption_number(values[4], values[4], 1.0, &absorption)
              && absorption == 0.0,
        "absorption underflow must preserve IEEE zero", 4662);
  require_condition(
        !ghl_m1_neutrino_absorption_number(values[5], values[5], values[4], &absorption),
        "nonrepresentable absorption must fail", 4663);
}

/* The scaled-positive helpers reject every NULL/invalid operand arm
 * transactionally. Load volatile values into plain locals and keep the pointer
 * operands volatile so reject paths cannot be constant-folded without
 * discarding any pointee qualifiers. */
static void test_scaled_helper_boundary_arms(void) {
  ghl_m1_scaled_positive a = { .mantissa = 0.5, .exponent = 1 };
  ghl_m1_scaled_positive b = { .mantissa = 0.25, .exponent = 3 };
  ghl_m1_scaled_positive zero = { .mantissa = 0.0, .exponent = 0 };
  ghl_m1_scaled_positive out = { .mantissa = -1.0, .exponent = -7 };
  ghl_m1_scaled_positive *volatile null_operand = NULL;
  ghl_m1_scaled_positive *volatile a_operand = &a;
  ghl_m1_scaled_positive *volatile b_operand = &b;
  ghl_m1_scaled_positive *volatile out_operand = &out;
  double quotient_sentinel = -3.0;
  double *volatile quotient_out = &quotient_sentinel;
  double *volatile null_quotient = NULL;
  const volatile double good_value = 2.0;
  const volatile double negative_value = -2.0;
  const volatile double nan_value = NAN;
  const volatile double tiny_pair[2] = { DBL_MIN, DBL_MIN };
  const volatile double unit_value = 1.0;
  const volatile double huge_pair[2] = { DBL_MAX, DBL_MAX };
  const volatile double min_value = DBL_MIN;
  const volatile double denominator_value = 4.0;
  /* The helpers receive ordinary const double storage after volatile reads. */
  const double good_operand = good_value;
  const double negative_operand = negative_value;
  const double nan_operand = nan_value;
  const double tiny_operands[2] = { tiny_pair[0], tiny_pair[1] };
  const double unit_operand = unit_value;
  const double huge_operands[2] = { huge_pair[0], huge_pair[1] };
  const double min_operand = min_value;
  const double denominator_operand = denominator_value;
  const double *volatile good_values = &good_operand;
  const double *volatile negative_values = &negative_operand;
  const double *volatile nan_values = &nan_operand;
  const double *volatile tiny_values = tiny_operands;
  const double *volatile unit_values = &unit_operand;
  const double *volatile huge_values = huge_operands;
  const double *volatile min_values = &min_operand;
  const double *volatile denominator_values = &denominator_operand;
  const double *volatile null_values = NULL;

  /* multiply: NULL operands leave the output untouched. */
  out = (ghl_m1_scaled_positive){ .mantissa = -1.0, .exponent = -7 };
  require_condition(
        !ghl_m1_scaled_positive_multiply(null_operand, b_operand, out_operand)
              && out.mantissa == -1.0 && out.exponent == -7,
        "scaled multiply accepted a NULL left operand", 4700);
  require_condition(
        !ghl_m1_scaled_positive_multiply(a_operand, null_operand, out_operand)
              && out.mantissa == -1.0 && out.exponent == -7,
        "scaled multiply accepted a NULL right operand", 4701);
  require_condition(
        !ghl_m1_scaled_positive_multiply(a_operand, b_operand, null_operand),
        "scaled multiply accepted a NULL output", 4702);
  require_condition(
        ghl_m1_scaled_positive_multiply(a_operand, b_operand, out_operand)
              && out.mantissa > 0.0,
        "scaled multiply rejected a productive pair", 4703);

  /* add: NULL operands, NULL output, and the zero-operand shortcuts. */
  out = (ghl_m1_scaled_positive){ .mantissa = -1.0, .exponent = -7 };
  require_condition(
        !ghl_m1_scaled_positive_add(a_operand, null_operand, out_operand)
              && out.mantissa == -1.0 && out.exponent == -7,
        "scaled add accepted a NULL right operand", 4710);
  require_condition(
        !ghl_m1_scaled_positive_add(null_operand, b_operand, out_operand)
              && out.mantissa == -1.0 && out.exponent == -7,
        "scaled add accepted a NULL left operand", 4711);
  require_condition(
        !ghl_m1_scaled_positive_add(a_operand, b_operand, null_operand),
        "scaled add accepted a NULL sum output", 4712);
  require_condition(
        ghl_m1_scaled_positive_add(a_operand, &zero, out_operand) && out.mantissa == 0.5
              && out.exponent == 1,
        "scaled add lost the zero-right shortcut", 4714);

  /* divide: NULL operands, a NULL output, a zero denominator, and a zero
   * numerator. */
  quotient_sentinel = -3.0;
  require_condition(
        !ghl_m1_scaled_positive_divide(null_operand, b_operand, quotient_out)
              && quotient_sentinel == -3.0,
        "scaled divide accepted a NULL numerator", 4720);
  require_condition(
        !ghl_m1_scaled_positive_divide(a_operand, null_operand, quotient_out)
              && quotient_sentinel == -3.0,
        "scaled divide accepted a NULL denominator", 4721);
  require_condition(
        !ghl_m1_scaled_positive_divide(a_operand, b_operand, null_quotient),
        "scaled divide accepted a NULL quotient output", 4722);
  require_condition(
        !ghl_m1_scaled_positive_divide(a_operand, &zero, quotient_out)
              && quotient_sentinel == -3.0,
        "scaled divide accepted a zero denominator", 4723);
  require_condition(
        ghl_m1_scaled_positive_divide(&zero, a_operand, quotient_out)
              && quotient_sentinel == 0.0,
        "scaled divide lost its zero-numerator shortcut", 4724);
  quotient_sentinel = -3.0;
  require_condition(
        ghl_m1_scaled_positive_divide(a_operand, b_operand, quotient_out)
              && quotient_sentinel > 0.0 && isfinite(quotient_sentinel),
        "scaled divide rejected a productive quotient", 4725);

  /* compare: the exponent, mantissa-less-than, greater, and equal arms. */
  const ghl_m1_scaled_positive same_mantissa_left = { .mantissa = 0.5, .exponent = 2 };
  const ghl_m1_scaled_positive same_mantissa_right = { .mantissa = 0.5, .exponent = 3 };
  require_condition(
        ghl_m1_scaled_positive_compare(&same_mantissa_left, &same_mantissa_right) < 0,
        "scaled compare missed the exponent order", 4730);
  const ghl_m1_scaled_positive mantissa_left = { .mantissa = 0.25, .exponent = 3 };
  const ghl_m1_scaled_positive mantissa_right = { .mantissa = 0.5, .exponent = 3 };
  require_condition(
        ghl_m1_scaled_positive_compare(&mantissa_left, &mantissa_right) < 0,
        "scaled compare missed the smaller mantissa", 4731);
  require_condition(
        ghl_m1_scaled_positive_compare(&mantissa_right, &mantissa_left) > 0,
        "scaled compare missed the larger mantissa", 4732);
  require_condition(
        ghl_m1_scaled_positive_compare(&mantissa_right, &mantissa_right) == 0,
        "scaled compare missed the equal pair", 4733);
  require_condition(
        ghl_m1_scaled_positive_compare(&zero, &zero) == 0,
        "scaled compare missed the equal zeros", 4734);
  require_condition(
        ghl_m1_scaled_positive_compare(&zero, &mantissa_right) < 0,
        "scaled compare missed the zero-left order", 4735);
  require_condition(
        ghl_m1_scaled_positive_compare(&mantissa_right, &zero) > 0,
        "scaled compare missed the zero-right order", 4736);

  /* product: NULL values, a nonpositive count, a NULL output, and invalid
   * factors. */
  require_condition(
        !ghl_m1_scaled_positive_product(null_values, 1, out_operand),
        "scaled product accepted NULL values", 4740);
  require_condition(
        !ghl_m1_scaled_positive_product(good_values, 0, out_operand),
        "scaled product accepted an empty factor list", 4741);
  require_condition(
        !ghl_m1_scaled_positive_product(good_values, 2, null_operand),
        "scaled product accepted a NULL output", 4742);
  require_condition(
        !ghl_m1_scaled_positive_product(negative_values, 1, out_operand),
        "scaled product accepted a negative factor", 4743);
  require_condition(
        !ghl_m1_scaled_positive_product(nan_values, 1, out_operand),
        "scaled product accepted a nonfinite factor", 4744);
  require_condition(
        ghl_m1_scaled_positive_product(good_values, 1, out_operand)
              && out.mantissa > 0.0,
        "scaled product rejected a productive factor", 4745);

  /* sqrt: NULL operands and the exact-zero shortcut. */
  require_condition(
        !ghl_m1_scaled_positive_sqrt(null_operand, out_operand),
        "scaled sqrt accepted a NULL value", 4750);
  require_condition(
        !ghl_m1_scaled_positive_sqrt(a_operand, null_operand),
        "scaled sqrt accepted a NULL root output", 4751);
  require_condition(
        ghl_m1_scaled_positive_sqrt(&zero, out_operand) && out.mantissa == 0.0
              && out.exponent == 0,
        "scaled sqrt lost its exact-zero shortcut", 4752);

  /* ratio_of_products: NULL quotient, a failing numerator product, a failing
   * denominator product, an underflowed zero quotient, and an unrepresentable
   * quotient. */
  double ratio_sentinel = -3.0;
  double *volatile ratio_out = &ratio_sentinel;
  double *volatile null_ratio = NULL;
  require_condition(
        !ghl_m1_neutrino_scaled_ratio_of_products(
              good_values, 1, denominator_values, 1, null_ratio),
        "scaled ratio accepted a NULL quotient", 4760);
  require_condition(
        !ghl_m1_neutrino_scaled_ratio_of_products(
              negative_values, 1, denominator_values, 1, ratio_out)
              && ratio_sentinel == -3.0,
        "scaled ratio accepted a failing numerator product", 4761);
  require_condition(
        !ghl_m1_neutrino_scaled_ratio_of_products(
              good_values, 1, nan_values, 1, ratio_out)
              && ratio_sentinel == -3.0,
        "scaled ratio accepted a failing denominator product", 4762);
  /* A true quotient below binary64's representable range must fail rather
   * than publish zero: (DBL_MIN * DBL_MIN) / 1 = 2^-2043. */
  require_condition(
        !ghl_m1_neutrino_scaled_ratio_of_products(
              tiny_values, 2, unit_values, 1, ratio_out)
              && ratio_sentinel == -3.0,
        "scaled ratio published an underflowed zero quotient", 4763);
  /* (DBL_MAX * DBL_MAX) / DBL_MIN = 2^3070 is beyond binary64's range. */
  require_condition(
        !ghl_m1_neutrino_scaled_ratio_of_products(
              huge_values, 2, min_values, 1, ratio_out)
              && ratio_sentinel == -3.0,
        "scaled ratio published an unrepresentable quotient", 4764);
  require_condition(
        ghl_m1_neutrino_scaled_ratio_of_products(
              good_values, 1, denominator_values, 1, ratio_out)
              && ratio_sentinel == 0.5,
        "scaled ratio rejected a productive quotient", 4765);
  /* The policy-selection helper rejects invalid product factors directly;
   * load the volatile operand storage into plain locals so the reject arms
   * stay executable without discarding the volatile qualifier. */
  {
    const volatile double threshold_value = 1.0;
    const volatile double negative_value = -1.0;
    const volatile double nan_value = NAN;
    const double threshold_operands[1] = { threshold_value };
    const double negative_operands[1] = { negative_value };
    const double nan_operands[1] = { nan_value };
    const double *volatile threshold_values = threshold_operands;
    const double *volatile negative_values = negative_operands;
    const double *volatile nan_values = nan_operands;
    const double *volatile null_values = NULL;
    require_condition(
          !ghl_m1_neutrino_scaled_product_meets_threshold(
                null_values, 1, threshold_values, 1, true),
          "scaled threshold accepted NULL factors", 4770);
    require_condition(
          !ghl_m1_neutrino_scaled_product_meets_threshold(
                negative_values, 1, threshold_values, 1, true),
          "scaled threshold accepted a negative factor", 4771);
    require_condition(
          !ghl_m1_neutrino_scaled_product_meets_threshold(
                nan_values, 1, threshold_values, 1, true),
          "scaled threshold accepted a nonfinite factor", 4772);
    require_condition(
          !ghl_m1_neutrino_scaled_product_meets_threshold(
                threshold_values, 1, null_values, 1, true),
          "scaled threshold accepted NULL threshold factors", 4773);
    require_condition(
          !ghl_m1_neutrino_scaled_product_meets_threshold(
                threshold_values, 1, negative_values, 1, true),
          "scaled threshold accepted a negative threshold factor", 4774);
    require_condition(
          ghl_m1_neutrino_scaled_product_meets_threshold(
                threshold_values, 1, threshold_values, 1, true),
          "scaled threshold rejected an equal inclusive pair", 4775);
    require_condition(
          !ghl_m1_neutrino_scaled_product_meets_threshold(
                threshold_values, 1, threshold_values, 1, false),
          "scaled threshold accepted an equal exclusive pair", 4776);
    require_condition(
          !ghl_m1_neutrino_scaled_product_meets_threshold(
                threshold_values, 0, threshold_values, 1, true),
          "scaled threshold accepted an empty product", 4777);

    /* The scaled-positive conversion rejects each invalid operand arm. */
    ghl_m1_scaled_positive scaled;
    ghl_m1_scaled_positive *volatile scaled_out = &scaled;
    ghl_m1_scaled_positive *volatile null_scaled = NULL;
    require_condition(
          !ghl_m1_scaled_positive_from_double(1.0, null_scaled),
          "scaled conversion accepted a NULL output", 4780);
    require_condition(
          !ghl_m1_scaled_positive_from_double(NAN, scaled_out),
          "scaled conversion accepted a nonfinite value", 4781);
    require_condition(
          !ghl_m1_scaled_positive_from_double(-1.0, scaled_out),
          "scaled conversion accepted a negative value", 4782);
    require_condition(
          ghl_m1_scaled_positive_from_double(2.0, scaled_out) && scaled.mantissa > 0.0,
          "scaled conversion rejected a productive value", 4783);
  }

  /* A zero numerator is a valid exact zero quotient, not an underflow. */
  const volatile double zero_value = 0.0;
  const double zero_operands[1] = { zero_value };
  const double *volatile zero_values = zero_operands;
  ratio_sentinel = -3.0;
  require_condition(
        ghl_m1_neutrino_scaled_ratio_of_products(
              zero_values, 1, unit_values, 1, ratio_out)
              && ratio_sentinel == 0.0,
        "scaled ratio rejected an exact zero quotient", 4766);
}

static void test_number_endpoint_validation_contracts(void) {
  ghl_m1_neutrino_parameters params;
  make_neutrino_parameters(&params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nux, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_state state = { .E = 1.0, .N = 1.0 };
  ghl_m1_neutrino_current current = { .J = 1.0, .Gamma_N = 1.0 };
  for(int i = 0; i < 5; ++i) {
    double result = 7.0;
    require_error(
          ghl_m1_neutrino_update_endpoint_number_with_policy(
                i == 0 ? NULL : &params, i == 1 ? NULL : &rates, 1.0, 1.0, -1.0,
                i == 2 ? NULL : &state, i == 3 ? NULL : &current,
                i == 4 ? NULL : &result, NULL),
          ghl_error_m1_null_pointer, "number policy null input", 4670 + i);
    require_condition(
          result == 7.0, "number policy changed output on null input", 4670 + i);
  }
  for(int i = 0; i < 5; ++i) {
    double result = 7.0;
    require_error(
          ghl_m1_neutrino_update_endpoint_number_with_policy(
                &params, &rates,
                i == 0   ? NAN
                : i == 1 ? -1.0
                         : 1.0,
                i == 2   ? NAN
                : i == 3 ? -1.0
                         : 1.0,
                i == 4 ? NAN : -1.0, &state, &current, &result, NULL),
          ghl_error_m1_invalid_state, "number policy invalid time", 4680 + i);
    require_condition(
          result == 7.0, "number policy changed output on invalid time", 4680 + i);
  }
  double result;
  require_error(
        ghl_m1_neutrino_update_endpoint_number_with_policy(
              &params, &rates, 1.0, 1.0, -1.0, &state, &current, &result, NULL),
        ghl_success, "number policy optional diagnostic", 4690);
  require_close(result, 1.0, 0.0, 0.0, "equilibrium BE endpoint", 4690);
  rates.mean_energy = 0.0;
  require_error(
        ghl_m1_neutrino_update_endpoint_number_with_policy(
              &params, &rates, 1.0, 1.0, 0.0, &state, &current, &result, NULL),
        ghl_success, "zero mean energy endpoint", 4691);
  require_close(result, 0.0, 0.0, 0.0, "zero mean energy endpoint", 4691);
  const double invalid_floors[] = { NAN, -1.0 };
  for(int i = 0; i < 2; ++i) {
    params.Gamma_N_floor = invalid_floors[i];
    require_error(
          ghl_m1_update_neutrino_number_backward_euler(
                &params, &rates, 1.0, 1.0, 1.0, &result),
          ghl_error_m1_invalid_state, "invalid Gamma floor", 4692 + i);
  }
  params.Gamma_N_floor = 0.0;
  rates.kappa_a_N = 1.0;
  rates.n_eq = rates.J_eq = rates.eta_N = DBL_MAX;
  rates.mean_energy = 1.0;
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &params, &rates, 2.0, 4.0, 0.0, &result),
        ghl_error_m1_invalid_state, "unrepresentable emission endpoint", 4694);
  require_error(
        ghl_m1_update_neutrino_number_backward_euler(
              &params, &rates, 2.0, 4.0, -1.0, &result),
        ghl_error_m1_invalid_state, "negative scaled number endpoint", 4695);
}

static void test_implicit_context_validation(const ghl_m1_parameters *m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nux, 0.0, 0.1, 0.0, 1.0, 1.0, 0.0, 0.0, &rates);
  const double U[4] = { 1.0, 0.0, 0.0, 0.0 };
  const double residual_zero[4] = { 0.0 };
  double residual[4], jacobian[4][4];
  for(int i = 0; i < 9; ++i) {
    const ghl_m1_neutrino_implicit_context context
          = { .m1_params = i == 1 ? NULL : m1_params,
              .metric = i == 2 ? NULL : &metric,
              .prims_frozen = i == 3 ? NULL : &prims,
              .rates = i == 4 ? NULL : &rates };
    if(i != 8) {
      require_error(
            ghl_m1_neutrino_compute_implicit_residual_validated(
                  i == 0 ? NULL : &context, 1.0, i == 5 ? NULL : U, i == 6 ? NULL : U,
                  NULL, i == 7 ? NULL : residual),
            ghl_error_m1_null_pointer, "residual context null input", 4700 + i);
    }
    require_error(
          ghl_m1_neutrino_compute_implicit_jacobian_validated(
                i == 0 ? NULL : &context, 1.0, i == 5 ? NULL : U, i == 6 ? NULL : U,
                i == 7 ? NULL : residual_zero, i == 8 ? NULL : jacobian),
          ghl_error_m1_null_pointer, "Jacobian context null input", 4710 + i);
  }
  const ghl_m1_neutrino_implicit_context context
        = { m1_params, &metric, &prims, &rates };
  for(int i = 0; i < 2; ++i) {
    require_error(
          ghl_m1_neutrino_compute_implicit_residual_validated(
                &context, i == 0 ? NAN : -1.0, U, U, NULL, residual),
          ghl_error_m1_invalid_state, "residual invalid time", 4720 + i);
  }
  double bad_base[4] = { 1.0, NAN, 0.0, 0.0 };
  require_error(
        ghl_m1_neutrino_compute_implicit_residual_validated(
              &context, 1.0, bad_base, U, NULL, residual),
        ghl_error_m1_invalid_state, "residual nonfinite base", 4722);
  for(int i = 0; i < 6; ++i) {
    ghl_m1_sources sources = { .S_E = 7.0 };
    require_error(
          ghl_m1_neutrino_compute_EF_interaction_sources_validated(
                i == 0 ? NULL : m1_params, i == 1 ? NULL : &metric,
                i == 2 ? NULL : &prims, i == 3 ? NULL : &(ghl_m1_rad_state){ .E = 1.0 },
                i == 4 ? NULL : &rates, NULL, i == 5 ? NULL : &sources),
          ghl_error_m1_null_pointer, "EF source null input", 4730 + i);
    require_close(
          sources.S_E, 7.0, 0.0, 0.0, "EF null input output preservation", 4730 + i);
  }
  ghl_m1_neutrino_parameters params;
  make_neutrino_parameters(&params);
  ghl_m1_neutrino_state state = { .E = 1.0, .N = 1.0 };
  ghl_m1_neutrino_current current = { .J = 1.0, .Gamma_N = 1.0 };
  for(int i = 0; i < 5; ++i) {
    ghl_m1_neutrino_diagnostics diagnostics = { .mean_energy_diag = 7.0 };
    ghl_m1_neutrino_populate_mean_energy_diagnostics(
          i == 0 ? NULL : &state, i == 1 ? NULL : &current, i == 2 ? NULL : &params,
          i == 3 ? NULL : &rates, i == 4 ? NULL : &diagnostics);
    require_close(
          diagnostics.mean_energy_diag, 7.0, 0.0, 0.0,
          "mean diagnostic null input preservation", 4740 + i);
  }
  for(int i = 0; i < 4; ++i) {
    ghl_m1_neutrino_state changed_state = state;
    ghl_m1_neutrino_current changed_current = current;
    ghl_m1_neutrino_rates changed_rates = rates;
    if(i == 0) {
      changed_state.N = NAN;
    }
    if(i == 1) {
      changed_current.Gamma_N = NAN;
    }
    if(i == 2) {
      changed_rates.mean_energy = 0.0;
    }
    if(i == 3) {
      changed_rates.n_eq = DBL_MIN;
      changed_rates.J_eq = DBL_MAX;
    }
    ghl_m1_neutrino_diagnostics diagnostics = { 0 };
    ghl_m1_neutrino_populate_mean_energy_diagnostics(
          &changed_state, &changed_current, &params, &changed_rates, &diagnostics);
    require_condition(
          i < 2 ? diagnostics.mean_energy_diag_invalid
                : !diagnostics.Jeq_over_neq_consistent,
          "invalid mean diagnostic input was accepted", 4750 + i);
  }
}

static void test_source_validation_and_overflow(const ghl_m1_parameters *m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters params;
  make_neutrino_parameters(&params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nux, 1.0, 0.1, 0.0, 1.0, 1.0, 0.0, 0.0, &rates);
  const ghl_m1_neutrino_state state = { .E = 1.0, .N = 1.0 };
  const ghl_m1_rad_state rad = { .E = 1.0 };
  ghl_m1_closure closure;
  require_error(
        ghl_m1_compute_closure_with_primitives(
              m1_params, &metric, &prims, &rad, &closure),
        ghl_success, "source boundary reference closure", 4760);
  for(int i = 0; i < 4; ++i) {
    ghl_m1_neutrino_parameters bad_params = params;
    ghl_m1_neutrino_state bad_state = state;
    if(i == 0) {
      bad_params.N_floor = -1.0;
    }
    if(i == 1) {
      bad_params.N_floor = NAN;
    }
    if(i == 2) {
      bad_state.N = NAN;
    }
    if(i == 3) {
      bad_state.N = 0.0;
    }
    ghl_m1_sources sources = { .S_E = 7.0 };
    double number_source = 7.0;
    require_error(
          ghl_m1_compute_neutrino_interaction_sources(
                m1_params, &bad_params, &metric, &prims, &bad_state, &rates, &sources,
                &number_source),
          ghl_error_m1_invalid_state, "source invalid number state", 4761 + i);
    require_error(
          ghl_m1_compute_neutrino_interaction_sources_from_closure(
                m1_params, &bad_params, &metric, &prims, &bad_state, &closure, &rates,
                &sources, &number_source),
          ghl_error_m1_invalid_state, "closure source invalid number state", 4761 + i);
    require_condition(
          sources.S_E == 7.0 && number_source == 7.0, "invalid source changed output",
          4761 + i);
  }
  for(int i = 0; i < 2; ++i) {
    ghl_m1_neutrino_rates changed_rates = rates;
    ghl_primitive_quantities changed_prims = prims;
    if(i == 0) {
      changed_rates.kappa_a_E = -1.0;
    }
    if(i == 1) {
      changed_prims.vU[0] = NAN;
    }
    ghl_m1_sources sources = { .S_E = 7.0 };
    require_error(
          ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
                m1_params, &metric, &changed_prims, &rad, &changed_rates, NULL,
                &sources),
          i == 0 ? ghl_error_m1_microphysics_failure : ghl_error_m1_invalid_state,
          "EF source invalid rates or primitive", 4770 + i);
    require_close(sources.S_E, 7.0, 0.0, 0.0, "EF failure preservation", 4770 + i);
  }
  for(int i = 0; i < 5; ++i) {
    ghl_m1_neutrino_state changed_state = state;
    if(i == 4) {
      changed_state.E = NAN;
    }
    ghl_m1_neutrino_state result;
    ghl_m1_neutrino_exchange exchange;
    ghl_m1_implicit_solve_diagnostics solve;
    ghl_m1_neutrino_diagnostics diagnostics = { 0 };
    require_error(
          ghl_m1_solve_neutrino_implicit_homogeneous_update_with_number_policy(
                m1_params, &params, &metric, &prims, &rates, i == 0 ? NAN : 1.0,
                i == 1   ? NAN
                : i == 2 ? 0.0
                         : 1.0,
                i == 3 ? NAN : -1.0, &changed_state, &result, &exchange, &solve,
                &diagnostics),
          ghl_error_m1_invalid_state, "homogeneous invalid input", 4780 + i);
    require_condition(
          diagnostics.source_failures == 1,
          "homogeneous invalid input did not record failure", 4780 + i);
    require_zero_exchange(&exchange, "homogeneous invalid input exchange", 4780 + i);
  }
  /* Sources can be nonrepresentable even though each frozen rate is valid. */
  ghl_m1_neutrino_rates large_rates;
  make_rates(
        ghl_m1_neutrino_nux, DBL_MAX, DBL_MAX, 0.0, 1.0, 1.0, 0.0, 0.0, &large_rates);
  for(int i = 0; i < 2; ++i) {
    ghl_m1_neutrino_state large_state = state;
    if(i == 0) {
      large_state.E = 2.0;
    }
    if(i == 1) {
      large_state.N = 2.0;
    }
    ghl_m1_sources sources = { .S_E = 7.0 };
    double number_source = 7.0;
    require_error(
          ghl_m1_compute_neutrino_interaction_sources(
                m1_params, &params, &metric, &prims, &large_state, &large_rates,
                &sources, &number_source),
          ghl_error_m1_invalid_state, "nonrepresentable interaction source", 4790 + i);
    require_condition(
          sources.S_E == 7.0 && number_source == 7.0, "overflowed source changed output",
          4790 + i);
  }
}

static void test_source_residual_range_errors(const ghl_m1_parameters *m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters params;
  make_neutrino_parameters(&params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nux, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, &rates);
  const double U[4] = { 1.0, 0.0, 0.0, 0.0 };
  const double zero[4] = { 0.0 };
  double residual[4], jacobian[4][4];
  ghl_m1_parameters no_fd = *m1_params;
  no_fd.fd_epsilon_abs = no_fd.fd_epsilon_rel = 0.0;
  require_error(
        ghl_m1_neutrino_compute_implicit_jacobian_with_base(
              &no_fd, &metric, &prims, &rates, 0.0, U, U, zero, jacobian),
        ghl_error_m1_invalid_state, "zero finite difference delta", 4800);
  for(int i = 0; i < 3; ++i) {
    ghl_metric_quantities changed_metric = metric;
    ghl_m1_neutrino_rates changed_rates = rates;
    if(i == 0) {
      changed_metric.lapse = 2.0;
      changed_metric.lapseinv = 0.5;
      changed_metric.lapseinv2 = 0.25;
    }
    if(i == 1) {
      ghl_initialize_metric(
            1.0, 0.0, 0.0, 0.0, 1.0e100, 0.0, 0.0, 1.0e100, 0.0, 1.0e100,
            &changed_metric);
    }
    if(i == 2) {
      changed_rates.kappa_a_E = changed_rates.kappa_tr = DBL_MAX;
      changed_rates.J_eq = changed_rates.n_eq = 0.0;
    }
    double changed_U[4] = { changed_metric.sqrt_detgamma, 0.0, 0.0, 0.0 };
    require_error(
          ghl_m1_neutrino_compute_implicit_residual_with_base(
                m1_params, &changed_metric, &prims, &changed_rates,
                i < 2 ? DBL_MAX : 2.0, changed_U, changed_U, residual),
          ghl_error_m1_invalid_state, "nonrepresentable residual arithmetic", 4801 + i);
  }
  ghl_metric_quantities bad_metric = metric;
  bad_metric.sqrt_detgamma = 0.0;
  const ghl_m1_neutrino_state state = { .N = 1.0, .E = 1.0 };
  require_error(
        ghl_m1_neutrino_compute_implicit_residual(
              m1_params, &params, &bad_metric, &prims, &rates, &state, 0.0, U, residual),
        ghl_error_m1_invalid_metric, "zero residual volume", 4804);
  ghl_m1_neutrino_current current = { .J = 1.0, .Gamma_N = 1.0 };
  double number;
  require_error(
        ghl_m1_neutrino_update_endpoint_number_with_policy(
              &params, &rates, 1.0, 1.0, 1.0, &state, &current, &number, NULL),
        ghl_success, "inactive nonnegative number policy threshold", 4805);
  require_close(number, state.N, 0.0, 0.0, "inactive number policy identity", 4805);
  double output = 7.0;
  volatile double tiny = DBL_MIN;
  require_condition(
        !ghl_m1_neutrino_absorption_number(4.0, 1.0, tiny, &output),
        "overflowed quotient accepted", 4810);
}

static void
test_pair_effective_rate_failure_contracts(const ghl_m1_parameters *m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters params[2];
  make_neutrino_parameters(&params[0]);
  make_neutrino_parameters(&params[1]);
  const ghl_m1_neutrino_state input[2]
        = { { .N = 1.0, .E = 1.0 }, { .N = 1.0, .E = 1.0 } };
  for(int i = 0; i < 4; ++i) {
    ghl_m1_neutrino_rates rates[2];
    make_rates(
          ghl_m1_neutrino_nue, 0.0, 0.0, 0.0, i == 1 ? 0.0 : 1.0, 1.0, 0.0, 0.0,
          &rates[0]);
    make_rates(
          ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, i == 0 ? 0.0 : 1.0, 1.0, 0.0, 0.0,
          &rates[1]);
    if(i == 0) {
      rates[0].eta_E_pair[0] = 1.0;
    }
    if(i == 1) {
      rates[1].eta_E_pair[0] = 1.0;
    }
    if(i == 2) {
      rates[0].eta_E_pair[0] = DBL_MAX;
      rates[0].eta_E_pair[1] = DBL_MAX;
    }
    if(i == 3) {
      for(int species = 0; species < 2; ++species) {
        rates[species].eta_E_pair[0] = 0.3 * DBL_MAX;
        rates[species].eta_E_pair[1] = 0.3 * DBL_MAX;
        rates[species].eta_E_pair[2] = 0.3 * DBL_MAX;
      }
    }
    ghl_m1_neutrino_state pair_input[2] = { input[0], input[1] };
    if(i == 3) {
      pair_input[0].N = 2.0;
      pair_input[1].N = 2.0;
    }
    ghl_m1_neutrino_state output[2];
    ghl_m1_neutrino_exchange exchange[2];
    ghl_m1_neutrino_source_diagnostics diagnostics[2];
    ghl_m1_neutrino_diagnostics neutrino_diagnostics[2] = { { 0 } };
    require_error(
          ghl_m1_solve_neutrino_pair_source_update(
                m1_params, params, &metric, &prims, rates, pair_input, pair_input, 1.0,
                1.0, output, exchange, diagnostics, neutrino_diagnostics),
          ghl_error_m1_invalid_state, "invalid effective pair rates", 4820 + i);
    for(int species = 0; species < 2; ++species) {
      require_condition(
            memcmp(&pair_input[species], &output[species], sizeof(input[species])) == 0,
            "effective pair failure changed state", 4820 + i);
      require_zero_exchange(
            &exchange[species], "effective pair failure exchange", 4820 + i);
    }
  }
}

static void test_pair_number_endpoint_failures(const ghl_m1_parameters *m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params[2];
  ghl_m1_neutrino_rates rates[2];
  for(int species = 0; species < 2; ++species) {
    make_neutrino_parameters(&nu_params[species]);
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0,
          DBL_MAX, 1.0, 1.0 / 3.0, 0.0, &rates[species]);
    nu_params[species].N_floor = 0.0;
  }
  const ghl_m1_neutrino_state transport[2]
        = { { .N = DBL_MAX, .E = 1.0 }, { .N = 0.0, .E = 1.0 } };
  ghl_m1_neutrino_state output[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics nd[2] = { { 0 }, { 0 } };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, transport, transport,
              DBL_MAX, 1.0, output, exchange, diagnostics, nd),
        ghl_error_m1_implicit_terminal_fallback, "pair number endpoint overflow", 4830);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&output[species], &transport[species], sizeof(output[species])) == 0
                && diagnostics[species].path
                         == ghl_m1_neutrino_source_path_terminal_no_update
                && diagnostics[species].terminal_no_update
                && nd[species].source_terminal_fallbacks == 1,
          "pair number overflow was not transactional", 4830);
    require_zero_exchange(&exchange[species], "pair number overflow exchange", 4830);
  }

  /* With Gamma_N=1.25, equal finite equilibrium populations produce an
   * unrepresentable smaller root too, not merely y + (N_e-N_a). */
  prims.vU[0] = 0.6;
  prims.u0 = 1.25;
  const ghl_m1_neutrino_state moving[2]
        = { { .N = 1.0, .E = 1.75, .F = { 1.25, 0.0, 0.0 } },
            { .N = 1.0, .E = 1.75, .F = { 1.25, 0.0, 0.0 } } };
  for(int species = 0; species < 2; ++species) {
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0,
          DBL_MAX, 1.0, DBL_MAX / 4.0, 0.0, &rates[species]);
    nd[species] = (ghl_m1_neutrino_diagnostics){ 0 };
    ghl_m1_neutrino_current current;
    require_error(
          ghl_m1_neutrino_derive_current(
                m1_params, &nu_params[species], &metric, &prims, &moving[species],
                &current),
          ghl_success, "moving pair overflow current", 4833);
    require_close(
          current.Gamma_N, 1.25, 64.0 * DBL_EPSILON, 0.0,
          "moving pair number normalization", 4833);
  }
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, moving, moving, DBL_MAX, 1.0,
              output, exchange, diagnostics, nd),
        ghl_error_m1_implicit_terminal_fallback, "both pair number endpoints overflow",
        4833);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&output[species], &moving[species], sizeof(output[species])) == 0
                && diagnostics[species].terminal_no_update
                && nd[species].source_terminal_fallbacks == 1,
          "equal pair overflow changed state or diagnostics", 4833);
    require_zero_exchange(&exchange[species], "equal pair overflow exchange", 4833);
  }
  make_primitives(&prims);

  for(int species = 0; species < 2; ++species) {
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 0.1,
          1.0, 100.0 / 3.0, 0.0, &rates[species]);
    nu_params[species].N_floor = 0.5;
  }
  const ghl_m1_neutrino_state above_floor[2]
        = { { .N = 1.0, .E = 1.0 }, { .N = 1.0, .E = 1.0 } };
  nd[0] = (ghl_m1_neutrino_diagnostics){ 0 };
  nd[1] = (ghl_m1_neutrino_diagnostics){ 0 };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, above_floor, above_floor,
              1.0, 1.0, output, exchange, diagnostics, nd),
        ghl_error_m1_implicit_terminal_fallback, "pair endpoint below number floor",
        4831);
  for(int species = 0; species < 2; ++species) {
    require_condition(
          memcmp(&output[species], &above_floor[species], sizeof(output[species])) == 0
                && diagnostics[species].path
                         == ghl_m1_neutrino_source_path_terminal_no_update
                && nd[species].source_terminal_fallbacks == 1,
          "pair number-floor failure was not transactional", 4831);
    require_zero_exchange(&exchange[species], "pair number-floor exchange", 4831);
  }
}

static void test_pair_scaled_opacity_fallback(const ghl_m1_parameters *m1_params) {
  /* Zero pair number emissivity leaves both number endpoints unchanged.
   * For nue, eta_E_pair/J_eq = 2^-67; multiplying the anti-nue occupancy
   * 2^-1074 rounds to zero before division by its n_eq = DBL_MIN.
   * The scaled quotient remains representable at 2^-119. */
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params[2];
  ghl_m1_neutrino_rates rates[2];
  for(int species = 0; species < 2; ++species) {
    make_neutrino_parameters(&nu_params[species]);
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0,
          species == 0 ? 1.0 : DBL_MIN, 1.0, 0.0, species == 0 ? 0x1p-67 : 0.0,
          &rates[species]);
    nu_params[species].N_floor = 0.0;
  }
  const double partner_n = 0x1p-1074;
  const ghl_m1_neutrino_state state[2]
        = { { .N = 1.0, .E = 1.0 }, { .N = partner_n, .E = 1.0 } };
  ghl_m1_neutrino_state output[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics nd[2] = { { 0 } };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, state, state, 0.25, 1.0,
              output, exchange, diagnostics, nd),
        ghl_success, "scaled pair effective-opacity fallback", 4846);
  require_condition(
        output[0].N == state[0].N && output[1].N == state[1].N,
        "zero pair emissivity changed number endpoints", 4846);
  require_condition(
        !diagnostics[0].terminal_no_update && !diagnostics[1].terminal_no_update,
        "scaled pair opacity fallback published terminal no-update", 4846);
  for(int species = 0; species < 2; ++species) {
    require_state_finite_and_admissible(
          m1_params, &nu_params[species], &metric, &output[species], 4846);
    require_condition(
          nd[species].source_failures == 0,
          "scaled pair opacity fallback did not converge cleanly", 4846);
  }
}

static void test_pair_zero_partner_occupancy(const ghl_m1_parameters *m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_parameters nu_params[2];
  ghl_m1_neutrino_rates rates[2];
  for(int species = 0; species < 2; ++species) {
    make_neutrino_parameters(&nu_params[species]);
    make_rates(
          species == 0 ? ghl_m1_neutrino_nue : ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 1.0,
          1.0, 0.0, 0.1, &rates[species]);
  }
  nu_params[1].N_floor = 0.0;
  const ghl_m1_neutrino_state input[2]
        = { { .N = 1.0, .E = 1.0 }, { .N = 0.0, .E = 1.0 } };
  ghl_m1_neutrino_state output[2];
  ghl_m1_neutrino_exchange exchange[2];
  ghl_m1_neutrino_source_diagnostics diagnostics[2];
  ghl_m1_neutrino_diagnostics nd[2] = { { 0 }, { 0 } };
  require_error(
        ghl_m1_solve_neutrino_pair_source_update(
              m1_params, nu_params, &metric, &prims, rates, input, input, 0.1, 1.0,
              output, exchange, diagnostics, nd),
        ghl_success, "pair emission with empty partner", 4832);
  require_condition(
        output[0].N == input[0].N && output[1].N == 0.0 && output[0].E > input[0].E
              && output[1].E >= m1_params->E_floor,
        "empty pair partner did not retain its number and emission", 4832);
  for(int species = 0; species < 2; ++species) {
    require_state_finite_and_admissible(
          m1_params, &nu_params[species], &metric, &output[species], 4832);
  }
}

/* Preserve outputs when only a momentum source (rather than the energy
 * source) overflows. */
static void test_partial_source_line_boundaries(const ghl_m1_parameters *m1_params) {
  ghl_m1_neutrino_parameters params;
  make_neutrino_parameters(&params);
  ghl_m1_neutrino_rates rates;
  double number = -7.0;
  const ghl_m1_neutrino_state state = { .E = 1.0, .N = 1.0 };

  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  make_rates(ghl_m1_neutrino_nux, 0.0, 0.0, DBL_MAX, 1.0, 1.0, 0.0, 0.0, &rates);
  const ghl_m1_neutrino_state flux_state
        = { .E = 2.0, .F = { 2.0, 0.0, 0.0 }, .N = 1.0 };
  const ghl_m1_sources sentinel = { .S_E = -7.0, .S = { -8.0, -9.0, -10.0 } };
  ghl_m1_sources sources = sentinel;
  number = -11.0;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &params, &metric, &prims, &flux_state, &rates, &sources,
              &number),
        ghl_error_m1_invalid_state, "momentum-only source overflow", 4902);
  require_condition(
        memcmp(&sources, &sentinel, sizeof(sources)) == 0 && number == -11.0,
        "momentum overflow changed outputs", 4902);

  /* The explicit E/F endpoint can fall below the number-current floor;
   * the shortcut must withdraw that otherwise admissible partial update. */
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  params.J_floor = 0.9;
  make_rates(ghl_m1_neutrino_nux, 0.5, 0.5, 0.0, 0.1, 1.0, 0.0, 0.0, &rates);
  ghl_m1_neutrino_source_options options = branched_options();
  ghl_m1_neutrino_state output;
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics nd = { 0 };
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &params, &metric, &prims, &rates, &state, &state, 1.0,
              1.0, &output, &exchange, &diagnostics, &nd),
        ghl_error_m1_invalid_state, "thin endpoint below current floor", 4904);
  require_condition(
        memcmp(&output, &state, sizeof(state)) == 0 && nd.source_failures == 1,
        "failed shortcut published partial endpoint", 4904);
  require_zero_exchange(&exchange, "failed shortcut exchange", 4904);

  prims.vU[0] = NAN;
  require_error(
        ghl_m1_compute_neutrino_interaction_sources(
              m1_params, &params, &metric, &prims, &flux_state, &rates, &sources,
              &number),
        ghl_error_m1_invalid_state, "interaction closure rejects invalid velocity",
        4903);
  require_condition(
        memcmp(&sources, &sentinel, sizeof(sources)) == 0 && number == -11.0,
        "closure failure changed outputs", 4903);

  /* The scaled closure remains representable at this energy, while the
   * subsequent unscaled comoving energy overflows for a moving fluid. */
  prims.vU[0] = 0.8;
  const double large_energy = 0.2 * DBL_MAX;
  const ghl_m1_rad_state large_rad
        = { .E = large_energy, .F = { 0.1 * large_energy, 0.0, 0.0 } };
  ghl_m1_closure large_closure;
  require_error(
        ghl_m1_compute_closure_with_primitives(
              m1_params, &metric, &prims, &large_rad, &large_closure),
        ghl_success, "large-energy source closure", 4906);
  ghl_m1_comoving large_comoving;
  require_error(
        ghl_m1_compute_comoving_moments(
              m1_params, &metric, &prims, &large_rad, &large_closure, &large_comoving),
        ghl_error_m1_invalid_state, "large-energy comoving source moments", 4906);
  sources = sentinel;
  require_error(
        ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
              m1_params, &metric, &prims, &large_rad, &rates, NULL, &sources),
        ghl_error_m1_invalid_state, "large-energy diagnostic source moments", 4906);
  require_condition(
        memcmp(&sources, &sentinel, sizeof(sources)) == 0,
        "comoving source failure changed outputs", 4906);
}

static void
test_branched_nonrepresentable_projected_number(const ghl_m1_parameters *m1_params) {
  ghl_m1_neutrino_parameters params;
  make_neutrino_parameters(&params);
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  const ghl_m1_neutrino_state transport = { .E = 1.0, .N = 1.0 };
  const double true_min = nextafter(0.0, 1.0);

  for(int thick = 0; thick < 2; ++thick) {
    ghl_m1_neutrino_rates rates;
    make_rates(
          ghl_m1_neutrino_nux, 0.0, thick ? 1.0 : 0.0, thick ? 1.0 : 0.0, 1.0, true_min,
          0.0, 0.0, &rates);
    require_error(
          ghl_m1_validate_neutrino_rates(&rates, NULL), ghl_success,
          "projected-number fixture rates", 4910 + thick);
    ghl_m1_neutrino_source_options options = branched_options();
    options.thermalized_number_threshold = 0.0;
    ghl_m1_neutrino_state output = { 0 };
    ghl_m1_neutrino_exchange exchange;
    ghl_m1_neutrino_source_diagnostics diagnostics;
    ghl_m1_neutrino_diagnostics nd = { 0 };
    require_error(
          ghl_m1_solve_neutrino_source_update(
                &options, m1_params, &params, &metric, &prims, &rates, &transport,
                &transport, 2.0, 1.0, &output, &exchange, &diagnostics, &nd),
          ghl_error_m1_invalid_state, "branched projected number range", 4910 + thick);
    require_condition(
          memcmp(&output, &transport, sizeof(output)) == 0
                && diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
                && nd.source_failures == 1,
          "failed projected number published a state", 4910 + thick);
    require_zero_exchange(&exchange, "failed projected number exchange", 4910 + thick);
  }
}

static void test_thin_endpoint_energy_overflow(const ghl_m1_parameters *m1_params) {
  ghl_m1_neutrino_parameters params;
  make_neutrino_parameters(&params);
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  prims.vU[0] = 0.99;
  ghl_m1_neutrino_rates rates;
  make_rates(
        ghl_m1_neutrino_nux, 0.0, 1.0e-300, 0.0, 1.0, 0.5 * DBL_MAX, 0.0, 0.0, &rates);
  require_error(
        ghl_m1_validate_neutrino_rates(&rates, NULL), ghl_success,
        "thin overflow fixture rates", 4912);
  const ghl_m1_neutrino_state transport = { .E = 1.0, .N = 1.0 };
  ghl_m1_neutrino_source_options options = branched_options();
  ghl_m1_neutrino_state output = { 0 };
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics nd = { 0 };
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &params, &metric, &prims, &rates, &transport,
              &transport, 5.0e299, 1.0, &output, &exchange, &diagnostics, &nd),
        ghl_error_m1_invalid_state, "thin endpoint energy overflow", 4912);
  require_condition(
        memcmp(&output, &transport, sizeof(output)) == 0
              && diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && nd.source_failures == 1,
        "thin endpoint overflow published a state", 4912);
  require_zero_exchange(&exchange, "thin endpoint overflow exchange", 4912);
}

static void
test_branched_large_coordinate_number_flux(const ghl_m1_parameters *m1_params) {
  ghl_m1_neutrino_parameters params;
  make_neutrino_parameters(&params);
  ghl_metric_quantities metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0, 1.0e-100, 0.0, 0.0, 1.0e100, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  prims.vU[0] = 5.0e49;
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nux, 0.0, 0.0, 0.0, 1.0, 1.0e-270, 0.0, 0.0, &rates);
  ghl_m1_neutrino_source_options options = branched_options();
  options.thermalized_number_threshold = 0.0;
  const ghl_m1_neutrino_state transport = { .E = 1.0, .N = 1.0 };
  ghl_m1_neutrino_state output = { 0 };
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics nd = { 0 };
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &params, &metric, &prims, &rates, &transport,
              &transport, 1.0, 1.0, &output, &exchange, &diagnostics, &nd),
        ghl_error_m1_invalid_state, "large-coordinate number flux", 4915);
  require_condition(
        memcmp(&output, &transport, sizeof(output)) == 0
              && diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && nd.source_failures == 1,
        "large-coordinate number flux published a state", 4915);
  require_zero_exchange(&exchange, "large-coordinate number flux exchange", 4915);
}

static void
test_stiff_moments_realizability_tolerance(const ghl_m1_parameters *m1_params) {
  ghl_m1_neutrino_parameters params;
  make_neutrino_parameters(&params);
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nux, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, &rates);
  const double cone = sqrt(1.0 - m1_params->epsilon_c);
  const ghl_m1_neutrino_state transport
        = { .E = 1.0, .F = { cone * (1.0 + 96.0 * DBL_EPSILON) }, .N = 1.0 };
  const ghl_m1_rad_state rad = ghl_m1_neutrino_project_rad_state(&transport);
  require_error(
        ghl_m1_validate_realizability_state(m1_params, &metric, &rad, 128.0, NULL),
        ghl_success, "dispatcher realizability tolerance", 4916);
  require_error(
        ghl_m1_validate_realizability_state(m1_params, &metric, &rad, 64.0, NULL),
        ghl_error_m1_invalid_state, "moments realizability tolerance", 4916);
  ghl_m1_closure closure;
  require_error(
        ghl_m1_compute_closure_with_primitives(
              m1_params, &metric, &prims, &rad, &closure),
        ghl_success, "stiff tolerance fixture closure", 4916);
  ghl_m1_neutrino_source_options options = branched_options();
  ghl_m1_neutrino_state output = { 0 };
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics nd = { 0 };
  require_error(
        ghl_m1_solve_neutrino_source_update(
              &options, m1_params, &params, &metric, &prims, &rates, &transport,
              &transport, 2.0, 1.0, &output, &exchange, &diagnostics, &nd),
        ghl_error_m1_invalid_state, "stiff predictor moments tolerance", 4916);
  require_condition(
        memcmp(&output, &transport, sizeof(output)) == 0
              && diagnostics.path == ghl_m1_neutrino_source_path_hard_failure
              && nd.source_failures == 1,
        "stiff moments failure published a state", 4916);
  require_zero_exchange(&exchange, "stiff moments failure exchange", 4916);
}

/* Geometry and interaction sources can each be representable even when
 * their required lapse-volume scaling is not. */
static void test_rhs_lapse_volume_range(const ghl_m1_parameters *m1_params) {
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);
  ghl_m1_neutrino_rates rates;
  make_rates(ghl_m1_neutrino_nux, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, &rates);
  ghl_primitive_quantities prims;
  make_primitives(&prims);
  const ghl_m1_neutrino_state state = { .E = 1.0, .N = 1.0 };
  const ghl_metric_quantities derivs = { 0 };
  const ghl_extrinsic_curvature curv = { 0 };
  for(int underflow = 0; underflow < 2; ++underflow) {
    const double g = underflow ? 1e-100 : cbrt(DBL_MAX);
    ghl_metric_quantities metric = { 0 };
    metric.lapse = underflow ? 1e-200 : sqrt(DBL_MAX);
    metric.detgamma = underflow ? g * g * g : DBL_MAX;
    /* The determinant consistency contract permits this two-ulp rounding
     * difference; the square is checked without overflowing double. */
    metric.sqrt_detgamma
          = underflow ? sqrt(metric.detgamma)
                      : nextafter(nextafter(sqrt(DBL_MAX), INFINITY), INFINITY);
    ghl_m1_closure closure = { 0 };
    for(int i = 0; i < 3; ++i) {
      metric.gammaDD[i][i] = g;
      metric.gammaUU[i][i] = 1.0 / g;
      closure.P[i][i] = 1.0 / (3.0 * g);
    }
    ghl_m1_sources individual;
    double number;
    require_error(
          ghl_m1_compute_neutrino_geometry_sources(
                m1_params, &metric, &derivs, &derivs, &derivs, &curv, &state, &closure,
                &individual),
          ghl_success, "range fixture geometry source", 4905 + underflow);
    require_error(
          ghl_m1_compute_neutrino_interaction_sources(
                m1_params, &nu_params, &metric, &prims, &state, &rates, &individual,
                &number),
          ghl_success, "range fixture interaction source", 4905 + underflow);
    const double alpha_volume = metric.lapse * metric.sqrt_detgamma;
    require_condition(
          underflow ? alpha_volume == 0.0 : isinf(alpha_volume),
          "lapse-volume fixture did not exceed arithmetic range", 4905 + underflow);
    double energy = 7.0, momentum[3] = { 8.0, 9.0, 10.0 }, count = 11.0;
    require_error(
          ghl_m1_compute_neutrino_explicit_rhs_sources(
                m1_params, &nu_params, &metric, &derivs, &derivs, &derivs, &curv, &prims,
                &state, &closure, true, &rates, &energy, momentum, &count),
          ghl_error_m1_invalid_metric, "unrepresentable lapse-volume RHS",
          4905 + underflow);
    require_condition(
          energy == 7.0 && momentum[0] == 8.0 && momentum[1] == 9.0
                && momentum[2] == 10.0 && count == 11.0,
          "lapse-volume rejection modified output", 4905 + underflow);
  }
}

static void
test_transverse_moving_stiff_branches(const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_m1_neutrino_parameters nu_params;
  make_neutrino_parameters(&nu_params);

  /* Rotate the reported v=0.5 scattering and weak-absorption cases onto y
   * and z. These directions must also select the physical implicit update. */
  for(int direction = 1; direction < 3; ++direction) {
    ghl_primitive_quantities prims;
    make_primitives(&prims);
    prims.vU[direction] = 0.5;
    for(int thick = 0; thick < 2; ++thick) {
      const int case_index = 5000 + 2 * direction + thick;
      ghl_m1_neutrino_rates rates;
      make_rates(
            ghl_m1_neutrino_nux, 0.0, thick ? 0.001 : 0.0, 1.0e6, 1.0, 1.0, 0.0, 0.0,
            &rates);
      const ghl_m1_neutrino_state input = { .N = 1.0, .E = 1.0 };
      ghl_m1_neutrino_state transport = input;
      transport.F[direction] = 0.1;
      ghl_m1_neutrino_source_options options = branched_options();
      options.thick_equilibrium_threshold = thick ? 1.0 : 0.0;
      options.scattering_threshold = 10.0;
      ghl_m1_neutrino_state output;
      ghl_m1_neutrino_exchange exchange;
      ghl_m1_neutrino_source_diagnostics diagnostics;
      ghl_m1_neutrino_diagnostics nd = { 0 };
      require_error(
            ghl_m1_solve_neutrino_source_update(
                  &options, m1_params, &nu_params, &metric, &prims, &rates, &input,
                  &transport, 1.0, 1.0, &output, &exchange, &diagnostics, &nd),
            ghl_success, "transverse moving stiff source", case_index);
      require_condition(
            diagnostics.path
                        == (thick ? ghl_m1_neutrino_source_path_thick_equilibrium
                                  : ghl_m1_neutrino_source_path_scattering_dominated)
                  && !diagnostics.closure_fallback_used
                  && !diagnostics.terminal_no_update
                  && diagnostics.implicit.newton_iterations > 0
                  && diagnostics.implicit.fallback_substeps == 1
                  && diagnostics.implicit.residual_scaled_norm <= 1.0
                  && nd.source_converged == 1 && nd.N_floor_repairs == 0
                  && nd.EF_repairs == 0,
            "transverse stiff branch lost physical solve diagnostics", case_index);
      require_state_finite_and_admissible(
            m1_params, &nu_params, &metric, &output, case_index);
      require_exchange_contract(
            &transport, &output, &exchange, &metric, 1.0, case_index);
      require_close(output.N, transport.N, 0.0, 0.0, "scattering number", case_index);
      require_close(
            exchange.dL_rad_cc, 0.0, 0.0, 0.0, "nu_x lepton exchange", case_index);
      require_close(exchange.dYe_matter, 0.0, 0.0, 0.0, "nu_x Ye exchange", case_index);

      ghl_m1_sources sources;
      double number_source;
      require_error(
            ghl_m1_compute_neutrino_interaction_sources(
                  m1_params, &nu_params, &metric, &prims, &output, &rates, &sources,
                  &number_source),
            ghl_success, "transverse endpoint interaction source", case_index);
      const double delta_contraction = output.E - 0.5 * output.F[direction]
                                       - (transport.E - 0.5 * transport.F[direction]);
      const double source_contraction = sources.S_E - 0.5 * sources.S[direction];
      /* With lapse=dt=1 and one successful step, the endpoint source must
       * equal the full contracted increment. Use the initialized solver's
       * mixed tolerances rather than introducing a new numerical threshold. */
      require_close(
            delta_contraction, source_contraction, m1_params->newton_tolerance,
            m1_params->newton_absolute_tolerance, "transverse contracted BE equation",
            case_index);
      if(!thick) {
        require_close(
              delta_contraction, 0.0, m1_params->newton_tolerance,
              m1_params->newton_absolute_tolerance, "elastic scattering conservation",
              case_index);
      }
    }
  }
}

#if defined(__linux__)
#include "m1_source_update_fault_injection.h"
#endif

int main(void) {
  test_arithmetic_failure_contracts();
  test_scaled_helper_boundary_arms();
  test_number_endpoint_validation_contracts();
  ghl_m1_parameters m1_params;
  ghl_error_codes_t error = ghl_m1_initialize_with_newton_tolerances(
        1.0e-10, 1.0e-10, 1.0e-6, 1.0e-8, 1.0e-10, 40, 1.0e-10, 1.0e-12, &m1_params);
  if(error != ghl_success) {
    ghl_error("ghl_m1_initialize_with_newton_tolerances returned %d\n", (int)error);
  }

  test_partial_source_line_boundaries(&m1_params);
  test_branched_nonrepresentable_projected_number(&m1_params);
  test_thin_endpoint_energy_overflow(&m1_params);
  test_branched_large_coordinate_number_flux(&m1_params);
  test_stiff_moments_realizability_tolerance(&m1_params);
  test_transverse_moving_stiff_branches(&m1_params);
  test_source_validation_and_overflow(&m1_params);
  test_rhs_lapse_volume_range(&m1_params);
  test_implicit_context_validation(&m1_params);
  test_source_residual_range_errors(&m1_params);
  test_pair_effective_rate_failure_contracts(&m1_params);
  test_pair_number_endpoint_failures(&m1_params);
  test_pair_scaled_opacity_fallback(&m1_params);
  test_pair_zero_partner_occupancy(&m1_params);
  m1_test_rng rng = { .state = UINT64_C(0x4d315f534f555243) };
  test_source_regimes(&m1_params, &rng);
  test_source_base_and_lapse_scaling(&m1_params, &rng);
  test_scaled_ratio_underflow_consumers();
  test_thermalized_number_projection(&m1_params);
  test_stiff_branch_arithmetic_boundaries(&m1_params);
  test_reachable_scaled_ratio_failure(&m1_params);
  test_branched_general_thermalized_number_projection(&m1_params);
  /* Coverage-only public API checks must not perturb the established corpus
   * consumed by the pre-existing randomized tests below. */
  m1_test_rng public_source_rng = rng;
  test_public_source_api_boundaries(&m1_params, &public_source_rng);
  m1_test_rng dispatcher_boundary_rng = rng;
  test_source_dispatcher_boundaries(&m1_params, &dispatcher_boundary_rng);
  m1_test_rng selector_boundary_rng = rng;
  test_source_compatibility_and_selector_failures(&m1_params, &selector_boundary_rng);
  test_pair_source_conservation(&m1_params, &rng);
  test_pair_effective_opacity_overflow(&m1_params);
  test_pair_source_independent_oracle(&m1_params);
  test_pair_final_mean_energy_bounds(&m1_params);
  test_pair_state_derivation_failures(&m1_params);
  test_pair_newton_only_closure_fallback(&m1_params);
  m1_test_rng pair_inactive_rng = rng;
  test_pair_inactive_delegation(&m1_params, &pair_inactive_rng);
  m1_test_rng pair_publication_rng = rng;
  test_pair_public_transaction_publication(&m1_params, &pair_publication_rng);
  test_repairs_and_transactional_failures(&m1_params, &rng);
  test_stiff_charged_current_conservation(&m1_params);
  test_number_floor_charged_current_accounting(&m1_params);
  test_lepton_policies(&m1_params, &rng);
  test_retry_and_terminal_no_update(&m1_params, &rng);
  test_implicit_predictor_recovery(&m1_params);
  test_direct_neutrino_validation_boundaries(&m1_params, &rng);
  test_implicit_kernel_boundaries(&m1_params, &rng);
  test_endpoint_failure_transactions(&m1_params);
  test_pair_moving_closure_diagnostics(&m1_params);
  test_pair_projected_backtracking(&m1_params);
  test_late_closure_fallback_transaction(&m1_params);
  test_pair_roundtrip_energy_floor(&m1_params);
  test_number_endpoint_and_policy_overflow(&m1_params);
  test_overflow_safe_number_source_terms(&m1_params);
  test_projected_charged_current_overflow(&m1_params);
  test_pair_dispatcher_boundaries(&m1_params, &rng);

#if defined(__linux__)
  test_pair_final_exchange_failures(&m1_params);
  test_source_update_fault_injections(&m1_params);
  test_neutrino_repair_delegate_failure(&m1_params);
#endif

  ghl_info("M1 neutrino source-update randomized/property tests passed\n");
  return 0;
}
