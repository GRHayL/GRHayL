#ifndef M1_SOURCE_UPDATE_FAULT_INJECTION_H
#define M1_SOURCE_UPDATE_FAULT_INJECTION_H

#include "ghl.h"
#include "ghl_m1.h"
#include "../GRHayL/Radiation/ghl_m1_utils.h"
#include "../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit.h"

#include <float.h>
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * The library's validated public inputs cannot force every defensive return
 * in the source update. Interpose selected helper results in this executable
 * so the actual shared-library code must reject them transactionally.
 */

typedef enum {
  fault_none,
  fault_late_closure,
  fault_second_repair,
  fault_repair_delegate,
  fault_hn,
  fault_ruu,
  fault_candidate
} fault_stage;

static fault_stage stage;
static bool route_thick;
static int closure_calls, moments_calls, repair_calls, realizability_calls, adm_calls;

static void fail_case(const char *reason, const fault_stage selected, const bool thick) {
  fprintf(stderr, "source-update fault %d (thick=%d): %s\n", selected, thick, reason);
  exit(EXIT_FAILURE);
}

ghl_error_codes_t ghl_m1_compute_closure_with_primitives(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      ghl_m1_closure *restrict closure) {
  ++closure_calls;
  if(stage == fault_late_closure && closure_calls == (route_thick ? 2 : 3)) {
    return ghl_error_m1_invalid_state;
  }
  if(stage == fault_hn || stage == fault_ruu || stage == fault_candidate) {
    /* The injected moments have no physical closure tensor. */
    *closure = (ghl_m1_closure){ 0 };
    return ghl_success;
  }
  typedef ghl_error_codes_t (*real_fn)(
        const ghl_m1_parameters *, const ghl_metric_quantities *,
        const ghl_primitive_quantities *, const ghl_m1_rad_state *, ghl_m1_closure *);
  real_fn real = (real_fn)dlsym(RTLD_NEXT, "ghl_m1_compute_closure_with_primitives");
  if(real == NULL) abort();
  return real(m1_params, metric, prims, rad_state, closure);
}

ghl_error_codes_t ghl_m1_compute_comoving_moments_with_velocity(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving,
      double V_con[3], double V_cov[3], double *restrict W_out) {
  ++moments_calls;
  if(stage == fault_hn || stage == fault_ruu || stage == fault_candidate) {
    /* All fields are finite; the predictor must detect the later overflow. */
    *comoving = (ghl_m1_comoving){ 0 };
    comoving->J = stage == fault_ruu ? 3.0e155 : 1.0;
    if(stage == fault_hn) {
      comoving->HD[0] = 100.0;
    }
    for(int i = 0; i < 3; ++i) {
      V_con[i] = 0.0;
      V_cov[i] = 0.0;
    }
    if(stage == fault_hn) {
      V_cov[0] = DBL_MAX;
    }
    *W_out = 1.0;
    return ghl_success;
  }
  typedef ghl_error_codes_t (*real_fn)(
        const ghl_m1_parameters *, const ghl_metric_quantities *,
        const ghl_primitive_quantities *, const ghl_m1_rad_state *,
        const ghl_m1_closure *, ghl_m1_comoving *, double[3], double[3], double *);
  real_fn real = (real_fn)dlsym(RTLD_NEXT, "ghl_m1_compute_comoving_moments_with_velocity");
  if(real == NULL) abort();
  return real(m1_params, metric, prims, rad_state, closure, comoving, V_con, V_cov, W_out);
}

ghl_error_codes_t ghl_m1_repair_neutrino_state(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      ghl_m1_neutrino_state *restrict state,
      ghl_m1_neutrino_diagnostics *restrict diagnostics) {
  ++repair_calls;
  if(stage == fault_second_repair && repair_calls == 2) {
    return ghl_error_m1_invalid_state;
  }
  typedef ghl_error_codes_t (*real_fn)(
        const ghl_m1_parameters *, const ghl_m1_neutrino_parameters *,
        const ghl_metric_quantities *, ghl_m1_neutrino_state *,
        ghl_m1_neutrino_diagnostics *);
  real_fn real = (real_fn)dlsym(RTLD_NEXT, "ghl_m1_repair_neutrino_state");
  if(real == NULL) abort();
  return real(m1_params, nu_params, metric, state, diagnostics);
}

ghl_error_codes_t ghl_m1_realizability_repair(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      ghl_m1_rad_state *restrict rad_state) {
  ++realizability_calls;
  if(stage == fault_repair_delegate) {
    /* A failed delegate may have changed its private candidate already. */
    rad_state->E = 3.0;
    rad_state->F[0] = 2.0;
    return ghl_error_m1_invalid_state;
  }
  typedef ghl_error_codes_t (*real_fn)(
        const ghl_m1_parameters *, const ghl_metric_quantities *, ghl_m1_rad_state *);
  real_fn real = (real_fn)dlsym(RTLD_NEXT, "ghl_m1_realizability_repair");
  if(real == NULL) abort();
  return real(m1_params, metric, rad_state);
}

void ghl_compute_ADM_auxiliaries(
      const ghl_metric_quantities *restrict metric,
      ghl_ADM_aux_quantities *restrict metric_aux) {
  ++adm_calls;
  if(stage == fault_candidate) {
    /* Model a finite but inconsistent ADM inverse from the geometry helper. */
    *metric_aux = (ghl_ADM_aux_quantities){ 0 };
    metric_aux->g4DD[0][0] = -4.0;
    metric_aux->g4UU[0][0] = 1.0e154;
    for(int i = 1; i < 4; ++i) {
      metric_aux->g4DD[i][i] = 1.0;
      metric_aux->g4UU[i][i] = 1.0;
    }
    return;
  }
  typedef void (*real_fn)(const ghl_metric_quantities *, ghl_ADM_aux_quantities *);
  real_fn real = (real_fn)dlsym(RTLD_NEXT, "ghl_compute_ADM_auxiliaries");
  if(real == NULL) abort();
  real(metric, metric_aux);
}

static void run_case(
      const ghl_m1_parameters *restrict m1_params,
      const fault_stage selected_stage,
      const bool thick) {
  stage = selected_stage;
  route_thick = thick;
  closure_calls = moments_calls = repair_calls = adm_calls = 0;
  ghl_metric_quantities metric;
  if(stage == fault_ruu) {
    ghl_initialize_metric(
          1.0, 0.0, 0.0, 0.0, 1.0e-154, 0.0, 0.0, 1.0, 0.0, 1.0e100, &metric);
  }
  else {
    ghl_initialize_metric(
          stage == fault_candidate ? 2.0 : 1.0, 0.0, 0.0, 0.0,
          1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  }
  ghl_primitive_quantities prims = { 0 };
  ghl_m1_neutrino_parameters params = {
    .N_floor = 1.0e-12, .J_floor = 1.0e-14, .Gamma_N_floor = 1.0e-12,
    .terminal_fallback_policy = ghl_m1_neutrino_terminal_fallback_no_update_all
  };
  ghl_m1_neutrino_rates rates = {
    .species = ghl_m1_neutrino_nux,
    .mean_energy = 1.0,
    .n_eq = stage == fault_candidate ? 0.0 : 1.0,
    .J_eq = stage == fault_candidate ? 0.0 : 1.0,
    .kappa_a_E = thick ? 1.0 : 0.0,
    .kappa_s = thick ? 1.0 : 0.0,
    .kappa_tr = thick ? 2.0 : 0.0,
    .eta_E = thick && stage != fault_candidate ? 1.0 : 0.0
  };
  ghl_m1_neutrino_source_options options = {
    .policy = ghl_m1_neutrino_source_branched_compatibility,
    .thick_equilibrium_threshold = 0.5,
    .scattering_threshold = 0.5,
    .thermalized_number_threshold = -1.0,
    .allow_closure_fallback = true,
    .ye_policy = ghl_m1_neutrino_ye_from_charged_current
  };
  const ghl_m1_neutrino_state input = { .E = 1.0, .N = 1.0 };
  ghl_m1_neutrino_state output = { 0 };
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics nd = { 0 };
  const ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
        &options, m1_params, &params, &metric, &prims, &rates, &input, &input,
        2.0, 1.0, &output, &exchange, &diagnostics, &nd);
  if(error != ghl_error_m1_invalid_state
     || diagnostics.path != ghl_m1_neutrino_source_path_hard_failure
     || diagnostics.terminal_no_update || nd.source_failures != 1) {
    fail_case("error or diagnostic contract failed", stage, thick);
  }
  if(output.E != input.E || output.N != input.N) {
    fail_case("partial state was published", stage, thick);
  }
  for(int i = 0; i < 3; ++i) {
    if(output.F[i] != input.F[i] || exchange.dF_rad[i] != 0.0
       || exchange.dS_matter[i] != 0.0) {
      fail_case("partial flux or momentum exchange was published", stage, thick);
    }
  }
  if(exchange.dN_rad_total != 0.0 || exchange.dL_rad_cc != 0.0
     || exchange.dE_rad != 0.0 || exchange.dTau_matter != 0.0
     || exchange.dYe_matter != 0.0) {
    fail_case("partial exchange was published", stage, thick);
  }
  if((stage == fault_late_closure && closure_calls != (thick ? 2 : 3))
     || (stage == fault_second_repair && repair_calls != 2)
     || ((stage == fault_hn || stage == fault_ruu || stage == fault_candidate)
         && moments_calls != 1)
     || (stage == fault_hn && adm_calls != 0)
     || ((stage == fault_ruu || stage == fault_candidate) && adm_calls != 1)) {
    fail_case("injection did not reach its intended stage", stage, thick);
  }
}

static void test_source_update_fault_injections(
      const ghl_m1_parameters *restrict m1_params) {
  run_case(m1_params, fault_late_closure, false);
  run_case(m1_params, fault_late_closure, true);
  run_case(m1_params, fault_second_repair, false);
  run_case(m1_params, fault_second_repair, true);
  run_case(m1_params, fault_hn, true);
  run_case(m1_params, fault_ruu, true);
  run_case(m1_params, fault_candidate, true);
  stage = fault_none;
}

static void test_neutrino_repair_delegate_failure(
      const ghl_m1_parameters *restrict m1_params) {
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  const ghl_m1_neutrino_parameters nu_params = { .N_floor = 2.0 };
  ghl_m1_neutrino_state state = { .N = 1.0, .E = 1.0 };
  ghl_m1_neutrino_diagnostics diagnostics = { .N_floor_repairs = 7, .repair_dN = 3.0 };
  unsigned char state_before[sizeof(state)], diagnostics_before[sizeof(diagnostics)];
  memcpy(state_before, &state, sizeof(state));
  memcpy(diagnostics_before, &diagnostics, sizeof(diagnostics));

  stage = fault_repair_delegate;
  realizability_calls = 0;
  const ghl_error_codes_t error = ghl_m1_repair_neutrino_state(
        m1_params, &nu_params, &metric, &state, &diagnostics);
  if(error != ghl_error_m1_invalid_state || realizability_calls != 1) {
    fail_case("repair delegate failure was not propagated", stage, false);
  }
  if(memcmp(&state, state_before, sizeof(state)) != 0
     || memcmp(&diagnostics, diagnostics_before, sizeof(diagnostics)) != 0) {
    fail_case("failed repair published state or diagnostics", stage, false);
  }
  stage = fault_none;
}

#endif
