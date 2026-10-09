#ifndef M1_SOURCE_UPDATE_FAULT_INJECTION_H
#define M1_SOURCE_UPDATE_FAULT_INJECTION_H

#include "../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit.h"
#include "../GRHayL/Radiation/ghl_m1_utils.h"
#include "ghl.h"
#include "ghl_m1.h"

#include <dlfcn.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * The library's validated public inputs cannot force every defensive return
 * in the source update. Interpose selected helper results in this executable
 * so the actual shared-library code must reject them transactionally.
 * These cases run only on Linux under __linux__: they require ELF symbol
 * interposition and dlsym(RTLD_NEXT). macOS builds omit the interposition
 * cases; ordinary source-update checks still execute there.
 * M1 exposes no public hook pointers for these cases: mutable/global hooks
 * would conflict with the caller-owned state in the public API. ELF
 * interposition is the accepted Linux-only injection mechanism; the ordinary
 * (non-interposition) source-update tests remain portable to macOS.
 */

typedef enum {
  fault_none,
  fault_late_closure,
  fault_second_repair,
  fault_repair_delegate,
  fault_hn,
  fault_ruu,
  fault_candidate,
  fault_predictor_closure,
  fault_rdd,
  fault_thin_exchange,
  fault_thick_exchange,
  fault_thick_current,
  fault_thick_lepton,
  fault_velocity,
  fault_pair_exchange_first,
  fault_pair_exchange_second
} fault_stage;

static fault_stage stage;
static bool route_thick;
static int closure_calls, moments_calls, repair_calls, realizability_calls, adm_calls;
static int assemble_calls, current_calls, velocity_calls;

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
  if(stage == fault_predictor_closure && closure_calls == 1) {
    return ghl_error_m1_invalid_state;
  }
  if(stage == fault_late_closure && closure_calls == (route_thick ? 2 : 3)) {
    return ghl_error_m1_invalid_state;
  }
  if(stage == fault_hn || stage == fault_ruu || stage == fault_candidate
     || stage == fault_rdd) {
    /* The injected moments have no physical closure tensor. */
    *closure = (ghl_m1_closure){ 0 };
    return ghl_success;
  }
  typedef ghl_error_codes_t (*real_fn)(
        const ghl_m1_parameters *, const ghl_metric_quantities *,
        const ghl_primitive_quantities *, const ghl_m1_rad_state *, ghl_m1_closure *);
  real_fn real = (real_fn)dlsym(RTLD_NEXT, "ghl_m1_compute_closure_with_primitives");
  if(real == NULL) {
    abort();
  }
  return real(m1_params, metric, prims, rad_state, closure);
}

ghl_error_codes_t ghl_m1_compute_comoving_moments_with_velocity(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_comoving *restrict comoving,
      double V_con[3],
      double V_cov[3],
      double *restrict W_out) {
  ++moments_calls;
  if(stage == fault_hn || stage == fault_ruu || stage == fault_candidate
     || stage == fault_rdd) {
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
  real_fn real
        = (real_fn)dlsym(RTLD_NEXT, "ghl_m1_compute_comoving_moments_with_velocity");
  if(real == NULL) {
    abort();
  }
  return real(
        m1_params, metric, prims, rad_state, closure, comoving, V_con, V_cov, W_out);
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
  if(real == NULL) {
    abort();
  }
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
  if(real == NULL) {
    abort();
  }
  return real(m1_params, metric, rad_state);
}

void ghl_compute_ADM_auxiliaries(
      const ghl_metric_quantities *restrict metric,
      ghl_ADM_aux_quantities *restrict metric_aux) {
  ++adm_calls;
  if(stage == fault_rdd) {
    /* Finite inconsistent ADM data make u_0 finite but u_0*u_0 overflow.
     * The predictor must reject R_DD before raising either tensor index. */
    *metric_aux = (ghl_ADM_aux_quantities){ 0 };
    metric_aux->g4DD[0][0] = -DBL_MAX;
    metric_aux->g4UU[0][0] = -1.0;
    for(int i = 1; i < 4; ++i) {
      metric_aux->g4DD[i][i] = 1.0;
      metric_aux->g4UU[i][i] = 1.0;
    }
    return;
  }
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
  if(real == NULL) {
    abort();
  }
  real(metric, metric_aux);
}

/* Selected-call counters keep the forced failures inside the source-update and
 * pair final assembly windows. An untargeted call count trips fail_case so the
 * guards below are proven at the intended stages rather than anywhere the
 * shared helper is invoked. */
ghl_error_codes_t ghl_m1_neutrino_assemble_exchange(
      const ghl_m1_neutrino_state *restrict state_in,
      const ghl_m1_neutrino_state *restrict state_out,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dL_rad_cc,
      const double sqrt_detgamma,
      const double baryon_density_conserved,
      ghl_m1_neutrino_exchange *restrict exchange) {
  ++assemble_calls;
  if(stage == fault_thin_exchange || stage == fault_thick_exchange
     || (stage == fault_pair_exchange_first && assemble_calls == 3)
     || (stage == fault_pair_exchange_second && assemble_calls == 4)) {
    /* Densitized exchange need not be representable after endpoint validation;
     * this forced failure must roll the update or species packet back. */
    return ghl_error_m1_invalid_state;
  }
  typedef ghl_error_codes_t (*real_fn)(
        const ghl_m1_neutrino_state *, const ghl_m1_neutrino_state *,
        const ghl_m1_neutrino_rates *, const double, const double, const double,
        ghl_m1_neutrino_exchange *);
  real_fn real = (real_fn)dlsym(RTLD_NEXT, "ghl_m1_neutrino_assemble_exchange");
  if(real == NULL) {
    abort();
  }
  return real(
        state_in, state_out, rates, dL_rad_cc, sqrt_detgamma, baryon_density_conserved,
        exchange);
}

ghl_error_codes_t ghl_m1_neutrino_derive_current_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_neutrino_current *restrict current) {
  ++current_calls;

  /* The thick branch derives the endpoint current once before the number
   * policy and once inside the lepton delta. Both readings are defensive: the
   * stiff predictor and independent per-species stage already passed validated
   * finite endpoints. */
  if((stage == fault_thick_current && current_calls == 1)
     || (stage == fault_thick_lepton && current_calls == 2)) {
    return ghl_error_m1_invalid_state;
  }
  typedef ghl_error_codes_t (*real_fn)(
        const ghl_m1_parameters *, const ghl_m1_neutrino_parameters *,
        const ghl_metric_quantities *, const ghl_primitive_quantities *,
        const ghl_m1_neutrino_state *, const ghl_m1_closure *,
        ghl_m1_neutrino_current *);
  real_fn real
        = (real_fn)dlsym(RTLD_NEXT, "ghl_m1_neutrino_derive_current_from_closure");
  if(real == NULL) {
    abort();
  }
  return real(m1_params, nu_params, metric, prims, state, closure, current);
}

static ghl_error_codes_t (*ghl_m1_real_finish_scaled_norm_ratio)(
      const double,
      const double,
      const double,
      const double,
      double *) = NULL;

ghl_error_codes_t ghl_m1_finish_scaled_norm_ratio(
      const double x_scale,
      const double A_scale,
      const double denom,
      const double scaled_norm,
      double *restrict ratio) {
  if(ghl_m1_real_finish_scaled_norm_ratio == NULL) {
    typedef ghl_error_codes_t (*real_fn)(
          const double, const double, const double, const double, double *);
    ghl_m1_real_finish_scaled_norm_ratio
          = (real_fn)dlsym(RTLD_NEXT, "ghl_m1_finish_scaled_norm_ratio");
    if(ghl_m1_real_finish_scaled_norm_ratio == NULL) {
      abort();
    }
  }
  /* ghl_m1_compute_eulerian_velocity is inlined, so the scaled-vector norm is
   * its only interposable error conduit. The stiff predictor reads the same
   * velocity through the moments path, so the first norm failure under this
   * stage can only be the thick/scattering dispatch recheck. */
  ++velocity_calls;
  if(stage == fault_velocity && velocity_calls == 1) {
    /* ghl_m1_compute_eulerian_velocity maps any scaled-norm failure to
     * ghl_error_u0_singular unconditionally. */
    return ghl_error_m1_invalid_state;
  }
  return ghl_m1_real_finish_scaled_norm_ratio(
        x_scale, A_scale, denom, scaled_norm, ratio);
}

static void run_case(
      const ghl_m1_parameters *restrict m1_params,
      const fault_stage selected_stage,
      const bool thick) {
  stage = selected_stage;
  route_thick = thick;
  closure_calls = moments_calls = repair_calls = adm_calls = 0;
  assemble_calls = current_calls = velocity_calls = 0;
  ghl_metric_quantities metric;
  if(stage == fault_ruu) {
    ghl_initialize_metric(
          1.0, 0.0, 0.0, 0.0, 1.0e-154, 0.0, 0.0, 1.0, 0.0, 1.0e100, &metric);
  }
  else {
    ghl_initialize_metric(
          stage == fault_candidate ? 2.0 : 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
          1.0, &metric);
  }
  ghl_primitive_quantities prims = { 0 };
  if(stage == fault_velocity) {
    /* A subluminal velocity keeps the physics valid; the interposed scaled-norm
     * failure below is the only forced condition at this selected call. */
    prims.vU[0] = 0.25;
    prims.u0 = 1.0;
  }
  ghl_m1_neutrino_parameters params
        = { .N_floor = 1.0e-12, .J_floor = 1.0e-14, .Gamma_N_floor = 1.0e-12 };
  ghl_m1_neutrino_rates rates
        = { .species = ghl_m1_neutrino_nux,
            .mean_energy = 1.0,
            .n_eq = stage == fault_candidate ? 0.0 : 1.0,
            .J_eq = stage == fault_candidate ? 0.0 : 1.0,
            .kappa_a_E = thick ? 1.0 : 0.0,
            .kappa_s = thick ? 1.0 : 0.0,
            .kappa_tr = thick ? 2.0 : 0.0,
            .eta_E = thick && stage != fault_candidate ? 1.0 : 0.0 };
  ghl_m1_neutrino_source_options options
        = { .policy = ghl_m1_neutrino_source_branched_compatibility,
            .thick_equilibrium_threshold = 0.5,
            .scattering_threshold = 0.5,
            .thermalized_number_threshold = -1.0,
            .allow_closure_fallback = true,
            .ye_policy = ghl_m1_neutrino_ye_from_charged_current };
  const ghl_m1_neutrino_state input = { .E = 1.0, .N = 1.0 };
  ghl_m1_neutrino_state output = { 0 };
  ghl_m1_neutrino_exchange exchange;
  ghl_m1_neutrino_source_diagnostics diagnostics;
  ghl_m1_neutrino_diagnostics nd = { 0 };
  const ghl_error_codes_t error = ghl_m1_solve_neutrino_source_update(
        &options, m1_params, &params, &metric, &prims, &rates, &input, &input, 2.0, 1.0,
        &output, &exchange, &diagnostics, &nd);
  const ghl_error_codes_t expected_error
        = stage == fault_velocity ? ghl_error_u0_singular : ghl_error_m1_invalid_state;
  if(error != expected_error
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
  if(exchange.dN_rad_total != 0.0 || exchange.dL_rad_cc != 0.0 || exchange.dE_rad != 0.0
     || exchange.dTau_matter != 0.0 || exchange.dYe_matter != 0.0) {
    fail_case("partial exchange was published", stage, thick);
  }
  if((stage == fault_late_closure && closure_calls != (thick ? 2 : 3))
     || (stage == fault_second_repair && repair_calls != 2)
     || ((stage == fault_hn || stage == fault_ruu || stage == fault_candidate
          || stage == fault_rdd)
         && moments_calls != 1)
     || (stage == fault_hn && adm_calls != 0)
     || ((stage == fault_ruu || stage == fault_candidate || stage == fault_rdd)
         && adm_calls != 1)
     || (stage == fault_predictor_closure
         && (closure_calls != 1 || moments_calls != 0 || adm_calls != 0))
     || ((stage == fault_predictor_closure || stage == fault_rdd) && repair_calls != 0)
     || (stage == fault_thin_exchange && assemble_calls != 1)
     || ((stage == fault_thick_current || stage == fault_thick_lepton)
         && current_calls != (stage == fault_thick_current ? 1 : 2))
     || (stage == fault_thick_exchange && (assemble_calls != 1 || current_calls != 2))
     || (stage == fault_velocity
         && (velocity_calls != 1 || adm_calls != 0 || repair_calls != 0))) {
    fail_case("injection did not reach its intended stage", stage, thick);
  }
}

static void
test_source_update_fault_injections(const ghl_m1_parameters *restrict m1_params) {
  run_case(m1_params, fault_late_closure, false);
  run_case(m1_params, fault_late_closure, true);
  run_case(m1_params, fault_second_repair, false);
  run_case(m1_params, fault_second_repair, true);
  run_case(m1_params, fault_hn, true);
  run_case(m1_params, fault_ruu, true);
  run_case(m1_params, fault_candidate, true);
  run_case(m1_params, fault_predictor_closure, true);
  run_case(m1_params, fault_rdd, true);
  run_case(m1_params, fault_thin_exchange, false);
  run_case(m1_params, fault_thick_current, true);
  run_case(m1_params, fault_thick_lepton, true);
  run_case(m1_params, fault_thick_exchange, true);
  run_case(m1_params, fault_velocity, true);
  stage = fault_none;
}

/*
 * The four counts below are the two independent per-species exchanges plus the
 * two final pair exchanges. The injection fires only on the final pair reads,
 * after both independent updates and the pair schedule all completed, so the
 * forced failure can only exercise the pair final rollback. */
static void run_pair_case(
      const ghl_m1_parameters *restrict m1_params,
      const fault_stage selected_stage,
      const int expected_error_calls) {
  stage = selected_stage;
  assemble_calls = 0;
  ghl_metric_quantities metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &metric);
  ghl_primitive_quantities prims = { 0 };
  ghl_m1_neutrino_parameters nu_params[2]
        = { { .N_floor = 1.0e-12, .J_floor = 1.0e-14, .Gamma_N_floor = 1.0e-12 },
            { .N_floor = 1.0e-12, .J_floor = 1.0e-14, .Gamma_N_floor = 1.0e-12 } };
  ghl_m1_neutrino_rates rates[2];
  make_rates(ghl_m1_neutrino_nue, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, &rates[0]);
  make_rates(ghl_m1_neutrino_anue, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, &rates[1]);
  for(int process = 0; process < 3; ++process) {
    rates[0].eta_N_pair[process] = 0.1;
    rates[1].eta_N_pair[process] = 0.1;
    rates[0].eta_E_pair[process] = 0.2;
    rates[1].eta_E_pair[process] = 0.2;
  }
  const ghl_m1_neutrino_state state_transport[2]
        = { { .N = 0.5, .E = 1.0 }, { .N = 1.5, .E = 3.0 } };
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
        m1_params, nu_params, &metric, &prims, rates, state_input, state_transport, 0.25,
        3.0, state_out, exchange, diagnostics, neutrino_diagnostics);
  if(error != ghl_error_m1_invalid_state) {
    fail_case("pair final exchange failure did not propagate", stage, false);
  }
  for(int species = 0; species < 2; ++species) {
    if(diagnostics[species].path != ghl_m1_neutrino_source_path_hard_failure
       || diagnostics[species].terminal_no_update
       || neutrino_diagnostics[species].source_failures != 1) {
      fail_case("pair rollback diagnostics contract failed", stage, false);
    }
    if(state_out[species].N != state_transport[species].N
       || state_out[species].E != state_transport[species].E) {
      fail_case("pair failure lifted a partial state", stage, false);
    }
    for(int direction = 0; direction < 3; ++direction) {
      if(state_out[species].F[direction] != state_transport[species].F[direction]
         || exchange[species].dF_rad[direction] != 0.0
         || exchange[species].dS_matter[direction] != 0.0) {
        fail_case(
              "pair failure published a partial flux or momentum exchange", stage,
              false);
      }
    }
    if(exchange[species].dN_rad_total != 0.0 || exchange[species].dL_rad_cc != 0.0
       || exchange[species].dE_rad != 0.0 || exchange[species].dTau_matter != 0.0
       || exchange[species].dYe_matter != 0.0) {
      fail_case("pair failure published a partial exchange", stage, false);
    }
  }
  if(assemble_calls != expected_error_calls) {
    fail_case(
          "pair final exchange injection did not reach the intended species", stage,
          false);
  }
  stage = fault_none;
}

static void
test_pair_final_exchange_failures(const ghl_m1_parameters *restrict m1_params) {
  run_pair_case(m1_params, fault_pair_exchange_first, 3);
  run_pair_case(m1_params, fault_pair_exchange_second, 4);
}

static void
test_neutrino_repair_delegate_failure(const ghl_m1_parameters *restrict m1_params) {
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
