#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"
#include "../ghl_m1_utils.h"

#include <float.h>

/*
 * Host-neutral source-policy dispatcher for one grey neutrino species.
 *
 * The established GRHayL implicit solver remains the zero-valued/default
 * policy. The opt-in branched policy is deliberately a thin dispatch layer:
 * branch selection is based only on the frozen stage rates and
 * dt_alpha = alpha * dt; all candidates are kept private until their endpoint
 * current and exchange packet have
 * validated, and the source base is always state_transport.  In particular,
 * state_input is retained for the public pre-transport validation contract;
 * every branched source candidate is evaluated from state_transport. This
 * uses the post-transport predictor as both its old and source-base radiation
 * state.
 *
 * No host schedule, grid, matter update, or stage-global state belongs here.
 */

static void ghl_m1_neutrino_zero_exchange(
      ghl_m1_neutrino_exchange *restrict exchange) {
  if(exchange != NULL)
    *exchange = (ghl_m1_neutrino_exchange){0};
}

static void ghl_m1_neutrino_initialize_source_diagnostics(
      ghl_m1_neutrino_source_diagnostics *restrict diagnostics) {
  if(diagnostics == NULL)
    return;
  *diagnostics = (ghl_m1_neutrino_source_diagnostics){0};
  ghl_m1_initialize_implicit_solve_diagnostics(&diagnostics->implicit);
}

static bool ghl_m1_neutrino_state_is_finite(
      const ghl_m1_neutrino_state *restrict state) {
  if(state == NULL || !isfinite(state->N) || !isfinite(state->E))
    return false;
  for(int i = 0; i < 3; ++i)
    if(!isfinite(state->F[i]))
      return false;
  return true;
}

static ghl_error_codes_t ghl_m1_neutrino_validate_state_input(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL || state == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(nu_params->N_floor) || nu_params->N_floor < 0.0 ||
     !isfinite(nu_params->J_floor) || nu_params->J_floor < 0.0 ||
     !isfinite(nu_params->Gamma_N_floor) || nu_params->Gamma_N_floor < 0.0)
    return ghl_error_m1_invalid_state;
  if(!ghl_m1_neutrino_state_is_finite(state) ||
     state->N < nu_params->N_floor)
    return ghl_error_m1_invalid_state;
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  return ghl_m1_validate_realizability(m1_params, metric, &rad_state, 128.0,
                                       NULL);
}

static bool ghl_m1_neutrino_closure_fallback_status(
      const ghl_m1_closure *restrict closure) {
  return closure != NULL &&
      (closure->solve_status == ghl_m1_closure_solve_endpoint_fallback ||
       closure->solve_status == ghl_m1_closure_solve_iteration_exhausted);
}

static ghl_error_codes_t ghl_m1_neutrino_evaluate_closure_and_sources(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_closure *restrict closure,
      bool *restrict closure_fallback,
      ghl_m1_sources *restrict EF_sources,
      double *restrict N_source) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL || prims == NULL ||
     state == NULL || rates == NULL || closure == NULL ||
     closure_fallback == NULL || EF_sources == NULL || N_source == NULL)
    return ghl_error_m1_null_pointer;

  *closure_fallback = false;
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  ghl_error_codes_t error = ghl_m1_compute_closure_with_primitives(
      m1_params, metric, prims, &rad_state, closure);
  if(error != ghl_success)
    return error;
  *closure_fallback = ghl_m1_neutrino_closure_fallback_status(closure);

  return ghl_m1_compute_neutrino_interaction_sources_from_closure(
      m1_params, nu_params, metric, prims, state, closure, rates,
      EF_sources, N_source);
}

static ghl_error_codes_t ghl_m1_neutrino_compute_endpoint_lepton_delta(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const ghl_m1_neutrino_state *restrict state_out,
      bool *restrict closure_fallback,
      double *restrict dL_rad_cc) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL || prims == NULL ||
     rates == NULL || state_out == NULL || closure_fallback == NULL ||
     dL_rad_cc == NULL)
    return ghl_error_m1_null_pointer;

  *closure_fallback = false;
  *dL_rad_cc = 0.0;
  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state_out);
  ghl_m1_closure closure;
  const ghl_error_codes_t closure_error = ghl_m1_compute_closure_with_primitives(
      m1_params, metric, prims, &rad_state, &closure);
  if(closure_error != ghl_success)
    return closure_error;
  *closure_fallback = ghl_m1_neutrino_closure_fallback_status(&closure);

  ghl_m1_neutrino_current current;
  const ghl_error_codes_t current_error =
      ghl_m1_neutrino_derive_current_from_closure(
          m1_params, nu_params, metric, prims, state_out, &closure, &current);
  if(current_error != ghl_success)
    return current_error;

  const ghl_error_codes_t bounds_error =
      ghl_m1_neutrino_check_EN_bounds(state_out, nu_params, &current);
  if(bounds_error != ghl_success)
    return bounds_error;

  const double dt_alpha = metric->lapse * dt;
  if(!isfinite(dt_alpha) || dt_alpha < 0.0 ||
     !isfinite(current.Gamma_N) || current.Gamma_N <= 0.0)
    return ghl_error_m1_invalid_state;

  const double dN_cc = dt_alpha *
      (rates->eta_N_cc - rates->kappa_a_N_cc * state_out->N /
       current.Gamma_N);
  const double candidate = rates->lepton_weight * dN_cc;
  if(!isfinite(candidate))
    return ghl_error_m1_invalid_state;
  *dL_rad_cc = candidate;
  return ghl_success;
}

static ghl_error_codes_t ghl_m1_neutrino_repair_candidate(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      ghl_m1_neutrino_state *restrict state,
      ghl_m1_neutrino_diagnostics *restrict diagnostics) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL || state == NULL ||
     diagnostics == NULL)
    return ghl_error_m1_null_pointer;
  return ghl_m1_repair_neutrino_state(
      m1_params, nu_params, metric, state, diagnostics);
}

static ghl_error_codes_t ghl_m1_neutrino_build_exchange(
      const ghl_m1_neutrino_state *restrict state_transport,
      const ghl_m1_neutrino_state *restrict state_candidate,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dL_rad_cc,
      const ghl_metric_quantities *restrict metric,
      const double n_b_cons,
      ghl_m1_neutrino_exchange *restrict exchange) {
  if(state_transport == NULL || state_candidate == NULL || rates == NULL ||
     metric == NULL || exchange == NULL)
    return ghl_error_m1_null_pointer;
  return ghl_m1_neutrino_assemble_exchange(
      state_transport, state_candidate, rates, dL_rad_cc,
      metric->sqrt_detgamma, n_b_cons, exchange);
}

static ghl_error_codes_t ghl_m1_neutrino_apply_ye_policy(
      const ghl_m1_neutrino_ye_policy_t ye_policy,
      const ghl_m1_neutrino_rates *restrict rates,
      const double n_b_cons,
      ghl_m1_neutrino_exchange *restrict exchange) {
  if(rates == NULL || exchange == NULL)
    return ghl_error_m1_null_pointer;
  if(ye_policy == ghl_m1_neutrino_ye_from_charged_current)
    return ghl_success;
  if(ye_policy != ghl_m1_neutrino_ye_from_signed_total_number)
    return ghl_error_m1_invalid_state;

  /* The three-species composition source is +DrN_nue-DrN_anue.  The
   * exchange stores an unsigned per-species number increment, so the
   * validated species lepton weight supplies the sign here. */
  const double dL_rad_total = rates->lepton_weight * exchange->dN_rad_total;
  return ghl_m1_compute_neutrino_lepton_increment(
      rates, dL_rad_total, n_b_cons, &exchange->dYe_matter);
}

static ghl_error_codes_t ghl_m1_neutrino_update_endpoint_number(
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double dt_alpha,
      const double thermalized_number_threshold,
      const ghl_m1_neutrino_state *restrict state_base,
      const ghl_m1_neutrino_current *restrict endpoint_current,
      double *restrict N_out) {
  if(nu_params == NULL || rates == NULL || state_base == NULL ||
     endpoint_current == NULL || N_out == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(dt) || dt < 0.0 || !isfinite(dt_alpha) || dt_alpha < 0.0 ||
     !isfinite(thermalized_number_threshold))
    return ghl_error_m1_invalid_state;

  /* The endpoint neutrino mean energy may be treated as thermalized.  A
   * negative threshold disables that option, preserving the
   * ordinary endpoint-Gamma backward-Euler update. */
  const double number_stiffness = dt_alpha * rates->kappa_a_N;
  if(thermalized_number_threshold >= 0.0 &&
     number_stiffness >= thermalized_number_threshold) {
    const double candidate = rates->mean_energy > 0.0
        ? endpoint_current->Gamma_N * endpoint_current->J / rates->mean_energy
        : 0.0;
    if(!isfinite(candidate))
      return ghl_error_m1_invalid_state;
    *N_out = candidate;
    return ghl_success;
  }

  return ghl_m1_update_neutrino_number_backward_euler(
      nu_params, rates, dt_alpha, endpoint_current->Gamma_N,
      state_base->N, N_out);
}

/*
 * This is the default local predictor for the thick/scattering shortcuts. It
 * first advances transport to Estar/Fstar,
 * transforms that state to the fluid frame, applies backward-Euler damping to
 * J and H, assumes chi=1/3, and boosts the tensor back to the Eulerian frame.
 * source_update then returns THICK/SCAT without replacing this predictor by a
 * second equilibrium or explicit source update. In this API state_transport is
 * that post-transport predictor.
 */
static ghl_error_codes_t ghl_m1_neutrino_build_stiff_predictor(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_transport,
      const double dt,
      ghl_m1_neutrino_diagnostics *restrict diagnostics,
      bool *restrict closure_fallback_used,
      ghl_m1_neutrino_state *restrict state) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL ||
     prims == NULL || rates == NULL || state_transport == NULL ||
     diagnostics == NULL || closure_fallback_used == NULL || state == NULL)
    return ghl_error_m1_null_pointer;

  const double dt_alpha = metric->lapse * dt;
  if(!isfinite(dt) || dt < 0.0 || !isfinite(dt_alpha) || dt_alpha < 0.0)
    return ghl_error_m1_invalid_state;

  *closure_fallback_used = false;
  const ghl_m1_rad_state transport_rad =
      ghl_m1_neutrino_project_rad_state(state_transport);
  ghl_m1_closure transport_closure;
  const ghl_error_codes_t closure_error = ghl_m1_compute_closure_with_primitives(
      m1_params, metric, prims, &transport_rad, &transport_closure);
  if(closure_error != ghl_success)
    return closure_error;
  *closure_fallback_used = ghl_m1_neutrino_closure_fallback_status(
      &transport_closure);

  ghl_m1_comoving comoving;
  double V_con[3], V_cov[3], W;
  const ghl_error_codes_t moments_error =
      ghl_m1_compute_comoving_moments_with_velocity(
          m1_params, metric, prims, &transport_rad, &transport_closure,
          &comoving, V_con, V_cov, &W);
  if(moments_error != ghl_success)
    return moments_error;

  const double dtau = dt_alpha / W;
  const double J_denominator = 1.0 + dtau * rates->kappa_a_E;
  const double H_denominator = 1.0 + dtau * rates->kappa_tr;
  if(!isfinite(dtau) || dtau < 0.0 || !isfinite(J_denominator) ||
     J_denominator <= 0.0 || !isfinite(H_denominator) ||
     H_denominator <= 0.0)
    return ghl_error_m1_invalid_state;

  const double J_new = (comoving.J + dtau * rates->eta_E) /
                       J_denominator;
  if(!isfinite(J_new) || J_new < 0.0)
    return ghl_error_m1_invalid_state;

  double HD_new[3], HU_new[3];
  for(int i = 0; i < 3; ++i)
    HD_new[i] = comoving.HD[i] / H_denominator;
  ghl_raise_lower_vector_3D(metric->gammaUU, HD_new, HU_new);
  double Hn_new = 0.0;
  for(int i = 0; i < 3; ++i)
    Hn_new -= V_cov[i] * HU_new[i];
  if(!isfinite(Hn_new))
    return ghl_error_m1_invalid_state;

  /* Form the same ADM 4-metric used by GRHayL's stress-energy kernels. */
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(metric, &metric_aux);

  const double alpha = metric->lapse;
  const double inv_alpha = 1.0 / alpha;
  double uU[4] = {W * inv_alpha, 0.0, 0.0, 0.0};
  for(int i = 0; i < 3; ++i)
    uU[i + 1] = W * (V_con[i] - metric->betaU[i] * inv_alpha);

  double uD[4];
  ghl_raise_lower_vector_4D(metric_aux.g4DD, uU, uD);
  const double hU[4] = {0.0, HU_new[0], HU_new[1], HU_new[2]};
  double hD[4];
  ghl_raise_lower_vector_4D(metric_aux.g4DD, hU, hD);
  const double nD[4] = {-alpha, 0.0, 0.0, 0.0};
  double H_D[4];
  for(int mu = 0; mu < 4; ++mu)
    H_D[mu] = hD[mu] - Hn_new * nD[mu];

  /* The boosted predictor's closure is fixed to chi=1/3. */
  double R_DD[4][4] = {{0.0}};
  for(int mu = 0; mu < 4; ++mu) {
    for(int nu = 0; nu < 4; ++nu) {
      const double K_DD = (J_new / 3.0) *
          (metric_aux.g4DD[mu][nu] + uD[mu] * uD[nu]);
      R_DD[mu][nu] = J_new * uD[mu] * uD[nu]
                   + H_D[mu] * uD[nu] + H_D[nu] * uD[mu]
                   + K_DD;
      if(!isfinite(R_DD[mu][nu]))
        return ghl_error_m1_invalid_state;
    }
  }

  double R_UU[4][4] = {{0.0}};
  for(int mu = 0; mu < 4; ++mu) {
    for(int nu = 0; nu < 4; ++nu) {
      for(int a = 0; a < 4; ++a) {
        for(int b = 0; b < 4; ++b)
          R_UU[mu][nu] += metric_aux.g4UU[mu][a]
                        * metric_aux.g4UU[nu][b] * R_DD[a][b];
      }
      if(!isfinite(R_UU[mu][nu]))
        return ghl_error_m1_invalid_state;
    }
  }

  double betaD[3];
  ghl_raise_lower_vector_3D(metric->gammaDD, metric->betaU, betaD);
  ghl_m1_neutrino_state candidate = *state_transport;
  candidate.E = alpha * alpha * R_UU[0][0];
  for(int i = 0; i < 3; ++i) {
    candidate.F[i] = alpha * (betaD[i] * R_UU[0][0]
        + metric->gammaDD[i][0] * R_UU[0][1]
        + metric->gammaDD[i][1] * R_UU[0][2]
        + metric->gammaDD[i][2] * R_UU[0][3]);
  }
  if(!ghl_m1_neutrino_state_is_finite(&candidate))
    return ghl_error_m1_invalid_state;

  *state = candidate;
  return ghl_m1_neutrino_repair_candidate(
      m1_params, nu_params, metric, state, diagnostics);
}

static ghl_error_codes_t ghl_m1_neutrino_try_thin_branch(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_transport,
      const double dt,
      const double n_b_cons,
      const double thermalized_number_threshold,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      bool *restrict closure_fallback_used) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL ||
     prims == NULL || rates == NULL || state_transport == NULL ||
     neutrino_diagnostics == NULL || state_out == NULL || exchange == NULL ||
     closure_fallback_used == NULL)
    return ghl_error_m1_null_pointer;

  *closure_fallback_used = false;
  ghl_m1_closure transport_closure;
  ghl_m1_sources EF_sources;
  double ignored_N_source;
  bool transport_fallback = false;
  ghl_error_codes_t error = ghl_m1_neutrino_evaluate_closure_and_sources(
      m1_params, nu_params, metric, prims, state_transport, rates,
      &transport_closure, &transport_fallback, &EF_sources,
      &ignored_N_source);
  if(error != ghl_success)
    return error;
  *closure_fallback_used = transport_fallback;

  const double dt_alpha = metric->lapse * dt;
  if(!isfinite(dt_alpha) || dt_alpha < 0.0)
    return ghl_error_m1_invalid_state;

  ghl_m1_neutrino_state candidate = *state_transport;
  candidate.E += dt_alpha * EF_sources.S_E;
  for(int i = 0; i < 3; ++i)
    candidate.F[i] += dt_alpha * EF_sources.S[i];
  error = ghl_m1_neutrino_repair_candidate(
      m1_params, nu_params, metric, &candidate, neutrino_diagnostics);
  if(error != ghl_success)
    return error;

  /* The implicit thin branch computes the number update after the explicit
   * E/F endpoint, using that endpoint's Gamma_N and backward Euler. */
  ghl_m1_neutrino_current endpoint_current;
  error = ghl_m1_neutrino_derive_current(
      m1_params, nu_params, metric, prims, &candidate, &endpoint_current);
  if(error != ghl_success)
    return error;
  error = ghl_m1_neutrino_update_endpoint_number(
      nu_params, rates, dt, dt_alpha, thermalized_number_threshold,
      state_transport, &endpoint_current, &candidate.N);
  if(error != ghl_success)
    return error;
  error = ghl_m1_neutrino_repair_candidate(
      m1_params, nu_params, metric, &candidate, neutrino_diagnostics);
  if(error != ghl_success)
    return error;

  bool endpoint_fallback = false;
  double dL_rad_cc = 0.0;
  error = ghl_m1_neutrino_compute_endpoint_lepton_delta(
      m1_params, nu_params, metric, prims, rates, dt, &candidate,
      &endpoint_fallback, &dL_rad_cc);
  if(error != ghl_success)
    return error;
  *closure_fallback_used |= endpoint_fallback;

  error = ghl_m1_neutrino_build_exchange(
      state_transport, &candidate, rates, dL_rad_cc, metric, n_b_cons,
      exchange);
  if(error != ghl_success)
    return error;
  *state_out = candidate;
  return ghl_success;
}

static ghl_error_codes_t ghl_m1_neutrino_try_thick_branch(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_transport,
      const double dt,
      const double n_b_cons,
      const double thermalized_number_threshold,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      bool *restrict closure_fallback_used) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL ||
     prims == NULL || rates == NULL || state_transport == NULL ||
     neutrino_diagnostics == NULL || state_out == NULL || exchange == NULL ||
     closure_fallback_used == NULL)
    return ghl_error_m1_null_pointer;

  /* The default source method uses the same local boost predictor for
   * both its THICK and SCAT status returns. */
  ghl_m1_neutrino_state candidate;
  ghl_error_codes_t error = ghl_m1_neutrino_build_stiff_predictor(
      m1_params, nu_params, metric, prims, rates, state_transport, dt,
      neutrino_diagnostics, closure_fallback_used, &candidate);
  if(error != ghl_success)
    return error;

  const double dt_alpha = metric->lapse * dt;
  if(!isfinite(dt_alpha) || dt_alpha < 0.0)
    return ghl_error_m1_invalid_state;

  /* Update N after the E/F predictor, using the endpoint current. */
  bool endpoint_fallback = false;
  const ghl_m1_rad_state endpoint_rad =
      ghl_m1_neutrino_project_rad_state(&candidate);
  ghl_m1_closure endpoint_closure;
  error = ghl_m1_compute_closure_with_primitives(
      m1_params, metric, prims, &endpoint_rad, &endpoint_closure);
  if(error != ghl_success)
    return error;
  endpoint_fallback = ghl_m1_neutrino_closure_fallback_status(&endpoint_closure);
  *closure_fallback_used |= endpoint_fallback;

  ghl_m1_neutrino_current endpoint_current;
  error = ghl_m1_neutrino_derive_current_from_closure(
      m1_params, nu_params, metric, prims, &candidate, &endpoint_closure,
      &endpoint_current);
  if(error != ghl_success)
    return error;

  double N_out = 0.0;
  error = ghl_m1_neutrino_update_endpoint_number(
      nu_params, rates, dt, dt_alpha, thermalized_number_threshold,
      state_transport, &endpoint_current, &N_out);
  if(error != ghl_success)
    return error;
  candidate.N = N_out;
  error = ghl_m1_neutrino_repair_candidate(
      m1_params, nu_params, metric, &candidate, neutrino_diagnostics);
  if(error != ghl_success)
    return error;

  bool lepton_fallback = false;
  double dL_rad_cc = 0.0;
  error = ghl_m1_neutrino_compute_endpoint_lepton_delta(
      m1_params, nu_params, metric, prims, rates, dt, &candidate,
      &lepton_fallback, &dL_rad_cc);
  if(error != ghl_success)
    return error;
  *closure_fallback_used |= lepton_fallback;

  error = ghl_m1_neutrino_build_exchange(
      state_transport, &candidate, rates, dL_rad_cc, metric, n_b_cons,
      exchange);
  if(error != ghl_success)
    return error;
  *state_out = candidate;
  return ghl_success;
}

static ghl_error_codes_t ghl_m1_neutrino_try_scattering_branch(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_transport,
      const double dt,
      const double n_b_cons,
      const double thermalized_number_threshold,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      bool *restrict closure_fallback_used) {
  if(m1_params == NULL || nu_params == NULL || metric == NULL ||
     prims == NULL || rates == NULL || state_transport == NULL ||
     neutrino_diagnostics == NULL || state_out == NULL || exchange == NULL ||
     closure_fallback_used == NULL)
    return ghl_error_m1_null_pointer;

  return ghl_m1_neutrino_try_thick_branch(
      m1_params, nu_params, metric, prims, rates, state_transport, dt,
      n_b_cons, thermalized_number_threshold, neutrino_diagnostics,
      state_out, exchange,
      closure_fallback_used);
}

static ghl_error_codes_t ghl_m1_neutrino_publish_dispatch_failure(
      const ghl_error_codes_t error,
      const ghl_m1_neutrino_state *restrict state_transport,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      ghl_m1_neutrino_source_diagnostics *restrict diagnostics,
      const bool closure_fallback_used,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics,
      const bool count_source_failure) {
  if(state_transport != NULL && state_out != NULL)
    *state_out = *state_transport;
  ghl_m1_neutrino_zero_exchange(exchange);
  if(diagnostics != NULL) {
    diagnostics->path = ghl_m1_neutrino_source_path_hard_failure;
    diagnostics->closure_fallback_used = closure_fallback_used;
    diagnostics->terminal_no_update = false;
  }
  if(neutrino_diagnostics != NULL && count_source_failure)
    neutrino_diagnostics->source_failures++;
  return error;
}

static bool ghl_m1_neutrino_thick_limit_selected(
      const double dt_alpha,
      const ghl_m1_neutrino_rates *restrict rates,
      const double threshold) {
  if(threshold <= 0.0)
    return false;
  const double opacity_product = rates->kappa_a_E * rates->kappa_tr;
  if(!isfinite(opacity_product) || opacity_product <= 0.0)
    return false;
  const double stiffness = dt_alpha * sqrt(opacity_product);
  return isfinite(stiffness) ? stiffness > threshold : true;
}

static bool ghl_m1_neutrino_scattering_limit_selected(
      const double dt_alpha,
      const ghl_m1_neutrino_rates *restrict rates,
      const double threshold) {
  if(threshold <= 0.0)
    return false;
  const double stiffness = dt_alpha * rates->kappa_s;
  return isfinite(stiffness) ? stiffness > threshold : true;
}

ghl_error_codes_t ghl_m1_solve_neutrino_source_update(
      const ghl_m1_neutrino_source_options *restrict options,
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_input,
      const ghl_m1_neutrino_state *restrict state_transport,
      const double dt,
      const double n_b_cons,
      ghl_m1_neutrino_state *restrict state_out,
      ghl_m1_neutrino_exchange *restrict exchange,
      ghl_m1_neutrino_source_diagnostics *restrict diagnostics,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics) {
  if(state_out == NULL || exchange == NULL || diagnostics == NULL ||
     neutrino_diagnostics == NULL || state_transport == NULL)
    return ghl_error_m1_null_pointer;

  *state_out = *state_transport;
  ghl_m1_neutrino_zero_exchange(exchange);
  ghl_m1_neutrino_initialize_source_diagnostics(diagnostics);

  ghl_m1_neutrino_source_options selected = {
      .policy = ghl_m1_neutrino_source_grhayl_implicit,
      .thick_equilibrium_threshold = 0.0,
      .scattering_threshold = 0.0,
      .thermalized_number_threshold = -1.0,
      .allow_closure_fallback = false,
      .interaction_sources_already_applied = false,
      .ye_policy = ghl_m1_neutrino_ye_from_charged_current};
  if(options != NULL)
    selected = *options;

  if(selected.interaction_sources_already_applied)
    return ghl_m1_neutrino_publish_dispatch_failure(
        ghl_error_m1_source_double_application, state_transport, state_out,
        exchange, diagnostics, false, neutrino_diagnostics, true);
  if(selected.policy != ghl_m1_neutrino_source_grhayl_implicit &&
     selected.policy != ghl_m1_neutrino_source_branched_compatibility)
    return ghl_m1_neutrino_publish_dispatch_failure(
        ghl_error_m1_invalid_state, state_transport, state_out, exchange,
        diagnostics, false, neutrino_diagnostics, true);
  if(selected.ye_policy != ghl_m1_neutrino_ye_from_charged_current &&
     selected.ye_policy != ghl_m1_neutrino_ye_from_signed_total_number)
    return ghl_m1_neutrino_publish_dispatch_failure(
        ghl_error_m1_invalid_state, state_transport, state_out, exchange,
        diagnostics, false, neutrino_diagnostics, true);
  if(!isfinite(selected.thick_equilibrium_threshold) ||
     !isfinite(selected.scattering_threshold) ||
     !isfinite(selected.thermalized_number_threshold) ||
     !isfinite(dt) || dt < 0.0 ||
     !isfinite(n_b_cons) || n_b_cons <= 0.0)
    return ghl_m1_neutrino_publish_dispatch_failure(
        ghl_error_m1_invalid_state, state_transport, state_out, exchange,
        diagnostics, false, neutrino_diagnostics, true);

  if(m1_params == NULL || nu_params == NULL || metric == NULL ||
     prims_frozen == NULL || rates == NULL || state_input == NULL)
    return ghl_m1_neutrino_publish_dispatch_failure(
        ghl_error_m1_null_pointer, state_transport, state_out, exchange,
        diagnostics, false, neutrino_diagnostics, true);

  ghl_error_codes_t error = ghl_m1_validate_configuration(m1_params, metric);
  if(error != ghl_success)
    return ghl_m1_neutrino_publish_dispatch_failure(
        error, state_transport, state_out, exchange, diagnostics, false,
        neutrino_diagnostics, true);
  const double dt_alpha = metric->lapse * dt;
  if(!isfinite(dt_alpha) || dt_alpha < 0.0)
    return ghl_m1_neutrino_publish_dispatch_failure(
        ghl_error_m1_invalid_state, state_transport, state_out, exchange,
        diagnostics, false, neutrino_diagnostics, true);
  error = ghl_m1_neutrino_validate_state_input(
      m1_params, nu_params, metric, state_input);
  if(error != ghl_success)
    return ghl_m1_neutrino_publish_dispatch_failure(
        error, state_transport, state_out, exchange, diagnostics, false,
        neutrino_diagnostics, true);
  error = ghl_m1_neutrino_validate_state_input(
      m1_params, nu_params, metric, state_transport);
  if(error != ghl_success)
    return ghl_m1_neutrino_publish_dispatch_failure(
        error, state_transport, state_out, exchange, diagnostics, false,
        neutrino_diagnostics, true);

  ghl_m1_neutrino_diagnostics candidate_neutrino_diagnostics =
      *neutrino_diagnostics;
  error = ghl_m1_neutrino_validate_single_species_rates(
      rates, &candidate_neutrino_diagnostics);
  if(error != ghl_success) {
    candidate_neutrino_diagnostics.source_failures++;
    *neutrino_diagnostics = candidate_neutrino_diagnostics;
    return ghl_m1_neutrino_publish_dispatch_failure(
        error, state_transport, state_out, exchange, diagnostics, false,
        neutrino_diagnostics, false);
  }

  if(selected.policy == ghl_m1_neutrino_source_grhayl_implicit) {
    ghl_m1_implicit_solve_diagnostics implicit_diagnostics;
    ghl_m1_initialize_implicit_solve_diagnostics(&implicit_diagnostics);
    error = ghl_m1_solve_neutrino_implicit_homogeneous_update(
        m1_params, nu_params, metric, prims_frozen, rates, dt, n_b_cons,
        state_transport, state_out, exchange, &implicit_diagnostics,
        &candidate_neutrino_diagnostics);
    diagnostics->implicit = implicit_diagnostics;
    diagnostics->closure_fallback_used =
        (implicit_diagnostics.solution_path_flags &
         ghl_m1_solution_path_closure_fallback) != 0u;
    if(error == ghl_success) {
      error = ghl_m1_neutrino_apply_ye_policy(
          selected.ye_policy, rates, n_b_cons, exchange);
      if(error != ghl_success) {
        *neutrino_diagnostics = candidate_neutrino_diagnostics;
        return ghl_m1_neutrino_publish_dispatch_failure(
            error, state_transport, state_out, exchange, diagnostics,
            diagnostics->closure_fallback_used, neutrino_diagnostics, true);
      }
      diagnostics->path = ghl_m1_neutrino_source_path_general_implicit;
      *neutrino_diagnostics = candidate_neutrino_diagnostics;
      return ghl_success;
    }
    if(error == ghl_error_m1_implicit_terminal_fallback) {
      diagnostics->path = ghl_m1_neutrino_source_path_terminal_no_update;
      diagnostics->terminal_no_update = true;
      *state_out = *state_transport;
      ghl_m1_neutrino_zero_exchange(exchange);
      *neutrino_diagnostics = candidate_neutrino_diagnostics;
      return error;
    }
    diagnostics->path = ghl_m1_neutrino_source_path_hard_failure;
    *state_out = *state_transport;
    ghl_m1_neutrino_zero_exchange(exchange);
    *neutrino_diagnostics = candidate_neutrino_diagnostics;
    return error;
  }

  /* Compatibility branches are fail-closed on closure fallback unless
   * the host explicitly opts into the fallback candidate.  The pre-check
   * catches the input closure; the general path checks the full Newton trace
   * below because trial states can select a different endpoint status. */
  bool preclosure_fallback = false;
  if(!selected.allow_closure_fallback) {
    ghl_m1_closure preclosure;
    const ghl_m1_rad_state transport_rad =
        ghl_m1_neutrino_project_rad_state(state_transport);
    error = ghl_m1_compute_closure_with_primitives(
        m1_params, metric, prims_frozen, &transport_rad, &preclosure);
    if(error != ghl_success)
      return ghl_m1_neutrino_publish_dispatch_failure(
          error, state_transport, state_out, exchange, diagnostics, false,
          neutrino_diagnostics, true);
    preclosure_fallback = ghl_m1_neutrino_closure_fallback_status(&preclosure);
    if(preclosure_fallback)
      return ghl_m1_neutrino_publish_dispatch_failure(
          ghl_error_m1_implicit_solve_failure, state_transport, state_out,
          exchange, diagnostics, true, neutrino_diagnostics, true);
  }

  const bool thin_selected = dt_alpha * rates->kappa_a_E < 1.0 &&
                             dt_alpha * rates->kappa_s < 1.0;
  bool closure_fallback_used = preclosure_fallback;
  ghl_m1_neutrino_source_path_t path =
      ghl_m1_neutrino_source_path_general_implicit;

  if(thin_selected) {
    path = ghl_m1_neutrino_source_path_thin_explicit;
    error = ghl_m1_neutrino_try_thin_branch(
        m1_params, nu_params, metric, prims_frozen, rates,
        state_transport, dt, n_b_cons, selected.thermalized_number_threshold,
        &candidate_neutrino_diagnostics,
        state_out, exchange, &closure_fallback_used);
  } else if(ghl_m1_neutrino_thick_limit_selected(
                 dt_alpha, rates, selected.thick_equilibrium_threshold)) {
    path = ghl_m1_neutrino_source_path_thick_equilibrium;
    error = ghl_m1_neutrino_try_thick_branch(
        m1_params, nu_params, metric, prims_frozen, rates, state_transport,
        dt, n_b_cons, selected.thermalized_number_threshold,
        &candidate_neutrino_diagnostics, state_out, exchange,
        &closure_fallback_used);
  } else if(ghl_m1_neutrino_scattering_limit_selected(
                 dt_alpha, rates, selected.scattering_threshold)) {
    path = ghl_m1_neutrino_source_path_scattering_dominated;
    error = ghl_m1_neutrino_try_scattering_branch(
        m1_params, nu_params, metric, prims_frozen, rates, state_transport,
        dt, n_b_cons, selected.thermalized_number_threshold,
        &candidate_neutrino_diagnostics, state_out, exchange,
        &closure_fallback_used);
  } else {
    ghl_m1_implicit_solve_diagnostics implicit_diagnostics;
    ghl_m1_initialize_implicit_solve_diagnostics(&implicit_diagnostics);
    error = ghl_m1_solve_neutrino_implicit_homogeneous_update(
        m1_params, nu_params, metric, prims_frozen, rates, dt, n_b_cons,
        state_transport, state_out, exchange, &implicit_diagnostics,
        &candidate_neutrino_diagnostics);
    diagnostics->implicit = implicit_diagnostics;
    closure_fallback_used |=
        (implicit_diagnostics.solution_path_flags &
         ghl_m1_solution_path_closure_fallback) != 0u;
    if(error == ghl_error_m1_implicit_terminal_fallback) {
      diagnostics->path = ghl_m1_neutrino_source_path_terminal_no_update;
      diagnostics->terminal_no_update = true;
      *state_out = *state_transport;
      ghl_m1_neutrino_zero_exchange(exchange);
      diagnostics->closure_fallback_used = closure_fallback_used;
      *neutrino_diagnostics = candidate_neutrino_diagnostics;
      return error;
    }
    if(error != ghl_success) {
      /* The legacy solver accounts for its own hard failure through the
       * diagnostics object passed above. Preserve that single increment when
       * the dispatcher converts the result to its transactional failure form. */
      *neutrino_diagnostics = candidate_neutrino_diagnostics;
      return ghl_m1_neutrino_publish_dispatch_failure(
          error, state_transport, state_out, exchange, diagnostics,
          closure_fallback_used, neutrino_diagnostics, false);
    }
  }

  diagnostics->path = path;
  diagnostics->closure_fallback_used = closure_fallback_used;
  diagnostics->terminal_no_update = false;
  if(error != ghl_success)
    return ghl_m1_neutrino_publish_dispatch_failure(
        error, state_transport, state_out, exchange, diagnostics,
        closure_fallback_used, neutrino_diagnostics, true);

  if(closure_fallback_used && !selected.allow_closure_fallback)
    return ghl_m1_neutrino_publish_dispatch_failure(
        ghl_error_m1_implicit_solve_failure, state_transport, state_out,
        exchange, diagnostics, true, neutrino_diagnostics, true);

  error = ghl_m1_neutrino_apply_ye_policy(
      selected.ye_policy, rates, n_b_cons, exchange);
  if(error != ghl_success) {
    *neutrino_diagnostics = candidate_neutrino_diagnostics;
    return ghl_m1_neutrino_publish_dispatch_failure(
        error, state_transport, state_out, exchange, diagnostics,
        closure_fallback_used, neutrino_diagnostics, true);
  }

  if(path != ghl_m1_neutrino_source_path_general_implicit) {
    candidate_neutrino_diagnostics.source_converged++;
    diagnostics->implicit.solution_path_flags =
        ghl_m1_solution_path_primary_convergence |
        ghl_m1_solution_path_endpoint_acceptance |
        (closure_fallback_used ? ghl_m1_solution_path_closure_fallback : 0u);
  }
  *neutrino_diagnostics = candidate_neutrino_diagnostics;
  return ghl_success;
}
