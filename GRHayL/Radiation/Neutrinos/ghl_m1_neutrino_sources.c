#include "../ghl_m1_utils.h"
#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"
#include <float.h>

static bool scaled_positive_divide_to_scaled(
      const ghl_m1_scaled_positive *restrict numerator,
      const ghl_m1_scaled_positive *restrict denominator,
      ghl_m1_scaled_positive *restrict quotient) {
  if(numerator == NULL || denominator == NULL || quotient == NULL
     || denominator->mantissa == 0.0) {
    return false;
  }
  if(numerator->mantissa == 0.0) {
    *quotient = (ghl_m1_scaled_positive){ .mantissa = 0.0, .exponent = 0 };
    return true;
  }
  int normalization = 0;
  quotient->mantissa
        = frexp(numerator->mantissa / denominator->mantissa, &normalization);
  quotient->exponent = numerator->exponent - denominator->exponent + normalization;
  return true;
}

/* Compute the complete positive BE ratio without materializing any
 * range-limited product or sum. The ordinary input path keeps the established
 * direct expression in ghl_m1_update_neutrino_number_backward_euler. */
static bool scaled_backward_euler_number_endpoint(
      const double dt_alpha,
      const double kappa_a_N,
      const double eta_N,
      const double Gamma_N,
      const double N_in,
      double *restrict N_out) {
  if(N_out == NULL || N_in < 0.0) {
    return false;
  }

  const double emission_factors[2] = { dt_alpha, eta_N };
  const double absorption_factors[2] = { dt_alpha, kappa_a_N };
  ghl_m1_scaled_positive initial_number;
  ghl_m1_scaled_positive emission;
  ghl_m1_scaled_positive numerator;
  ghl_m1_scaled_positive absorption_product;
  ghl_m1_scaled_positive gamma;
  ghl_m1_scaled_positive absorption_term;
  ghl_m1_scaled_positive one;
  ghl_m1_scaled_positive denominator;
  double candidate = 0.0;
  if(!ghl_m1_scaled_positive_from_double(N_in, &initial_number)
     || !ghl_m1_scaled_positive_product(emission_factors, 2, &emission)
     || !ghl_m1_scaled_positive_add(&initial_number, &emission, &numerator)
     || !ghl_m1_scaled_positive_product(absorption_factors, 2, &absorption_product)
     || !ghl_m1_scaled_positive_from_double(Gamma_N, &gamma)
     || !scaled_positive_divide_to_scaled(&absorption_product, &gamma, &absorption_term)
     || !ghl_m1_scaled_positive_from_double(1.0, &one)
     || !ghl_m1_scaled_positive_add(&one, &absorption_term, &denominator)
     || !ghl_m1_scaled_positive_divide(&numerator, &denominator, &candidate)
     || !isfinite(candidate) || candidate < 0.0) {
    return false;
  }
  *N_out = candidate;
  return true;
}

/*
 * Aggregate neutrino E/F_i and N interaction sources plus the shared
 * endpoint-number policy (ordinary endpoint-Gamma_N backward Euler or the
 * optional thermalized mean-energy projection).
 *
 *   Q   = eta_E - kappa_a_E * J
 *   S_E = Q * W + kappa_tr * h_n
 *   S_i = Q * W * V_i - kappa_tr * H_perp_i
 *
 *   N_source = eta_N - kappa_a_N * N / Gamma_N
 *
 *   Backward-Euler N update at the E/F endpoint:
 *     dt_alpha = dt * alpha
 *     N_new = (N_old + dt_alpha * eta_N)
 *             / (1 + dt_alpha * kappa_a_N/Gamma_N)
 *
 * Scattering is number-conserving and comoving-energy-conserving. Pair and
 * thermal electron-flavor reactions require the paired source operation.
 */

static ghl_error_codes_t compute_EF_sources_from_moments(
      const ghl_m1_comoving *restrict comoving,
      const double V_cov[3],
      const double W,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_sources *restrict EF_sources) {

  if(!isfinite(W) || W < 1.0) {
    return ghl_error_u0_singular;
  }

  ghl_m1_sources candidate = { 0 };
  const double Q = rates->eta_E - rates->kappa_a_E * comoving->J;
  candidate.S_E = Q * W + rates->kappa_tr * comoving->Hn;
  for(int i = 0; i < 3; ++i) {
    candidate.S[i] = Q * W * V_cov[i] - rates->kappa_tr * comoving->HD[i];
  }
  if(!isfinite(candidate.S_E)) {
    return ghl_error_m1_invalid_state;
  }
  for(int i = 0; i < 3; ++i) {
    if(!isfinite(candidate.S[i])) {
      return ghl_error_m1_invalid_state;
    }
  }
  *EF_sources = candidate;
  return ghl_success;
}

static ghl_error_codes_t compute_EF_sources_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_sources *restrict EF_sources,
      const bool configuration_validated) {

  ghl_m1_comoving comoving;
  double V_con[3], V_cov[3], W;
  const ghl_error_codes_t error = configuration_validated
                                        ? ghl_m1_compute_comoving_moments_validated(
                                                m1_params, metric, prims, rad_state,
                                                closure, &comoving, V_con, V_cov, &W)
                                        : ghl_m1_compute_comoving_moments_with_velocity(
                                                m1_params, metric, prims, rad_state,
                                                closure, &comoving, V_con, V_cov, &W);
  if(error != ghl_success) {
    return error;
  }
  return compute_EF_sources_from_moments(&comoving, V_cov, W, rates, EF_sources);
}

static ghl_error_codes_t compute_interaction_sources_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_closure *restrict closure,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_sources *restrict EF_sources,
      double *restrict N_source) {

  ghl_m1_comoving comoving;
  double V_con[3], V_cov[3], W;
  ghl_error_codes_t error = ghl_m1_compute_comoving_moments_with_velocity(
        m1_params, metric, prims, rad_state, closure, &comoving, V_con, V_cov, &W);
  if(error != ghl_success) {
    return error;
  }

  ghl_m1_sources candidate_EF = { 0 };
  error = compute_EF_sources_from_moments(&comoving, V_cov, W, rates, &candidate_EF);
  if(error != ghl_success) {
    return error;
  }

  ghl_m1_neutrino_current current;
  error = ghl_m1_neutrino_build_current_from_moments(
        metric, nu_params, state, &comoving, V_con, W, &current);
  if(error != ghl_success) {
    return error;
  }
  double absorption_number = 0.0;
  if(!ghl_m1_neutrino_scaled_absorption_number(
           rates->kappa_a_N, state->N, current.Gamma_N, &absorption_number)) {
    return ghl_error_m1_invalid_state;
  }
  const double candidate_N = rates->eta_N - absorption_number;
  if(!isfinite(candidate_N)) {
    return ghl_error_m1_invalid_state;
  }

  *EF_sources = candidate_EF;
  *N_source = candidate_N;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_neutrino_rates *restrict rates,
      bool *restrict closure_fallback_observed,
      ghl_m1_sources *restrict EF_sources) {

  if(m1_params == NULL || metric == NULL || prims == NULL || rad_state == NULL
     || rates == NULL || EF_sources == NULL) {
    return ghl_error_m1_null_pointer;
  }

  const ghl_error_codes_t rates_error
        = ghl_m1_neutrino_validate_single_species_rates(rates, NULL);
  if(rates_error != ghl_success) {
    return rates_error;
  }

  ghl_m1_closure closure;
  ghl_error_codes_t error = ghl_m1_compute_closure_with_primitives(
        m1_params, metric, prims, rad_state, &closure);
  if(error != ghl_success) {
    return error;
  }
  if(closure_fallback_observed != NULL
     && (closure.solve_status == ghl_m1_closure_solve_endpoint_fallback
         || closure.solve_status == ghl_m1_closure_solve_iteration_exhausted)) {
    *closure_fallback_observed = true;
  }
  return compute_EF_sources_from_closure(
        m1_params, metric, prims, rad_state, &closure, rates, EF_sources, false);
}

ghl_error_codes_t ghl_m1_neutrino_compute_EF_interaction_sources_validated(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_neutrino_rates *restrict rates,
      bool *restrict closure_fallback_observed,
      ghl_m1_sources *restrict EF_sources) {

  if(m1_params == NULL || metric == NULL || prims == NULL || rad_state == NULL
     || rates == NULL || EF_sources == NULL) {
    return ghl_error_m1_null_pointer;
  }

  ghl_m1_closure closure;
  const ghl_error_codes_t error = ghl_m1_compute_closure_minerbo_validated(
        m1_params, metric, prims, rad_state, &closure);
  if(error != ghl_success) {
    return error;
  }
  if(closure_fallback_observed != NULL
     && (closure.solve_status == ghl_m1_closure_solve_endpoint_fallback
         || closure.solve_status == ghl_m1_closure_solve_iteration_exhausted)) {
    *closure_fallback_observed = true;
  }
  return compute_EF_sources_from_closure(
        m1_params, metric, prims, rad_state, &closure, rates, EF_sources, true);
}

ghl_error_codes_t ghl_m1_compute_neutrino_interaction_sources(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_sources *restrict EF_sources,
      double *restrict N_source) {

  if(m1_params == NULL || nu_params == NULL || metric == NULL || prims == NULL
     || state == NULL || rates == NULL || EF_sources == NULL || N_source == NULL) {
    return ghl_error_m1_null_pointer;
  }

  const ghl_error_codes_t rates_error
        = ghl_m1_neutrino_validate_single_species_rates(rates, NULL);
  if(rates_error != ghl_success) {
    return rates_error;
  }
  if(!isfinite(nu_params->N_floor) || nu_params->N_floor < 0.0 || !isfinite(state->N)
     || state->N < nu_params->N_floor) {
    return ghl_error_m1_invalid_state;
  }

  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  ghl_m1_closure computed_closure;
  ghl_error_codes_t error = ghl_m1_compute_closure_with_primitives(
        m1_params, metric, prims, &rad_state, &computed_closure);
  if(error != ghl_success) {
    return error;
  }

  return compute_interaction_sources_from_closure(
        m1_params, nu_params, metric, prims, state, &rad_state, &computed_closure, rates,
        EF_sources, N_source);
}

ghl_error_codes_t ghl_m1_compute_neutrino_interaction_sources_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_sources *restrict EF_sources,
      double *restrict N_source) {

  if(m1_params == NULL || nu_params == NULL || metric == NULL || prims == NULL
     || state == NULL || closure == NULL || rates == NULL || EF_sources == NULL
     || N_source == NULL) {
    return ghl_error_m1_null_pointer;
  }

  const ghl_error_codes_t rates_error
        = ghl_m1_neutrino_validate_single_species_rates(rates, NULL);
  if(rates_error != ghl_success) {
    return rates_error;
  }
  if(!isfinite(nu_params->N_floor) || nu_params->N_floor < 0.0 || !isfinite(state->N)
     || state->N < nu_params->N_floor) {
    return ghl_error_m1_invalid_state;
  }

  const ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(state);
  return compute_interaction_sources_from_closure(
        m1_params, nu_params, metric, prims, state, &rad_state, closure, rates,
        EF_sources, N_source);
}

ghl_error_codes_t ghl_m1_update_neutrino_number_backward_euler(
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt_alpha,
      const double Gamma_N,
      const double N_in,
      double *restrict N_out) {

  if(nu_params == NULL || rates == NULL || N_out == NULL) {
    return ghl_error_m1_null_pointer;
  }

  if(!isfinite(dt_alpha)) {
    return ghl_error_m1_invalid_state;
  }
  if(dt_alpha < 0.0) {
    return ghl_error_m1_invalid_state;
  }
  /* Gamma_N is a current normalization, not a fluid Lorentz factor. */
  const double gamma_floor = nu_params->Gamma_N_floor == 0.0 ? 64.0 * DBL_EPSILON
                                                             : nu_params->Gamma_N_floor;
  if(!isfinite(gamma_floor) || gamma_floor <= 0.0 || !isfinite(Gamma_N)
     || Gamma_N <= gamma_floor) {
    return ghl_error_m1_invalid_state;
  }
  if(!isfinite(N_in)) {
    return ghl_error_m1_invalid_state;
  }
  const ghl_error_codes_t rates_error
        = ghl_m1_neutrino_validate_single_species_rates(rates, NULL);
  if(rates_error != ghl_success) {
    return rates_error;
  }

  /* Endpoint-Gamma_N backward-Euler discretization of
   * dN/dt = eta_N - kappa_a_N*N/Gamma_N. */
  const double absorption_product = dt_alpha * rates->kappa_a_N;
  const double absorption_term = absorption_product / Gamma_N;
  const double denom = 1.0 + absorption_term;
  const double emission = dt_alpha * rates->eta_N;
  const double numer = N_in + emission;
  const double direct_candidate = numer / denom;
  const bool absorption_underflow
        = dt_alpha > 0.0 && rates->kappa_a_N > 0.0 && absorption_term == 0.0;
  const bool emission_underflow
        = dt_alpha > 0.0 && rates->eta_N > 0.0 && emission == 0.0;
  double candidate_N_out = 0.0;
  if(isfinite(absorption_product) && isfinite(absorption_term) && isfinite(denom)
     && denom > 0.0 && isfinite(emission) && isfinite(numer)
     && isfinite(direct_candidate) && !absorption_underflow && !emission_underflow) {
    candidate_N_out = direct_candidate;
  }
  else if(!scaled_backward_euler_number_endpoint(
                dt_alpha, rates->kappa_a_N, rates->eta_N, Gamma_N, N_in,
                &candidate_N_out)) {
    return ghl_error_m1_invalid_state;
  }

  *N_out = candidate_N_out;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_neutrino_update_endpoint_number_with_policy(
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double dt_alpha,
      const double thermalized_number_threshold,
      const ghl_m1_neutrino_state *restrict state_base,
      const ghl_m1_neutrino_current *restrict endpoint_current,
      double *restrict N_out,
      bool *restrict number_projected) {

  if(nu_params == NULL || rates == NULL || state_base == NULL || endpoint_current == NULL
     || N_out == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(!isfinite(dt) || dt < 0.0 || !isfinite(dt_alpha) || dt_alpha < 0.0
     || !isfinite(thermalized_number_threshold)) {
    return ghl_error_m1_invalid_state;
  }

  /* Keep the policy comparison and endpoint formula in one production helper
   * so shortcut and implicit callers cannot silently diverge. */
  const double number_factors[2] = { dt_alpha, rates->kappa_a_N };
  const double threshold_factor[1] = { thermalized_number_threshold };
  if(thermalized_number_threshold >= 0.0
     && ghl_m1_neutrino_scaled_product_meets_threshold(
           number_factors, 2, threshold_factor, 1, true)) {
    double candidate = 0.0;
    if(rates->mean_energy > 0.0) {
      /* Retain the established evaluation order while every intermediate is
       * representable, so normal-range endpoints keep their exact bits.  The
       * scaled fallback also checks finite zero results: with nonnegative
       * inputs, zero is either exact or an underflowed positive quotient.  A
       * genuinely nonrepresentable endpoint still fails below. */
      const double direct
            = endpoint_current->Gamma_N * endpoint_current->J / rates->mean_energy;
      if(isfinite(direct) && direct != 0.0) {
        candidate = direct;
      }
      else {
        const double numerator[2] = { endpoint_current->Gamma_N, endpoint_current->J };
        const double denominator[1] = { rates->mean_energy };
        if(!ghl_m1_neutrino_scaled_ratio_of_products(
                 numerator, 2, denominator, 1, &candidate)) {
          return ghl_error_m1_invalid_state;
        }
      }
    }
    if(!isfinite(candidate)) {
      return ghl_error_m1_invalid_state;
    }
    if(number_projected != NULL) {
      *number_projected = true;
    }
    *N_out = candidate;
    return ghl_success;
  }

  if(number_projected != NULL) {
    *number_projected = false;
  }
  return ghl_m1_update_neutrino_number_backward_euler(
        nu_params, rates, dt_alpha, endpoint_current->Gamma_N, state_base->N, N_out);
}
