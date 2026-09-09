#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"
#include "../ghl_m1_utils.h"
#include <float.h>

ghl_error_codes_t ghl_m1_neutrino_build_current_from_moments(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_comoving *restrict comoving,
      const double V_con[3],
      const double W,
      ghl_m1_neutrino_current *restrict current) {

  if(metric == NULL || nu_params == NULL || state == NULL ||
     comoving == NULL || V_con == NULL || current == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(nu_params->N_floor) || nu_params->N_floor < 0.0 ||
     !isfinite(nu_params->J_floor) || nu_params->J_floor < 0.0 ||
     !isfinite(nu_params->Gamma_N_floor) || nu_params->Gamma_N_floor < 0.0 ||
     !isfinite(state->N) || state->N < nu_params->N_floor)
    return ghl_error_m1_invalid_state;

  if(!isfinite(W) || W < 1.0 ||
     !isfinite(comoving->J) || !(comoving->J > nu_params->J_floor) ||
     !isfinite(comoving->Hn))
    return ghl_error_m1_invalid_state;
  const double gamma_floor = nu_params->Gamma_N_floor == 0.0
                           ? 64.0 * DBL_EPSILON : nu_params->Gamma_N_floor;
  const double Gamma_N = W - comoving->Hn / comoving->J;
  ghl_m1_neutrino_current candidate = {0};
  candidate.J = comoving->J;
  candidate.h_n = comoving->Hn;
  for(int i = 0; i < 3; i++) {
    if(!isfinite(comoving->HU[i]))
      return ghl_error_m1_invalid_state;
    candidate.HU[i] = comoving->HU[i];
  }

  /* N=0 carries no number current.  Gamma_N is only needed to define the
   * velocity of a nonzero number current; in the zero-density limit that
   * velocity is undefined and a singular Gamma_N must not reject an otherwise
   * valid E/F state.  Keep the nonzero-N path strict and fail closed. */
  if(!isfinite(Gamma_N) || !(Gamma_N > gamma_floor)) {
    if(state->N != 0.0)
      return ghl_error_m1_invalid_state;
    /* W is a finite positive normalization for the zero-density limit.  It
     * keeps downstream source formulas such as N/Gamma_N defined while the
     * physical number current remains exactly zero. */
    candidate.Gamma_N = W;
    candidate.n_com = 0.0;
    for(int i = 0; i < 3; i++) {
      candidate.number_transport_velocity[i] = 0.0;
      candidate.number_flux[i] = 0.0;
    }
    *current = candidate;
    return ghl_success;
  }

  const double n_com = state->N / Gamma_N;
  if(!isfinite(n_com))
    return ghl_error_m1_invalid_state;
  candidate.Gamma_N = Gamma_N;
  candidate.n_com = n_com;
  for(int i = 0; i < 3; i++) {
    candidate.number_transport_velocity[i] =
        (W * V_con[i] + comoving->HU[i] / comoving->J) / Gamma_N;
    candidate.number_flux[i] = state->N * candidate.number_transport_velocity[i];
    if(!isfinite(candidate.number_transport_velocity[i]) ||
       !isfinite(candidate.number_flux[i]))
      return ghl_error_m1_invalid_state;
  }
  const ghl_error_codes_t error = ghl_m1_validate_transport_velocity(
      metric, candidate.number_transport_velocity);
  if(error != ghl_success)
    return error;
  *current = candidate;
  return ghl_success;
}

static ghl_error_codes_t derive_current(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure_in,
      ghl_m1_neutrino_current *restrict current) {

  if(m1_params == NULL || nu_params == NULL || metric == NULL || prims == NULL ||
     state == NULL || current == NULL)
    return ghl_error_m1_null_pointer;
  if(!ghl_m1_metric_is_symmetric_spd(metric))
    return ghl_error_m1_invalid_metric;
  if(!isfinite(nu_params->N_floor) || nu_params->N_floor < 0.0 ||
     !isfinite(nu_params->J_floor) || nu_params->J_floor < 0.0 ||
     !isfinite(nu_params->Gamma_N_floor) || nu_params->Gamma_N_floor < 0.0 ||
     !isfinite(state->N) || state->N < nu_params->N_floor)
    return ghl_error_m1_invalid_state;

  const ghl_m1_rad_state rad = ghl_m1_neutrino_project_rad_state(state);
  ghl_error_codes_t error = ghl_m1_validate_realizability(
      m1_params, metric, &rad, 0.0, NULL);
  if(error != ghl_success)
    return error;

  ghl_m1_closure computed_closure;
  const ghl_m1_closure *closure = closure_in;
  if(closure == NULL) {
    error = ghl_m1_compute_closure_with_primitives(
        m1_params, metric, prims, &rad, &computed_closure);
    if(error != ghl_success)
      return error;
    closure = &computed_closure;
  }

  ghl_m1_comoving comoving;
  double VU[3], VD[3], W;
  error = ghl_m1_compute_comoving_moments_with_velocity(
      m1_params, metric, prims, &rad, closure, &comoving, VU, VD, &W);
  if(error != ghl_success)
    return error;
  return ghl_m1_neutrino_build_current_from_moments(
      metric, nu_params, state, &comoving, VU, W, current);
}

ghl_error_codes_t ghl_m1_neutrino_derive_current(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      ghl_m1_neutrino_current *restrict current) {
  return derive_current(
      m1_params, nu_params, metric, prims, state, NULL, current);
}

ghl_error_codes_t ghl_m1_neutrino_derive_current_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_neutrino_current *restrict current) {
  if(closure == NULL || current == NULL)
    return ghl_error_m1_null_pointer;
  return derive_current(
      m1_params, nu_params, metric, prims, state, closure, current);
}

ghl_error_codes_t ghl_m1_neutrino_physical_number_flux_from_current(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_neutrino_current *restrict current,
      const ghl_m1_direction_t direction,
      double *restrict physical_number_flux) {
  if(metric == NULL || state == NULL || current == NULL ||
     physical_number_flux == NULL)
    return ghl_error_m1_null_pointer;
  const ghl_error_codes_t direction_error =
      ghl_m1_validate_direction(direction);
  if(direction_error != ghl_success)
    return direction_error;
  if(!ghl_m1_metric_is_symmetric_spd(metric))
    return ghl_error_m1_invalid_metric;
  if(!isfinite(state->N))
    return ghl_error_m1_invalid_state;
  for(int i = 0; i < 3; ++i) {
    if(!isfinite(current->number_flux[i]) ||
       !isfinite(current->number_transport_velocity[i]))
      return ghl_error_m1_invalid_state;
  }
  const ghl_error_codes_t velocity_error = ghl_m1_validate_transport_velocity(
      metric, current->number_transport_velocity);
  if(velocity_error != ghl_success)
    return velocity_error;
  const double candidate = metric->lapse * current->number_flux[direction]
                         - metric->betaU[direction] * state->N;
  if(!isfinite(candidate))
    return ghl_error_m1_invalid_state;
  *physical_number_flux = candidate;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_neutrino_number_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      double number_flux[3],
      double number_transport_velocity[3]) {
  if(number_flux == NULL || number_transport_velocity == NULL)
    return ghl_error_m1_null_pointer;
  ghl_m1_neutrino_current current;
  const ghl_error_codes_t error = ghl_m1_neutrino_derive_current(
      m1_params, nu_params, metric, prims, state, &current);
  if(error != ghl_success)
    return error;
  for(int i = 0; i < 3; i++) {
    number_flux[i] = current.number_flux[i];
    number_transport_velocity[i] = current.number_transport_velocity[i];
  }
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_neutrino_number_flux_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      double number_flux[3],
      double number_transport_velocity[3]) {
  if(closure == NULL || number_flux == NULL ||
     number_transport_velocity == NULL)
    return ghl_error_m1_null_pointer;
  ghl_m1_neutrino_current current;
  const ghl_error_codes_t error = ghl_m1_neutrino_derive_current_from_closure(
      m1_params, nu_params, metric, prims, state, closure, &current);
  if(error != ghl_success)
    return error;
  for(int i = 0; i < 3; i++) {
    number_flux[i] = current.number_flux[i];
    number_transport_velocity[i] = current.number_transport_velocity[i];
  }
  return ghl_success;
}

ghl_error_codes_t ghl_m1_compute_neutrino_physical_number_flux(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_direction_t direction,
      double *restrict physical_number_flux) {
  if(metric == NULL || state == NULL || physical_number_flux == NULL)
    return ghl_error_m1_null_pointer;
  ghl_m1_neutrino_current current;
  const ghl_error_codes_t error = ghl_m1_neutrino_derive_current(
      m1_params, nu_params, metric, prims, state, &current);
  if(error != ghl_success)
    return error;
  return ghl_m1_neutrino_physical_number_flux_from_current(
      metric, state, &current, direction, physical_number_flux);
}

ghl_error_codes_t ghl_m1_compute_neutrino_physical_number_flux_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      const ghl_m1_direction_t direction,
      double *restrict physical_number_flux) {
  if(closure == NULL || physical_number_flux == NULL)
    return ghl_error_m1_null_pointer;
  ghl_m1_neutrino_current current;
  const ghl_error_codes_t error = ghl_m1_neutrino_derive_current_from_closure(
      m1_params, nu_params, metric, prims, state, closure, &current);
  if(error != ghl_success)
    return error;
  return ghl_m1_neutrino_physical_number_flux_from_current(
      metric, state, &current, direction, physical_number_flux);
}
