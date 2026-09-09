#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"
#include "../ghl_m1_utils.h"
#include <float.h>

ghl_error_codes_t ghl_m1_apply_neutrino_number_floor(
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const double N_in,
      double *restrict N_out,
      bool *restrict floor_applied) {
  if(nu_params == NULL || N_out == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(nu_params->N_floor) || nu_params->N_floor < 0.0 ||
     !isfinite(N_in))
    return ghl_error_m1_invalid_state;
  const bool applied = N_in < nu_params->N_floor;
  *N_out = applied ? nu_params->N_floor : N_in;
  if(floor_applied != NULL)
    *floor_applied = applied;
  return ghl_success;
}

void ghl_m1_neutrino_diagnostics_initialize(
      ghl_m1_neutrino_diagnostics *restrict diagnostics) {
  if(diagnostics != NULL)
    *diagnostics = (ghl_m1_neutrino_diagnostics){0};
}

static void ghl_m1_accumulate_neutrino_repair_magnitudes(
      ghl_m1_neutrino_diagnostics *restrict diagnostics,
      const ghl_m1_neutrino_state *restrict state_before,
      const ghl_m1_neutrino_state *restrict state_after) {
  if(diagnostics == NULL || state_before == NULL || state_after == NULL)
    return;

  /* These fields are magnitudes, not signed conservation deltas. Keep the
   * accumulation in one helper so positive and negative repairs cannot
   * cancel across calls. */
  diagnostics->repair_dN += fabs(state_after->N - state_before->N);
  diagnostics->repair_dE += fabs(state_after->E - state_before->E);
  for(int i = 0; i < 3; i++)
    diagnostics->repair_dF[i] +=
        fabs(state_after->F[i] - state_before->F[i]);

  const double dL_e = diagnostics->repair_lepton_weight *
      (state_after->N - state_before->N);
  diagnostics->repair_dL_e += fabs(dL_e);
}

ghl_error_codes_t ghl_m1_neutrino_check_EN_bounds(
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_current *restrict current) {
  if(state == NULL || nu_params == NULL || current == NULL)
    return ghl_error_m1_null_pointer;
  if(!nu_params->enforce_mean_energy_bounds)
    return ghl_success;

  if(!isfinite(nu_params->mean_energy_min) ||
     !isfinite(nu_params->mean_energy_max))
    return ghl_error_m1_invalid_state;
  if(nu_params->mean_energy_min > 0.0 &&
     nu_params->mean_energy_max > 0.0 &&
     nu_params->mean_energy_min > nu_params->mean_energy_max)
    return ghl_error_m1_invalid_state;
  if(!isfinite(state->N) || state->N < nu_params->N_floor)
    return ghl_error_m1_invalid_state;
  if(state->N == 0.0)
    return ghl_success;
  if(!isfinite(current->J) || !isfinite(current->Gamma_N))
    return ghl_error_m1_invalid_state;

  const double mean_energy = current->J * current->Gamma_N / state->N;
  if(!isfinite(mean_energy))
    return ghl_error_m1_invalid_state;
  if(nu_params->mean_energy_min > 0.0 &&
     mean_energy < nu_params->mean_energy_min)
    return ghl_error_m1_invalid_state;
  if(nu_params->mean_energy_max > 0.0 &&
     mean_energy > nu_params->mean_energy_max)
    return ghl_error_m1_invalid_state;
  return ghl_success;
}

/*
 * Repair a neutrino state transactionally: apply the N floor, delegate E/F_i
 * repair to ghl_m1_realizability_repair, and publish diagnostics only after
 * the candidate state has been repaired successfully.
 */

ghl_error_codes_t ghl_m1_repair_neutrino_state(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      ghl_m1_neutrino_state *restrict state,
      ghl_m1_neutrino_diagnostics *restrict diagnostics) {

  if(m1_params == NULL || nu_params == NULL ||
     metric == NULL  || state == NULL)
    return ghl_error_m1_null_pointer;

  /* Preserve the public error order without performing a throwaway repair. */
  const ghl_error_codes_t configuration_error =
      ghl_m1_validate_configuration(m1_params, metric);
  if(configuration_error != ghl_success)
    return configuration_error;

  /* N_floor must be finite and nonnegative. */
  if(!isfinite(nu_params->N_floor) || nu_params->N_floor < 0.0)
    return ghl_error_m1_invalid_state;

  /* Reject nonfinite raw input before any repair is applied. */
  if(!isfinite(state->N) || !isfinite(state->E))
    return ghl_error_m1_invalid_state;
  for(int i = 0; i < 3; i++) {
    if(!isfinite(state->F[i]))
      return ghl_error_m1_invalid_state;
  }

  /* Repair is an explicit host-elected operation, but its publication is
   * transactional: a delegated E/F repair failure leaves the caller state and
   * every diagnostic counter exactly unchanged. */
  ghl_m1_neutrino_state candidate = *state;
  ghl_m1_neutrino_diagnostics candidate_diagnostics;
  ghl_m1_neutrino_diagnostics *candidate_diag = NULL;
  if(diagnostics != NULL) {
    candidate_diagnostics = *diagnostics;
    candidate_diag = &candidate_diagnostics;
  }

  /* Snapshot the complete state so all repair budgets are accounted for once,
   * after the transactional repair has succeeded. */
  const ghl_m1_neutrino_state state_pre_repair = candidate;

  /* N floor */
  bool N_floor_applied;
  double N_floored;
  ghl_error_codes_t error = ghl_m1_apply_neutrino_number_floor(
      nu_params, candidate.N, &N_floored, &N_floor_applied);
  if(error != ghl_success)
    return error;
  candidate.N = N_floored;
  if(N_floor_applied && candidate_diag != NULL) {
    candidate_diag->N_floor_repairs++;
  }

  /* E/F_i repair delegated to the shared M1 realizability repair. */
  ghl_m1_rad_state rad_state = ghl_m1_neutrino_project_rad_state(&candidate);

  error = ghl_m1_realizability_repair(m1_params, metric, &rad_state);
  if(error != ghl_success)
    return error;

  /* Inputs are finite. Every changed binary64 value is a published repair,
   * including a small floor/cone correction near zero. */
  candidate.E = rad_state.E;
  for(int i = 0; i < 3; i++)
    candidate.F[i] = rad_state.F[i];

  bool ef_mutated = candidate.E != state_pre_repair.E;
  for(int i = 0; i < 3; i++)
    ef_mutated |= candidate.F[i] != state_pre_repair.F[i];

  if(ef_mutated && candidate_diag != NULL)
    candidate_diag->EF_repairs++;

  if(candidate_diag != NULL)
    ghl_m1_accumulate_neutrino_repair_magnitudes(
        candidate_diag, &state_pre_repair, &candidate);

  *state = candidate;
  if(diagnostics != NULL)
    *diagnostics = candidate_diagnostics;

  return ghl_success;
}
