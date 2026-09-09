#ifndef GHL_M1_NEUTRINO_IMPLICIT_H
#define GHL_M1_NEUTRINO_IMPLICIT_H

/*
 * Generic private Radiation header for grey three-species neutrino M1
 * number-current, interaction-source, and implicit-solve helpers.
 *
 * This header declares helpers used by
 * ghl_m1_solve_neutrino_implicit_homogeneous_update across translation units.
 * Public residual, Jacobian, trial-state, and admissibility declarations live
 * in ghl_m1.h; the explicit-base variants and remaining helpers here remain
 * internal:
 *
 *   - The neutrino residual and Jacobian use frozen rates and frozen
 *     primitives.
 *   - The neutrino frozen-rate solve uses the shared Newton driver but keeps
 *     its rate residual, Jacobian, and solver helpers in this subdirectory.
 *
 * Its scope includes the shared number-current helper, interaction sources,
 * and the local implicit solve. It is listed in the neutrino build manifest's
 * private include line (#! INCS =).
 */

#include "ghl_m1.h"

/* Public single-species source boundaries additionally require that electron
 * rates contain no reaction requiring the missing partner state. */
ghl_error_codes_t ghl_m1_neutrino_validate_single_species_rates(
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_neutrino_diagnostics *restrict diagnostics);

typedef struct {
  double J;
  double h_n;
  double HU[3];
  double Gamma_N;
  double n_com;
  double number_flux[3];
  double number_transport_velocity[3];
} ghl_m1_neutrino_current;

/* Immutable inputs shared by all Newton residual and Jacobian evaluations in
 * one neutrino solve. The solve entry point validates this context before it
 * is handed to the private callbacks. */
typedef struct {
  const ghl_m1_parameters *restrict m1_params;
  const ghl_metric_quantities *restrict metric;
  const ghl_primitive_quantities *restrict prims_frozen;
  const ghl_m1_neutrino_rates *restrict rates;
} ghl_m1_neutrino_implicit_context;

/* Build the number current from moments and velocity quantities that have
 * already been evaluated for the same radiation state. */
ghl_error_codes_t ghl_m1_neutrino_build_current_from_moments(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_comoving *restrict comoving,
      const double V_con[3],
      const double W,
      ghl_m1_neutrino_current *restrict current);

/* Radiation-private transactional exchange assembly. */
ghl_error_codes_t ghl_m1_neutrino_assemble_exchange(
      const ghl_m1_neutrino_state *restrict state_in,
      const ghl_m1_neutrino_state *restrict state_out,
      const ghl_m1_neutrino_rates *restrict rates,
      double dL_rad_cc,
      double sqrt_detgamma,
      double baryon_density_conserved,
      ghl_m1_neutrino_exchange *restrict exchange);

ghl_error_codes_t ghl_m1_neutrino_derive_current(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      ghl_m1_neutrino_current *restrict current);

/* Private final-endpoint validation shared by all local source routes. */
ghl_error_codes_t ghl_m1_neutrino_check_EN_bounds(
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_current *restrict current);

/* Closure-supplied localization layer. Production callers use
 * ghl_m1_neutrino_derive_current, which first evaluates the configured
 * primitive-aware closure. Direct verifiers may use this boundary to compare
 * the downstream current algebra from one independently generated tensor. */
ghl_error_codes_t ghl_m1_neutrino_derive_current_from_closure(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_closure *restrict closure,
      ghl_m1_neutrino_current *restrict current);

ghl_error_codes_t ghl_m1_neutrino_physical_number_flux_from_current(
      const ghl_metric_quantities *restrict metric,
      const ghl_m1_neutrino_state *restrict state,
      const ghl_m1_neutrino_current *restrict current,
      const ghl_m1_direction_t direction,
      double *restrict physical_number_flux);

/* Solver-observability variant. This private entry point only reports whether
 * the selected closure used a finite nonordinary fallback candidate. */
ghl_error_codes_t ghl_m1_neutrino_compute_EF_interaction_sources_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_neutrino_rates *restrict rates,
      bool *restrict closure_fallback_observed,
      ghl_m1_sources *restrict EF_sources);

/* Internal solve path. The metric and immutable M1 configuration are checked
 * once by the solve entry point; this variant retains trial-state and closure
 * output checks without repeating that boundary validation. */
ghl_error_codes_t ghl_m1_neutrino_compute_EF_interaction_sources_validated(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims,
      const ghl_m1_rad_state *restrict rad_state,
      const ghl_m1_neutrino_rates *restrict rates,
      bool *restrict closure_fallback_observed,
      ghl_m1_sources *restrict EF_sources);

/**
 * Compute the same residual using an explicit densitized substep base U_base.
 *
 * This is used by the fallback substepping driver so substep k solves
 * U_{k+1} - U_k - dt_sub*S(U_{k+1}) = 0 instead of repeatedly subtracting
 * the original full-step input state. The public-like wrapper above passes
 * U_base derived from state_in for single-step tests.
 */
ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_with_base(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double U_base[4],
      const double U[4],
      double residual[4]);

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_with_base_diagnostics(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double U_base[4],
      const double U[4],
      bool *restrict closure_fallback_observed,
      double residual[4]);

/* Private callback core. The context's metric/configuration/rates have
 * already crossed their public validation boundary. */
ghl_error_codes_t ghl_m1_neutrino_compute_implicit_residual_validated(
      const ghl_m1_neutrino_implicit_context *restrict context,
      const double dt,
      const double U_base[4],
      const double U[4],
      bool *restrict closure_fallback_observed,
      double residual[4]);

/**
 * Finite-difference Jacobian for the explicit-base residual variant.
 */
ghl_error_codes_t ghl_m1_neutrino_compute_implicit_jacobian_with_base(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double U_base[4],
      const double U[4],
      const double residual_0[4],
      double jacobian[4][4]);

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_jacobian_validated(
      const ghl_m1_neutrino_implicit_context *restrict context,
      const double dt,
      const double U_base[4],
      const double U[4],
      const double residual_0[4],
      double jacobian[4][4]);

/* Run one validated E/F Newton substep. Pair-source code supplies an
 * internally constructed effective rate bundle after validating the original
 * provider rates; this bridge deliberately does not repeat public rate
 * validation. */
ghl_error_codes_t ghl_m1_neutrino_attempt_EF_newton_step(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt_sub,
      const double U_in[4],
      double U_out[4],
      ghl_m1_newton_diagnostics *restrict diagnostics,
      bool *restrict closure_fallback_observed);

/* Populate endpoint mean-energy observability fields after a paired solve. */
void ghl_m1_neutrino_populate_mean_energy_diagnostics(
      const ghl_m1_neutrino_state *restrict state_out,
      const ghl_m1_neutrino_current *restrict current,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_m1_neutrino_rates *restrict rates,
      ghl_m1_neutrino_diagnostics *restrict neutrino_diagnostics);

#endif // GHL_M1_NEUTRINO_IMPLICIT_H
