#include "ghl_m1.h"
#include "ghl_m1_neutrino_implicit.h"
#include "../ghl_m1_utils.h"

/*
 * Finite-difference Jacobian for the neutrino implicit E/F_i residual.
 * The step-size and one-sided fallback policies are shared by the M1 Jacobian.
 *
 * Uses the shared relative/absolute step-size policy
 * (m1_params->fd_epsilon_rel, m1_params->fd_epsilon_abs) and preserves the
 * one-sided fallback behavior for admissibility failures. The Jacobian is 4x4
 * and matches the residual ordering (tildeE, tildeF_i).
 *
 * Because the neutrino residual uses frozen rates and frozen primitives (no
 * Con2Prim, no opacity callback), the finite-difference perturbation only
 * changes the trial rad_state, closure, comoving moments, and source terms.
 * Hard failures (non-admissibility errors that are not
 * ghl_error_m1_implicit_admissibility) are propagated directly instead of
 * being re-labeled as Jacobian failures.
 */

static ghl_error_codes_t ghl_m1_neutrino_compute_implicit_jacobian_core(
      const ghl_m1_neutrino_implicit_context *restrict context,
      const double dt,
      const double U_base[4],
      const double U[4],
      const double residual_0[4],
      const bool use_checked_residual,
      double jacobian[4][4]) {

  if(context == NULL || context->m1_params == NULL || context->metric == NULL ||
     context->prims_frozen == NULL || context->rates == NULL ||
     U_base == NULL || U == NULL || residual_0 == NULL || jacobian == NULL)
    return ghl_error_m1_null_pointer;

  const ghl_m1_parameters *restrict m1_params = context->m1_params;
  const ghl_metric_quantities *restrict metric = context->metric;

  for(int i = 0; i < 4; i++) {
    if(!isfinite(U_base[i]) || !isfinite(U[i]) || !isfinite(residual_0[i]))
      return ghl_error_m1_invalid_state;
  }

  for(int n = 0; n < 4; n++) {
    const double delta = ghl_m1_compute_fd_delta(m1_params, metric, U_base[n], U[n]);
    if(!isfinite(delta) || delta <= 0.0) {
      return ghl_error_m1_invalid_state;
    }

    double U_perturbed[4] = { U[0], U[1], U[2], U[3] };
    U_perturbed[n] = U[n] + delta;

    double residual_perturbed[4] = { 0.0, 0.0, 0.0, 0.0 };
    ghl_error_codes_t fd_error = use_checked_residual
        ? ghl_m1_neutrino_compute_implicit_residual_with_base(
              m1_params, metric, context->prims_frozen, context->rates, dt,
              U_base, U_perturbed, residual_perturbed)
        : ghl_m1_neutrino_compute_implicit_residual_validated(
              context, dt, U_base, U_perturbed, NULL, residual_perturbed);

    double used_delta = delta;
    if(fd_error != ghl_success) {
      if(!ghl_m1_fd_error_allows_one_sided_fallback(fd_error)) {
        return fd_error;
      }

      /* One-sided backward difference when the forward perturbation leaves the
       * admissible domain. */
      U_perturbed[n] = U[n] - delta;
      fd_error = use_checked_residual
          ? ghl_m1_neutrino_compute_implicit_residual_with_base(
                m1_params, metric, context->prims_frozen, context->rates, dt,
                U_base, U_perturbed, residual_perturbed)
          : ghl_m1_neutrino_compute_implicit_residual_validated(
                context, dt, U_base, U_perturbed, NULL, residual_perturbed);
      used_delta = -delta;
      if(fd_error != ghl_success) {
        if(ghl_m1_fd_error_allows_one_sided_fallback(fd_error)) {
          return ghl_error_m1_invalid_implicit_jacobian;
        }
        return fd_error;
      }
    }

    for(int i = 0; i < 4; i++) {
      jacobian[i][n] = (residual_perturbed[i] - residual_0[i]) / used_delta;
      if(!isfinite(jacobian[i][n]))
        return ghl_error_m1_invalid_state;
    }
  }

  return ghl_success;
}

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_jacobian_with_base(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const double dt,
      const double U_base[4],
      const double U[4],
      const double residual_0[4],
      double jacobian[4][4]) {

  if(m1_params == NULL || metric == NULL ||
     prims_frozen == NULL || rates == NULL)
    return ghl_error_m1_null_pointer;
  const ghl_m1_neutrino_implicit_context context = {
        .m1_params = m1_params,
        .metric = metric,
        .prims_frozen = prims_frozen,
        .rates = rates };
  return ghl_m1_neutrino_compute_implicit_jacobian_core(
      &context, dt, U_base, U, residual_0, true, jacobian);
}

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_jacobian_validated(
      const ghl_m1_neutrino_implicit_context *restrict context,
      const double dt,
      const double U_base[4],
      const double U[4],
      const double residual_0[4],
      double jacobian[4][4]) {
  return ghl_m1_neutrino_compute_implicit_jacobian_core(
      context, dt, U_base, U, residual_0, false, jacobian);
}

ghl_error_codes_t ghl_m1_neutrino_compute_implicit_jacobian(
      const ghl_m1_parameters *restrict m1_params,
      const ghl_m1_neutrino_parameters *restrict nu_params,
      const ghl_metric_quantities *restrict metric,
      const ghl_primitive_quantities *restrict prims_frozen,
      const ghl_m1_neutrino_rates *restrict rates,
      const ghl_m1_neutrino_state *restrict state_in,
      const double dt,
      const double U[4],
      const double residual_0[4],
      double jacobian[4][4]) {

  if(metric == NULL || state_in == NULL)
    return ghl_error_m1_null_pointer;

  if(!ghl_m1_metric_is_symmetric_spd(metric))
    return ghl_error_m1_invalid_metric;

  if(nu_params == NULL)
    return ghl_error_m1_null_pointer;

  const double sqrt_detgamma = metric->sqrt_detgamma;
  const double U_base[4] = {
    state_in->E * sqrt_detgamma,
    state_in->F[0] * sqrt_detgamma,
    state_in->F[1] * sqrt_detgamma,
    state_in->F[2] * sqrt_detgamma
  };

  return ghl_m1_neutrino_compute_implicit_jacobian_with_base(
      m1_params, metric, prims_frozen, rates, dt, U_base, U, residual_0,
      jacobian);
}
