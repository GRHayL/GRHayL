#include "../../utils_Noble.h"

/**
 * @ingroup c2p_internal
 * @brief Runs Newton iteration with a separate budget for each invocation.
 * @details n_iter accumulates attempted iterations, including nested solves
 * and retries. Residual callbacks may set harm_aux->error to propagate a
 * failure before their outputs are used. Nonfinite outputs or iterates return
 * ghl_error_c2p_singular. A nonpositive budget performs no iterations.
 */
ghl_error_codes_t ghl_general_newton_raphson(
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double indep_var_in,
      double x[],
      void (*validate_x)(
            const harm_aux_vars_struct *restrict harm_aux,
            const double[],
            const double[],
            double[]),
      void (*funcd)(
            const ghl_eos_parameters *restrict,
            harm_aux_vars_struct *restrict,
            const double,
            const double[],
            double[],
            double *restrict,
            double *restrict)) {

  double dx[ndim], dx_resid[ndim], x_old[ndim];

  harm_aux->error = ghl_success;
  if(!isfinite(indep_var_in)) {
    return ghl_error_c2p_singular;
  }
  for(int id = 0; id < ndim; id++) {
    if(!isfinite(x[id])) {
      return ghl_error_c2p_singular;
    }
  }

  /* Start the Newton-Raphson iterations: */
  for(int n_iter = 0; n_iter < harm_aux->max_iterations; n_iter++) {
    harm_aux->n_iter++;
    double f, df;
    funcd(eos, harm_aux, indep_var_in, x, dx, &f, &df);
    if(harm_aux->error != ghl_success) {
      return harm_aux->error;
    }
    if(!isfinite(f) || !isfinite(df)) {
      return ghl_error_c2p_singular;
    }

    /* Save old values and make the newton step: */
    for(int id = 0; id < ndim; id++) {
      if(!isfinite(dx[id])) {
        return ghl_error_c2p_singular;
      }
      x_old[id] = x[id];
      x[id] += dx[id];
      if(!isfinite(x[id])) {
        return ghl_error_c2p_singular;
      }
    }

    /******************************************/
    /* Make sure that the new x[] is physical */
    /******************************************/
    validate_x(harm_aux, dx, x_old, x);

    /***************************************/
    /* Calculate the convergence criterion */
    /***************************************/
    double errx = 0.0;
    for(int id = 0; id < ndim; id++) {
      if(!isfinite(x[id])) {
        return ghl_error_c2p_singular;
      }
      // Only hybrid Noble2D's bounded v^2 uses the applied step. Z, density,
      // and tabulated recovery need the raw correction: validation can hide a
      // large Newton update behind a small absolute residual.
      const double abs_dx
            = fabs((funcd == ghl_func_2D && id == 1) ? x[id] - x_old[id] : dx[id]);
      const double x_scale = fabs(x[id]);
      const double errx_id = (x_scale == 0.0) ? abs_dx : abs_dx / x_scale;
      errx = fmax(errx, errx_id);
    }
    funcd(eos, harm_aux, indep_var_in, x, dx_resid, &f, &df);
    if(harm_aux->error != ghl_success) {
      return harm_aux->error;
    }
    if(!isfinite(f) || !isfinite(df)) {
      return ghl_error_c2p_singular;
    }
    for(int id = 0; id < ndim; id++) {
      if(!isfinite(dx_resid[id])) {
        return ghl_error_c2p_singular;
      }
    }
    if(errx <= harm_aux->solver_tolerance && fabs(f) <= harm_aux->solver_tolerance) {
      return ghl_success;
    }
  }

  return ghl_error_c2p_max_iter;
}
