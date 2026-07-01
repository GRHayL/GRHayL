#include "../../utils_Noble.h"

ghl_error_codes_t ghl_general_newton_raphson(
      const ghl_eos_parameters *restrict eos,
      harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double indep_var_in,
      double x[],
      void (*validate_x)(
            const harm_aux_vars_struct *restrict harm_aux,
            const double [],
            const double [],
            double []),
      void (*funcd)(
            const ghl_eos_parameters *restrict,
            harm_aux_vars_struct *restrict,
            const double,
            const double [],
            double [],
            double *restrict,
            double *restrict)) {

  double dx[ndim], dx_resid[ndim], x_old[ndim];

  // Initialize various parameters and variables:
  double errx = 1.0;
  double residual = 1.0;
  double df = 1.0;
  double f = 1.0;
  for(int id = 0; id < ndim; id++)
    x_old[id] = x[id];

  /* Start the Newton-Raphson iterations: */
  int keep_iterating = 1;
  while(keep_iterating) {
    funcd(eos, harm_aux, indep_var_in, x, dx, &f, &df);

    /* Save old values and make the newton step: */
    for(int id = 0; id < ndim; id++) {
      x_old[id] = x[id];
      x[id] += dx[id];
    }

    /******************************************/
    /* Make sure that the new x[] is physical */
    /******************************************/
    validate_x(harm_aux, dx, x_old, x);

    /***************************************/
    /* Calculate the convergence criterion */
    /***************************************/
    errx = 0.0;
    for(int id = 0; id < ndim; id++) {
      const double x_scale = fabs(x[id]);
      const double errx_id = (x_scale == 0.0) ? fabs(dx[id]) : fabs(dx[id])/x_scale;
      errx = fmax(errx, errx_id);
    }
    funcd(eos, harm_aux, indep_var_in, x, dx_resid, &f, &df);
    residual = fabs(f);

    harm_aux->n_iter++;

    if( ((errx <= harm_aux->solver_tolerance) && (residual <= harm_aux->solver_tolerance))
        || (harm_aux->n_iter >= harm_aux->max_iterations) ) {
      keep_iterating = 0;
    }
  }   // END of while(keep_iterating)

  /*  Check for bad untrapped divergences : */
  if(!isfinite(f) ||  !isfinite(df)) {
    return ghl_error_c2p_singular;
  } else if(errx <= harm_aux->solver_tolerance && residual <= harm_aux->solver_tolerance) {
    return ghl_success;
  } else {
    return ghl_error_c2p_max_iter;
  }
}
