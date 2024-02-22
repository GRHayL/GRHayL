#include "../../utils_Noble.h"

int ghl_general_newton_raphson(
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

  double dx[ndim], x_old[ndim];

  // Initialize various parameters and variables:
  double errx = 1.0;
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
    errx = (x[0]==0.0) ? fabs(dx[0]) : fabs(dx[0]/x[0]);

    harm_aux->n_iter++;

    if( (fabs(errx) <= harm_aux->solver_tolerance) || (harm_aux->n_iter >= harm_aux->max_iterations) ) {
      keep_iterating = 0;
    }
  }   // END of while(keep_iterating)

  /*  Check for bad untrapped divergences : */
  if(!isfinite(f) ||  !isfinite(df)) {
    return 3;
  } else if(fabs(errx) <= harm_aux->solver_tolerance) {
    return 0;
  } else {
    return 2;
  }
}
