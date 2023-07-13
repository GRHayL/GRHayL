#include "../../harm_u2p_util.h"

int ghl_general_newton_raphson(
      const ghl_eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double indep_var_in,
      int *restrict n_iter_ptr,
      double x[],
      void (*funcd)(const ghl_eos_parameters *restrict, const harm_aux_vars_struct *restrict, const int, const double, const double [], double [],
                    double [], double [][ndim], double *restrict, double *restrict, int *restrict))
{
  double dx[ndim], x_old[ndim];
  double resid[ndim], jac[ndim][ndim];

  // Initialize various parameters and variables:
  int n_iter = 0;
  double errx = 1.0;
  double df = 1.0;
  double f = 1.0;
  for( int id = 0; id < ndim; id++)  x_old[id] = x[id];

  /* Start the Newton-Raphson iterations : */
  int keep_iterating = 1;
  while( keep_iterating ) {
    funcd(eos, harm_aux, ndim, indep_var_in, x, dx, resid, jac, &f, &df, n_iter_ptr);

    /* Save old values before calculating the new: */
    errx = 0.0;
    for( int id = 0; id < ndim; id++) {
      x_old[id] = x[id];
    }

    /* Make the newton step: */
    for( int id = 0; id < ndim; id++) {
      x[id] += dx[id];
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx = (x[0]==0.0) ? fabs(dx[0]) : fabs(dx[0]/x[0]);

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    ghl_validate_x(x, x_old);

    if( ((fabs(errx) <= NEWT_TOL)) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)

  *n_iter_ptr = n_iter;

  /*  Check for bad untrapped divergences : */
  if( !isfinite(f) ||  !isfinite(df) ) {
    return 3;
  } else if( fabs(errx) <= NEWT_TOL ) {
    return 0;
  } else {
    return 2;
  }
}
