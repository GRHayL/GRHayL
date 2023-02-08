#include "../harm_u2p_util.h"

#define NEWT_DIM (2)

void validate_x(double x[2], const double x0[2] );

int general_newton_raphson(
      const eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      double x[],
      const int n,
      int *restrict n_iter_ptr,

      void (*funcd)(const eos_parameters *restrict, const harm_aux_vars_struct *restrict, const double [], double [], double [],
                    double [][NEWT_DIM], double *, double *, int) )
{
  double dx[n], x_old[n];
  double resid[n], jac[n][n];
  double x_orig[n];

  // Initialize various parameters and variables:
  int n_iter = 0;
  int i_extra = 0;
  int doing_extra = 0;
  double errx = 1.0;
  double df = 1.0;
  double f = 1.0;
  for( int id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  /* Start the Newton-Raphson iterations : */
  int keep_iterating = 1;
  while( keep_iterating ) {
    (*funcd) (eos, harm_aux, x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */


    /* Save old values before calculating the new: */
    errx = 0.;
    for( int id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    for( int id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.0) ? fabs(dx[0]) : fabs(dx[0]/x[0]);

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    validate_x( x, x_old ) ;

    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0))
        || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)

  *n_iter_ptr = n_iter;

  /*  Check for bad untrapped divergences : */
//TODO: error-checking
//  if( (robust_isfinite(f)==0) ||  (robust_isfinite(df)==0) ) {
//    return(2);
//  }

  if( fabs(errx) <= NEWT_TOL ) {
    return(0);
  } else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ) {
    return(0);
  } else {
    return(1);
  }

  return(0);

}

void validate_x(double x[2], const double x0[2] ) {
  const double dv = 1.e-15;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */
  x[0] = fabs(x[0]);
  x[0] = (x[0] > W_TOO_BIG) ?  x0[0] : x[0];


  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] > 1.) ?  (1. - dv) : x[1];  /* if it's too big   */

  return;
}
