#include "../../../utils_Noble.h"

void ghl_validate_1D(
      const double x0[1],
      double x[1]) {
  x[0] = fabs(x[0]);
}

void ghl_func_1D(
      const ghl_eos_parameters *restrict eos,
      const harm_aux_vars_struct *restrict harm_aux,
      const int ndim,
      const double dummy,
      const double x[],
      double dx[],
      double resid[],
      double jac[][1],
      double *restrict f,
      double *restrict df,
      int *restrict n_iter);

/* Function    :  Hybrid_Noble1D()
 * Description :  Unpacks the ghl_primitive_quantities struct into the variables
 *                needed by the Newton-Rapson solver, then repacks the  primitives.
 *                This function is adapted from the HARM function provided by IllinoisGRMHD.
 * Documentation: 
 */

/*************************************************************************************

utoprim_1d.c:
---------------

    Uses the 1D_W method:
       -- solves for one independent variable (W) via a 1D
          Newton-Raphson method
       -- can be used (in principle) with a general equation of state.

******************************************************************************/

/**********************************************************************************

  ghl_hybrid_Noble1D():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the
        Newton-Raphson routine.

  -- assumes that
             /   rho gamma   \
         U = | alpha T^t_\mu |
             \   alpha B^i   /



             /     rho     \
      prim = |     uu      |
             | \tilde{u}^i |
             \  alpha B^i  /


return: i where
        i = 0 -> success
            1 -> initial v^2 < 0 with initial primitive guess;
            2 -> Newton-Raphson solver did not converge to a solution with the
                 given tolerances;
            3 -> Newton-Raphson procedure encountered a numerical divergence
                 (occurrence of "nan" or "+/-inf");
            4 -> Z<0 or Z>Z_TOO_BIG
            5 -> v^2 < 0 returned by the Newton-Raphson solver;
            6 -> computed pressure is negative

**********************************************************************************/

int ghl_hybrid_Noble1D(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics ) {

  double gnr_out[1];

  harm_aux_vars_struct harm_aux;

  // Calculate Z:
  if( ghl_initialize_Noble(eos, ADM_metric, metric_aux, cons_undens,
                           prims, &harm_aux, &gnr_out[0]) )
    return 1;

  // To be consistent with entropy variants, unused argument 0.0 is needed
  const int retval = ghl_general_newton_raphson(eos, &harm_aux, 1, 0.0, params->con2prim_max_iterations,
                                                params->con2prim_solver_tolerance, &diagnostics->n_iter,
                                                gnr_out, ghl_validate_1D, ghl_func_1D);

  const double Z = gnr_out[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if(retval != 0) {
    return retval;
  } else if(Z <= 0. || Z > Z_TOO_BIG) {
    return 4;
  }

  // Calculate v^2:
  double vsq = ghl_vsq_calc(&harm_aux, Z);
  if( vsq >= 1. ) {
    vsq = 1.-2.e-16;
  } else if(vsq < 0.0) {
    //v should be real!
    return 5;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  ghl_finalize_Noble(params, eos, ADM_metric, metric_aux, cons_undens, &harm_aux, Z, vsq, prims);
  if(!params->ignore_negative_pressure && prims->press <= 0.0) {
    return 6;
  }

  /* Done! */
  diagnostics->which_routine = Noble1D;
  return 0;
}
