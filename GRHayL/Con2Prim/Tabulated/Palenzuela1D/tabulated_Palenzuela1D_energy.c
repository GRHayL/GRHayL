#include "../utils.h"

/*
 * Function   : compute_rho_P_eps_T_W_energy
 * Author     : Leo Werneck
 *
 * This function is used to compute rho, P, eps, T, and W from x and the input
 * conservative variables, using the specific internal energy to recover the
 * temperature.
 *
 * Parameters : x       - Current guess for the root of f(x).
 *            : fparams - Parameters used to evaluate f(x) (see palenzuela.h).
 *            : rho_ptr - Stores the result for rho.
 *            : P_ptr   - Stores the result for P.
 *            : eps_ptr - Stores the result for eps.
 *            : T_ptr   - Stores the result for T.
 *            : W_ptr   - Stores the result for W.
 *
 * Returns    : Nothing.
 *
 * References : [1] Palenzuela et al. (2015) arXiv:1505.01607
 *            : [2] Siegel et al. (2018) arXiv:1712.07538
 */
static void
compute_rho_P_eps_T_W_energy(
      const double x,
      fparams_struct *restrict fparams,
      double *restrict rho_ptr,
      double *restrict P_ptr,
      double *restrict eps_ptr,
      double *restrict T_ptr,
      double *restrict W_ptr ) {

  // Step 0: Unpack the fparams struct
  const double Y_e = fparams->Y_e;
  double T = fparams->temp_guess;
  const eos_parameters *eos = fparams->eos;

  // Step 1: First compute rho and W
  double rho, W;
  compute_rho_W_from_x_and_conservatives( x, fparams, &rho, &W );

  // Step 2: Now compute P and eps
  double P, eps;
  if( fparams->evolve_T ) {
    const double q = fparams->q;
    const double s = fparams->s;
    const double t = fparams->t;
    // Eq. (43) of https://arxiv.org/pdf/1712.07538.pdf
    eps = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W) );
    ghl_tabulated_compute_P_T_from_eps( eos, rho, Y_e, eps, &P, &T );
  }
  else {
    // If the temperature is not evolved, use the input guess to determine
    // the remaining primitives. Note that in this case one must provide
    // the appropriate temperature instead of the default guess T = T_min.
    ghl_tabulated_compute_P_eps_from_T( eos, rho, Y_e, T, &P, &eps );
  }

  // Step 3: Set the output
  *rho_ptr = rho;
  *P_ptr   = P;
  *eps_ptr = eps;
  *T_ptr   = T;
  *W_ptr   = W;
}

/*
 * Function : ghl_tabulated_Palenzuela1D_energy
 * Author   : Leo Werneck
 *
 * This is a wrapper function around the ghl_tabulated_Palenzuela1D that performs
 * a primitive recovery using the Palenzuela et al. scheme using the specific
 * energy to recover the temperature. See file ghl_tabulated_Palenzuela1D.c for
 * further details.
 */
int ghl_tabulated_Palenzuela1D_energy(
      const grhayl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const ADM_aux_quantities *restrict metric_aux,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  return ghl_tabulated_Palenzuela1D(compute_rho_P_eps_T_W_energy,
                                params,
                                eos,
                                ADM_metric,
                                cons_undens,
                                prims,
                                diagnostics);
}
