#include "../../utils_Palenzuela1D.h"

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
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_conservative_quantities *restrict cons_undens,
      fparams_struct *restrict fparams,
      ghl_primitive_quantities *restrict prims,
      double *restrict W_ptr) {

  // Step 1: First compute rho and W
  double W;
  compute_rho_W_from_x_and_conservatives(x, params, cons_undens, fparams, &prims->rho, &W);

  // Step 2: Now compute P and eps
  if(fparams->evolve_T) {
    const double q = fparams->q;
    const double s = fparams->s;
    const double t = fparams->t;
    // Eq. (43) of https://arxiv.org/pdf/1712.07538.pdf
    prims->eps = W - 1.0 + (1.0-W*W)*x/W + W*(q - s + t*t/(2*x*x) + s/(2*W*W) );

    // Enforce limits and call table interpolator
    ghl_tabulated_enforce_bounds_rho_Ye_eps(eos, &prims->rho, &prims->Y_e, &prims->eps);
    ghl_tabulated_compute_P_T_from_eps(eos, prims->rho, prims->Y_e, prims->eps, &prims->press, &prims->temperature);
  }
  else {
    // If the temperature is not evolved, use the input guess to determine
    // the remaining primitives. Note that in this case one must provide
    // the appropriate temperature instead of the default guess T = T_min.
    ghl_tabulated_enforce_bounds_rho_Ye_T(eos, &prims->rho, &prims->Y_e, &prims->temperature);
    ghl_tabulated_compute_P_eps_from_T(eos, prims->rho, prims->Y_e, prims->temperature, &prims->press, &prims->eps);
  }

  // Step 3: Set the output
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
ghl_error_codes_t ghl_tabulated_Palenzuela1D_energy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  diagnostics->which_routine = Palenzuela1D;
  return ghl_tabulated_Palenzuela1D(compute_rho_P_eps_T_W_energy,
                                params,
                                eos,
                                ADM_metric,
                                cons_undens,
                                prims,
                                diagnostics);
}
