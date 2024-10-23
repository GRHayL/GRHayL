#include "../../utils_Palenzuela1D.h"

/*
 * Function   : compute_rho_P_eps_W_entropy
 * Author     : Leo Werneck
 *
 * This function is used to compute rho, P, eps, and W from x and the input
 * conservative variables, using the specific internal energy to recover the
 * temperature.
 *
 * Parameters : x       - Current guess for the root of f(x).
 *            : fparams - Parameters used to evaluate f(x) (see palenzuela.h).
 *            : rho_ptr - Stores the result for rho.
 *            : P_ptr   - Stores the result for P.
 *            : eps_ptr - Stores the result for eps.
 *            : W_ptr   - Stores the result for W.
 *
 * Returns    : Nothing.
 *
 * References : [1] Palenzuela et al. (2015) arXiv:1505.01607
 *            : [2] Siegel et al. (2018) arXiv:1712.07538
 */
static void
compute_rho_P_eps_W_entropy(
      const double x,
      palenzuela_params *restrict params,
      double *restrict rho_ptr,
      double *restrict P_ptr,
      double *restrict eps_ptr,
      double *restrict W_ptr) {

  // Step 1: First compute rho and W
  double rho, W;
  compute_rho_W_from_x_and_conservatives(x, params, &rho, &W);

  // Step 2: Compute P (Eq. 19 of arXiv:0808.3140)
  const double ent   = params->cons_undens->entropy / W;
  const double Gamma = params->eos->Gamma_ppoly[ghl_hybrid_find_polytropic_index(params->eos, rho)];
  double P           = ent * pow(rho, Gamma - 1.0);

  // Step 3: Compute the pressure
  double P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(params->eos, rho, &P_cold, &eps_cold);
  // Don't let P be less than P_cold
  if(P < P_cold) P = P_cold;
  // Compute eps
  const double eps = eps_cold + (P - P_cold) / ( (params->eos->Gamma_th - 1.0) * rho);

  // Step 4: Set the output
  *rho_ptr = rho;
  *P_ptr   = P;
  *eps_ptr = eps;
  *W_ptr   = W;
}

/*
 * Function : ghl_hybrid_Palenzuela1D_energy
 * Author   : Leo Werneck
 *
 * This is a wrapper function around the ghl_hybrid_Palenzuela1D that performs
 * a primitive recovery using the Palenzuela et al. scheme using the specific
 * energy. See file ghl_hybrid_Palenzuela1D.c for further details.
 */
ghl_error_codes_t ghl_hybrid_Palenzuela1D_entropy(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  diagnostics->which_routine = Palenzuela1D;
  return ghl_hybrid_Palenzuela1D(compute_rho_P_eps_W_entropy,
                                 params,
                                 eos,
                                 ADM_metric,
                                 cons_undens,
                                 prims,
                                 diagnostics);
}
