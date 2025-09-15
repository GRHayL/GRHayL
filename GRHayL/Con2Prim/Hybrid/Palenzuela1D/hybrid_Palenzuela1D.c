#include "../../utils_Palenzuela1D.h"

/*
 * Function : froot_hybrid
 * Author   : Leo Werneck
 *
 * Computes Eq. (33) of Siegel et al., 2018 (arXiv: 1712.07538). The function
 * arguments follow the standards set by the roots.h file.
 *
 * Parameters : x        - The point at which f(x) is evaluated.
 *            : fparams  - Pointer to parameter structed containing auxiliary
 *                         variables needed by this function (see definition
 *                         above).
 *
 * Returns    : Nothing.
 */
static inline double
froot(
      const double x,
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_conservative_quantities *restrict cons_undens,
      fparams_struct *restrict fparams,
      ghl_primitive_quantities *restrict prims) {

  double W;
  fparams->compute_rho_P_eps_W(
    x, params, eos, cons_undens, fparams, prims, &W);

  // Eq: (33) of https://arxiv.org/pdf/1712.07538.pdf
  return x - (1.0 + prims->eps + prims->press/prims->rho)*W;
}

/* Function    : ghl_tabulated_Palenzuela1D
 * Author      : Leo Werneck
 *
 * Performs a primitive recovery using the Palenzuela et al. scheme.
 *
 * Parameters  : compute_rho_P_eps_W - Function pointer provided by the user
 *             : grhyal_params       - GRHayL parameters struct
 *             : eos                 - EOS parameters struct
 *             : metric              - ghl_metric_quantities struct containing local
 *                                     metric data (gamma_{ij}, g_{mu nu}, etc).
 *             : cons_undens         - ghl_conservative_quantities struct containing
 *                                     local undensitized conservative quantities.
 *             : prims               - ghl_primitive_quantities struct containing
 *                                     local hydrodynamic quantities. Stores the
 *                                     result if the root-finding succeeds.
 *             : diagnostics         - Pointer to con2prim diagnostics struct.
 *
 * References  : [1] Palenzuela et al. (2015) arXiv:1505.01607
 *             : [2] Siegel et al. (2018) arXiv:1712.07538
 */
ghl_error_codes_t ghl_hybrid_Palenzuela1D(
      void compute_rho_P_eps_W(
            const double x,
            const ghl_parameters *restrict params,
            const ghl_eos_parameters *restrict eos,
            const ghl_conservative_quantities *restrict cons_undens,
            fparams_struct *restrict fparams,
            ghl_primitive_quantities *restrict prims,
            double *restrict W_ptr),
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  double SU[3], B_squared, S_squared, BdotS;
  compute_SU_Bsq_Ssq_BdotS(ADM_metric, cons_undens, prims, SU, &B_squared, &S_squared, &BdotS);
  const double tau = MAX(cons_undens->tau, eos->tau_atm);

  // Set specific quantities for this routine (Eq. A7 of [1])
  fparams_struct fparams;
  const double invD           = 1.0/cons_undens->rho;
  fparams.compute_rho_P_eps_W = compute_rho_P_eps_W;
  fparams.q                   = tau * invD;
  fparams.r                   = S_squared * invD * invD;
  fparams.s                   = B_squared * invD;
  fparams.t                   = BdotS/pow(cons_undens->rho, 1.5);

  // Bracket x (Eq. A8 of [1])
  double xlow = 1 + fparams.q - fparams.s;
  double xup  = 2*(1 + fparams.q) - fparams.s;

  // Call the main function and perform the con2prim
  roots_params rparams;
  rparams.tol = 1e-15;
  rparams.max_iters = 300;
  ghl_error_codes_t error = ghl_brent(froot, params, eos, cons_undens, &fparams, prims, xlow, xup, &rparams);
  if(error)
    return error;

  diagnostics->n_iter = rparams.n_iters;

  // Set core primitives using the EOS and the root
  double x = rparams.root;
  double W;
  compute_rho_P_eps_W(x, params, eos, cons_undens, &fparams, prims, &W);

  // Set Z
  const double Z = x*prims->rho*W;

  // Compute utilde^{i}
  double utildeU[3] = {
    W*(SU[0] + BdotS*prims->BU[0]/Z)/(Z+B_squared),
    W*(SU[1] + BdotS*prims->BU[1]/Z)/(Z+B_squared),
    W*(SU[2] + BdotS*prims->BU[2]/Z)/(Z+B_squared)
  };

  // Set prims struct
  diagnostics->speed_limited = ghl_limit_utilde_and_compute_v(params, ADM_metric, utildeU, prims);
  prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);

  return ghl_success;
}
