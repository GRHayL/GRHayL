#include "../../utils_Palenzuela1D.h"
#include "ghl_roots.h"

/*
 * Function : froot
 * Author   : Leo Werneck
 *
 * Computes Eq. (33) of Siegel et al., 2018 (arXiv: 1712.07538). The function
 * arguments follow the standards set by the ghl_roots.h file.
 *
 * Parameters : x        - The point at which f(x) is evaluated.
 *            : fparams  - Pointer to parameter structed containing auxiliary
 *                         variables needed by this function (see definition
 *                         above).
 *
 * Returns    : Nothing.
 */
static inline double froot(const double x, void *restrict params_in) {

  palenzuela_params *params = (palenzuela_params *)params_in;

  double rho, P, eps, W;
  params->compute_rho_P_eps_W(x, params, &rho, &P, &eps, &W);

  // Eq: (33) of https://arxiv.org/pdf/1712.07538.pdf
  return x - (1.0 + eps + P / rho) * W;
}

/* Function    : ghl_hybrid_Palenzuela1D
 * Author      : Leo Werneck
 *
 * Performs a primitive recovery using the Palenzuela et al. scheme.
 *
 * Parameters  : compute_rho_P_eps_W - Function pointer provided by the user
 *             : grhyal_params       - GRHayL parameters struct
 *             : eos                 - EOS parameters struct
 *             : metric              - ghl_metric_quantities struct containing
 * local metric data (gamma_{ij}, g_{mu nu}, etc). : cons_undens         -
 * ghl_conservative_quantities struct containing local undensitized conservative
 * quantities. : prims               - ghl_primitive_quantities struct
 * containing local hydrodynamic quantities. Stores the result if the
 * root-finding succeeds. : diagnostics         - Pointer to con2prim
 * diagnostics struct.
 *
 * References  : [1] Palenzuela et al. (2015) arXiv:1505.01607
 *             : [2] Siegel et al. (2018) arXiv:1712.07538
 */
ghl_error_codes_t ghl_hybrid_Palenzuela1D(
      void compute_rho_P_eps_W(
            const double x,
            palenzuela_params *restrict params,
            double *restrict rho_ptr,
            double *restrict P_ptr,
            double *restrict eps_ptr,
            double *restrict W_ptr),
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  // Step 1: Compute S^{2} = gamma^{ij}S_{i}S_{j}
  double SD[3] = { cons_undens->SD[0], cons_undens->SD[1], cons_undens->SD[2] };
  double S_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, SD);

  // Step 2: Enforce ceiling on S^{2} (Eq. A5 of [1])
  // Step 2.1: Compute maximum allowed value for S^{2}
  const double S_squared_max = SQR(cons_undens->tau + cons_undens->rho);
  if(S_squared > S_squared_max) {
    // Step 2.2: Rescale S_{i}
    const double rescale_factor = sqrt(0.9999 * S_squared_max / S_squared);
    for(int i = 0; i < 3; i++) {
      SD[i] *= rescale_factor;
    }

    // Step 2.3: Recompute S^{2}
    S_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaUU, SD);
  }

  // Step 3: Compute B^{2} = gamma_{ij}B^{i}B^{j}
  const double B_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, prims->BU);

  // Step 4: Compute B.S = B^{i}S_{i}
  double BdotS = 0.0;
  for(int i = 0; i < 3; i++) {
    BdotS += prims->BU[i] * SD[i];
  }

  // Step 5: Set specific quantities for this routine (Eq. A7 of [1])
  palenzuela_params pparams;
  const double invD = 1.0 / cons_undens->rho;
  pparams.compute_rho_P_eps_W = compute_rho_P_eps_W;
  pparams.q = cons_undens->tau * invD;
  pparams.r = S_squared * invD * invD;
  pparams.s = B_squared * invD;
  pparams.t = BdotS / pow(cons_undens->rho, 1.5);
  pparams.cons_undens = cons_undens;
  pparams.ghl_params = params;
  pparams.eos = eos;

  // Step 6: Bracket x (Eq. A8 of [1])
  double xlow = 1 + pparams.q - pparams.s;
  double xup = 2 * (1 + pparams.q) - pparams.s;

  // Step 7: Call the main function and perform the con2prim
  roots_params rparams;
  rparams.tol = 1e-15;
  rparams.max_iters = 300;
  // ghl_toms748(froot, &fparams, xlow, xup, &rparams);
  ghl_error_codes_t error = ghl_brent(froot, &pparams, xlow, xup, &rparams);
  if(error)
    return error;

  diagnostics->n_iter = rparams.n_iters;

  // Step 8: Set core primitives using the EOS and the root
  double x = rparams.root;
  double rho, P, eps, W;
  compute_rho_P_eps_W(x, &pparams, &rho, &P, &eps, &W);

  // Step 9: Compute auxiliary quantities to obtain the velocities using
  //          Eq. (24) in [2]. Note, however, that GRHayL expects the velocity
  //          tilde(u)^{i} := W v^{i},
  // Step 9.a: Compute S^{i}
  double SU[3];
  ghl_raise_lower_vector_3D(ADM_metric->gammaUU, SD, SU);

  // Step 9.b: Set Z
  const double Z = x*prims->rho*W;

  // Step 9.c: Compute utilde^{i}
  double utildeU[3] = { W * (SU[0] + BdotS * prims->BU[0] / Z) / (Z + B_squared),
                        W * (SU[1] + BdotS * prims->BU[1] / Z) / (Z + B_squared),
                        W * (SU[2] + BdotS * prims->BU[2] / Z) / (Z + B_squared) };

  // Step 9.d: Set prims struct
  diagnostics->speed_limited = ghl_limit_utilde_and_compute_v(params, ADM_metric, utildeU, prims);
  prims->entropy = ghl_hybrid_compute_entropy_function(eos, prims->rho, prims->press);

  return ghl_success;
}
