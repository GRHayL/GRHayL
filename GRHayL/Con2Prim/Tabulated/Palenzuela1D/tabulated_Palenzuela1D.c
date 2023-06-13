#include "../utils.h"

/* Function    : ghl_tabulated_Palenzuela1D
 * Author      : Leo Werneck
 *
 * Performs a primitive recovery using the Palenzuela et al. scheme.
 *
 * Parameters  : compute_rho_P_eps_T_W - Function pointer provided by the user
 *             : grhyal_params         - GRHayL parameters struct
 *             : eos                   - EOS parameters struct
 *             : metric                - metric_quantities struct containing local
 *                                       metric data (gamma_{ij}, g_{mu nu}, etc).
 *             : cons_undens           - conservative_quantities struct containing
 *                                       local undensitized conservative quantities.
 *             : prims                 - primitive_quantities struct containing
 *                                       local hydrodynamic quantities. Stores the
 *                                       result if the root-finding succeeds.
 *             : diagnostics           - Pointer to con2prim diagnostics struct.
 *
 * Returns     : error_key (see GRHayL.h and roots.h).
 *
 * References  : [1] Palenzuela et al. (2015) arXiv:1505.01607
 *             : [2] Siegel et al. (2018) arXiv:1712.07538
 */
int ghl_tabulated_Palenzuela1D(
      void compute_rho_P_eps_T_W(
            const double x,
            fparams_struct *restrict fparams,
            double *restrict rho_ptr,
            double *restrict P_ptr,
            double *restrict eps_ptr,
            double *restrict T_ptr,
            double *restrict W_ptr ),
      const ghl_parameters *restrict params,
      const eos_parameters *restrict eos,
      const metric_quantities *restrict ADM_metric,
      const conservative_quantities *restrict cons_undens,
      primitive_quantities *restrict prims,
      con2prim_diagnostics *restrict diagnostics ) {

  // Step 1: Compute S^{2} = gamma^{ij}S_{i}S_{j}
  double SD[3] = {cons_undens->SD[0], cons_undens->SD[1], cons_undens->SD[2]};
  double S_squared = ghl_compute_vec2_from_vecD(ADM_metric->gammaUU, SD);

  // Step 2: Enforce ceiling on S^{2} (Eq. A5 of [1])
  // Step 2.1: Compute maximum allowed value for S^{2}
  const double S_squared_max = SQR(cons_undens->tau + cons_undens->rho);
  if( S_squared > S_squared_max ) {
    // Step 2.2: Rescale S_{i}
    const double rescale_factor = sqrt(0.9999*S_squared_max/S_squared);
    for(int i=0;i<3;i++)
      SD[i] *= rescale_factor;

    // Step 2.3: Recompute S^{2}
    S_squared = ghl_compute_vec2_from_vecD(ADM_metric->gammaUU, SD);
  }

  // Step 3: Compute B^{2} = gamma_{ij}B^{i}B^{j}
  const double BbarU[3] = {prims->BU[0] * ONE_OVER_SQRT_4PI,
                           prims->BU[1] * ONE_OVER_SQRT_4PI,
                           prims->BU[2] * ONE_OVER_SQRT_4PI};
  const double B_squared = ghl_compute_vec2_from_vecU(ADM_metric->gammaDD, BbarU);

  // Step 4: Compute B.S = B^{i}S_{i}
  double BdotS = 0.0;
  for(int i=0;i<3;i++) BdotS += BbarU[i]*SD[i];

  // Step 5: Set specific quantities for this routine (Eq. A7 of [1])
  fparams_struct fparams;
  const double invD             = 1.0/cons_undens->rho;
  fparams.compute_rho_P_eps_T_W = compute_rho_P_eps_T_W;
  fparams.evolve_T              = params->evolve_temp;
  fparams.Y_e                   = cons_undens->Y_e * invD;
  fparams.q                     = cons_undens->tau * invD;
  fparams.r                     = S_squared * invD * invD;
  fparams.s                     = B_squared * invD;
  fparams.t                     = BdotS/pow(cons_undens->rho, 1.5);
  fparams.eos                   = eos;
  fparams.cons_undens           = cons_undens;

  // Step 6: Bracket x (Eq. A8 of [1])
  double xlow = 1 + fparams.q - fparams.s;
  double xup  = 2*(1 + fparams.q) - fparams.s;

  // Step 7: Set initial guess for temperature
  fparams.temp_guess = prims->temperature;

  // Step 8: Call the main function and perform the con2prim
  roots_params rparams;
  rparams.tol = 1e-15;
  rparams.max_iters = 300;
  ghl_toms748(froot, &fparams, xlow, xup, &rparams);
  if( rparams.error_key != roots_success ) {
    // Adjust the temperature guess and try again
    fparams.temp_guess = eos->T_max;
    prims->temperature = eos->T_max;
    ghl_toms748(froot, &fparams, xlow, xup, &rparams);
    if( rparams.error_key != roots_success ) {
      ghl_warn("TOMS748 did not find a root (see error report below)\n");
      roots_info(&rparams);
      return rparams.error_key;
    }
  }

  // Step 9: Set core primitives using the EOS and the root
  double x = rparams.root;
  double rho, P, eps, T, W;
  compute_rho_P_eps_T_W(x, &fparams, &rho, &P, &eps, &T, &W);

  // Step 10: Compute auxiliary quantities to obtain the velocities using
  //          Eq. (24) in [2]. Note, however, that GRHayL expects the velocity
  //          tilde(u)^{i} := W v^{i},
  // Step 10.a: Compute S^{i}
  double SU[3];
  ghl_raise_vector_3D(ADM_metric->gammaUU, SD, SU);

  // Step 10.b: Set Z
  const double Z = x*rho*W;

  // Step 10.c: Set prims struct
  prims->rho         = rho;
  prims->Y_e         = fparams.Y_e;
  prims->temperature = T;
  prims->press       = P;
  prims->eps         = eps;
  prims->vU[0]          = W*(SU[0] + BdotS*BbarU[0]/Z)/(Z+B_squared);
  prims->vU[1]          = W*(SU[1] + BdotS*BbarU[1]/Z)/(Z+B_squared);
  prims->vU[2]          = W*(SU[2] + BdotS*BbarU[2]/Z)/(Z+B_squared);

  return ghl_success;
}
