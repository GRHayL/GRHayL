#include "../../utils_Palenzuela1D.h"

/**
 * @ingroup c2p_internal
 * @brief Evaluate the Palenzuela residual x - hW for the selected EOS callback.
 * @author Leo Werneck
 *
 * Implements Eq. (33) of Siegel et al. (2018), arXiv:1712.07538.
 * @param[in] x Trial value of hW
 * @param[in] params Recovery and velocity-limit parameters
 * @param[in] eos EOS parameters
 * @param[in] cons_undens Undensitized conservatives
 * @param[in,out] fparams Normalized contractions and thermodynamic callback
 * @param[in,out] prims Primitive state evaluated at x
 * @return Residual x - (1 + eps + P/rho) W
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
  fparams->compute_rho_P_eps_T_W(
    x, params, eos, cons_undens, fparams, prims, &W);

  // Eq: (33) of https://arxiv.org/pdf/1712.07538.pdf
  return x - (1.0 + prims->eps + prims->press/prims->rho)*W;
}

/**
 * @ingroup c2p_internal
 * @brief Recover tabulated primitives with a Palenzuela thermodynamic callback.
 * @author Leo Werneck
 *
 * @param[in] compute_rho_P_eps_T_W Internal energy/entropy callback for trial x
 * @param[in] params Recovery and velocity-limit parameters
 * @param[in] eos EOS parameters
 * @param[in] metric_adm ADM metric
 * @param[in] cons_undens Undensitized conservatives
 * @param[in,out] prims Supplies magnetic fields and receives recovered primitives
 * @param[in,out] diagnostics Initialized recovery diagnostics
 * @return ghl_success on recovery; otherwise a root-finding or EOS error.
 *
 * See Palenzuela et al. (2015), arXiv:1505.01607, and Siegel et al. (2018),
 * arXiv:1712.07538. This private callback interface is not an installed API.
 */
ghl_error_codes_t ghl_tabulated_Palenzuela1D(
      void compute_rho_P_eps_T_W(
            const double x,
            const ghl_parameters *restrict params,
            const ghl_eos_parameters *restrict eos,
            const ghl_conservative_quantities *restrict cons_undens,
            fparams_struct *restrict fparams,
            ghl_primitive_quantities *restrict prims,
            double *restrict W_ptr),
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims,
      ghl_con2prim_diagnostics *restrict diagnostics) {

  diagnostics->speed_limited = false;

  double SU[3], B_squared, S_squared, BdotS;
  ghl_compute_SU_Bsq_Ssq_BdotS(metric_adm, cons_undens, prims, SU, &B_squared, &S_squared, &BdotS);
  const double tau = fmax(cons_undens->tau, eos->tau_atm);

  // Set specific quantities for this routine (Eq. A7 of [1])
  fparams_struct fparams;
  const double invD             = 1.0/cons_undens->rho;
  fparams.compute_rho_P_eps_T_W = compute_rho_P_eps_T_W;
  fparams.evolve_T              = params->evolve_temp;
  fparams.q                     = tau * invD;
  fparams.r                     = S_squared * invD * invD;
  fparams.s                     = B_squared * invD;
  fparams.t                     = BdotS/pow(cons_undens->rho, 1.5);

  // Bracket x (Eq. A8 of [1])
  double xlow = 1 + fparams.q - fparams.s;
  double xup  = 2*(1 + fparams.q) - fparams.s;

  // Call the main function and perform the con2prim
  roots_params rparams;
  rparams.tol = 1e-15;
  rparams.max_iters = 300;
  ghl_error_codes_t error = ghl_brent(froot, params, eos, cons_undens, &fparams, prims, xlow, xup, &rparams);
  if(error) {
    // Adjust the temperature guess and try again
    prims->temperature = eos->T_min;
    error = ghl_brent(froot, params, eos, cons_undens, &fparams, prims, xlow, xup, &rparams);
    if(error)
      return error;
  }
  diagnostics->n_iter = rparams.n_iters;

  // Set core primitives using the EOS and the root
  double x = rparams.root;
  double W;
  compute_rho_P_eps_T_W(x, params, eos, cons_undens, &fparams, prims, &W);
  // Brent may probe points where the table inversion fails; only the
  // evaluation at the root sets the returned primitives
  if(fparams.eos_error) {
    return fparams.eos_error;
  }

  // Set Z
  const double Z = x*prims->rho*W;

  // Compute utilde^{i}
  double utildeU[3] = {
    W*(SU[0] + BdotS*prims->BU[0]/Z)/(Z+B_squared),
    W*(SU[1] + BdotS*prims->BU[1]/Z)/(Z+B_squared),
    W*(SU[2] + BdotS*prims->BU[2]/Z)/(Z+B_squared)
  };

  // Set prims struct
  ghl_tabulated_enforce_bounds_rho_Ye_T(eos, &prims->rho, &prims->Y_e, &prims->temperature);
  diagnostics->speed_limited = ghl_limit_utilde_and_compute_v(params, metric_adm, utildeU, prims);
  return ghl_tabulated_compute_P_eps_S_from_T(eos, prims->rho, prims->Y_e, prims->temperature,
                                              &prims->press, &prims->eps, &prims->entropy);
}
