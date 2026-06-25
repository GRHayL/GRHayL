#include "ghl_con2prim.h"

/**
 * @ingroup Con2Prim
 * @brief Computes an initial guess for the @ref Con2Prim solvers
 *
 * @details
 * This function sets a default initial guess for hybrid or simple EOSs. It can
 * be used whenever a better guess is not available. It sets
 *
 * \f[
 * \begin{aligned}
 * \rho &= \frac{\rho_*}{\sqrt{|\gamma|}} \\
 * u^0 &= 1 \\
 * v^i &= -\beta^i \\
 * Y_e &= \frac{\tilde{Y_e}}{\rho_*} \\
 * T &= T_\mathrm{max}
 * \end{aligned}
 * \f]
 *
 * This choice sets the transport/utilde velocity \f$ v^i+\beta^i \f$ to zero;
 * it is only an initial guess for the Con2Prim solve. We set the pressure and
 * specific internal energy \f$ \epsilon \f$ to the cold values.
 *
 * @param[in] params pointer to ghl_parameters struct
 *
 * @param[in] eos pointer to ghl_eos_parameters struct
 *
 * @param[in] metric_adm pointer to ghl_metric_quantities struct with ADM metric
 *
 * @param[in] cons_undens pointer to ghl_conservative_quantities struct with
 *                        **undensitized** conservative variables
 *
 * @param[out] prims pointer to ghl_primitive_quantities containing the initial guess
 */
static void ghl_guess_primitives_hybrid_simple(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims) {

  // params is not used by this function
  (void)params;

  // Use atmosphere as initial guess:
  prims->rho   = cons_undens->rho;
  prims->u0    = 1.0;
  prims->vU[0] = -metric_adm->betaU[0];
  prims->vU[1] = -metric_adm->betaU[1];
  prims->vU[2] = -metric_adm->betaU[2];

  // Compute remaining primitives
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &prims->press, &prims->eps);
}

/**
 * @ingroup Con2Prim
 * @brief Computes an initial guess for the @ref Con2Prim solvers
 *
 * @details
 * This function sets a default initial guess for tabulated EOSs. It can be
 * used whenever a better guess is not available. It sets primitives following
 * the primitive-recovery strategy of the Palenzuela et al. routine, as outlined
 * in Siegel et al. (2018; https://arxiv.org/pdf/1712.07538).
 *
 * The only required guess is the temperature. We use \f$ T = T_\mathrm{max} \f$,
 * which uses ghl_eos_parameters::T_max.
 *
 * @param[in] params pointer to ghl_parameters struct
 *
 * @param[in] eos pointer to ghl_eos_parameters struct
 *
 * @param[in] metric_adm pointer to ghl_metric_quantities struct with ADM metric
 *
 * @param[in] cons_undens pointer to ghl_conservative_quantities struct with
 *                        **undensitized** conservative variables
 *
 * @param[out] prims pointer to ghl_primitive_quantities containing the initial guess
 */
static void ghl_guess_primitives_tabulated(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims) {

  // Compute auxiliary variables
  double SU[3], B_squared, S_squared, BdotS;
  ghl_compute_SU_Bsq_Ssq_BdotS(metric_adm, cons_undens, prims,
                               SU, &B_squared, &S_squared, &BdotS);

  // Compute Palenzuela variables {q, r, s, t}
  const double invD = 1.0 / cons_undens->rho;
  const double q = cons_undens->tau * invD;
  const double r = S_squared * invD * invD;
  const double s = B_squared * invD;
  const double t = BdotS / pow(cons_undens->rho, 1.5);

  // Compute the lower bound of the x variable, Eq. (35) of 1712.07538
  const double x = 1.0 + q - s;

  // Compute W, Eq. (42) of 1712.07538
  const double tmp_numer = x * x * r + (2.0 * x + s) * t * t;
  const double tmp_denom = x * x * (x + s) * (x + s);
  double Wminus2 = 1.0 - tmp_numer / tmp_denom;
  Wminus2 = ghl_clamp(Wminus2, params->inv_sq_max_Lorentz_factor, 1.0);
  const double W = pow(Wminus2, -0.5);

  // Compute rho Y_e, and u^0
  prims->rho = cons_undens->rho / W;
  prims->Y_e = cons_undens->Y_e / cons_undens->rho;
  prims->u0  = W * metric_adm->lapseinv;

  // Compute eps, Eq. (43-44) of 1712.07538
  prims->eps = - 1.0 + (1.0 - W * W) * x / W
             + W * (1.0 + q - s + t * t / (2.0 * x * x) + s / (2.0 * W * W));

  // Enforce table bounds, then compute missing hydrodynamic quantities, using
  // T = T_max as an initial guess.
  prims->temperature = eos->T_max;
  ghl_tabulated_enforce_bounds_rho_Ye_eps(eos, &prims->rho, &prims->Y_e, &prims->eps);
  ghl_tabulated_compute_P_S_T_from_eps(eos, prims->rho, prims->Y_e, prims->eps,
                                       &prims->press, &prims->entropy, &prims->temperature);

  // Compute the velocities. We start by computing the Valencia velocity v_n^i
  // using Eq. (24) of 1712.07538, with z = x rho W. We use it to compute
  // \tilde{u}^i = W * v_n^i, following the discussion below Eq. (15) of
  // astro-ph/0512420. Finally, we use a GRHayL function to limit the velocity
  // and compute the appropriate primitive variable v^i = u^i / u^0.
  const double z = x * prims->rho * W;
  double utildeU[3] = {
    W * (SU[0] + BdotS * prims->BU[0] / z) / (z + B_squared),
    W * (SU[1] + BdotS * prims->BU[1] / z) / (z + B_squared),
    W * (SU[2] + BdotS * prims->BU[2] / z) / (z + B_squared),
  };
  ghl_limit_utilde_and_compute_v(params, metric_adm, utildeU, prims);
}

/**
 * @ingroup Con2Prim
 * @brief Computes an initial guess for the @ref Con2Prim solvers
 *
 * @param[in] params pointer to ghl_parameters struct
 *
 * @param[in] eos pointer to ghl_eos_parameters struct
 *
 * @param[in] metric_adm pointer to ghl_metric_quantities struct with ADM metric
 *
 * @param[in] cons_undens pointer to ghl_conservative_quantities struct with
 *                        **undensitized** conservative variables
 *
 * @param[out] prims pointer to ghl_primitive_quantities containing the initial guess
 */
void ghl_guess_primitives(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict metric_adm,
      const ghl_conservative_quantities *restrict cons_undens,
      ghl_primitive_quantities *restrict prims) {

  switch(eos->eos_type) {
    case ghl_eos_hybrid:
    case ghl_eos_simple:
      ghl_guess_primitives_hybrid_simple(params, eos, metric_adm, cons_undens, prims);
      break;
    case ghl_eos_tabulated:
      ghl_guess_primitives_tabulated(params, eos, metric_adm, cons_undens, prims);
      break;
  }
}
