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

  ghl_tabulated_primitive_guess_aux aux = { 0 };
  ghl_tabulated_compute_primitive_guess_auxiliaries(metric_adm, cons_undens, prims, &aux);

  // Compute the lower bound of the x variable, Eq. (35) of 1712.07538
  const double x = 1.0 + aux.q - aux.s;

  // Complete the primitive guess using Eqs. (24), (42), (43), and (44) of
  // 1712.07538. We set W = 0.0 as a guess to ensure it's computed from x
  // and the conserved quantities in the function below.
  const double W = 0.0;
  ghl_tabulated_primitive_guess_from_x_and_W(params, eos, metric_adm,
                                             cons_undens, &aux, x, W, prims);
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
