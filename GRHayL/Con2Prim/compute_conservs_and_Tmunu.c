#include "ghl_con2prim.h"

/**
 * @ingroup Con2Prim
 * Computes conservatives and \f$ T_{\mu\nu} \f$ from primitives
 *
 * @details
 * This function computes the conservative variables and stress-energy tensor
 * from the primitive variables, which is recommended after applying primitive
 * limits to ensure that the two variable sets are self-consistent.
 *
 * @param[in] ADM_metric: pointer to ghl_metric_quantities struct with ADM metric
 *
 * @param[in] metric_aux: pointer to ghl_ADM_aux_quantities struct
 *
 * @param[in] prims: pointer to ghl_primitive_quantities struct; note that
 *                   this function requires that ghl_primitive_quantities::eps
 *                   and ghl_primitive_quantities::u0 must be initialized for
 *                   this function.
 *
 * @param[out] cons: pointer to ghl_conservative_quantities struct
 *
 * @param[out] Tmunu: pointer to ghl_stress_energy struct
 *
 * @returns void
 */
void ghl_compute_conservs_and_Tmunu(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons,
      ghl_stress_energy *restrict Tmunu) {

  /**
   * ## Step-by-step Procedure
   *
   * To compute the conservatives, we first compute the enthalpy
   *
   * \f[
   * h = 1 + \epsilon + \frac{P}{\rho}
   * \f]
   */
  const double h_enthalpy = 1.0 + prims->eps + prims->press/prims->rho;

  // Compute u^i. u^0 is provided to the function.
  /**
   * and the 4-vector and its dual:
   *
   * \f[
   * \begin{aligned}
   * u^i = \frac{v^i}{u^0} \\
   * u_i = \gamma_{ij} u^j
   * \end{aligned}
   * \f]
   *
   * where the dual is computed using @ref ghl_raise_lower_vector_4D.
   */
  const double uU[4] = {prims->u0,
                        prims->u0*prims->vU[0],
                        prims->u0*prims->vU[1],
                        prims->u0*prims->vU[2]};

  double uD[4];
  ghl_raise_lower_vector_4D(metric_aux->g4DD, uU, uD);

  /**
   * We then compute \f$ b^\mu \f$, \f$ b_\mu \f$, and \f$ b^2 \f$ using
   * @ref ghl_compute_smallb_and_b2 and @ref ghl_raise_lower_vector_4D.
   * Next, we compute some intermediate quantities and calculate
   * the conservative variables:
   */
  double smallb[4], smallb2;
  ghl_compute_smallb_and_b2(ADM_metric, prims, uD, smallb, &smallb2);

  // Precompute some useful quantities, for later:
  const double alpha_sqrt_gamma = ADM_metric->lapse*ADM_metric->sqrt_detgamma;
  const double rho0_h_plus_b2 = (prims->rho*h_enthalpy + smallb2);
  const double P_plus_half_b2 = (prims->press+0.5*smallb2);

  double smallb_lower[4];
  ghl_raise_lower_vector_4D(metric_aux->g4DD, smallb, smallb_lower);

  /**
   * \f[
   * \begin{aligned}
   * \rho_*       &= \alpha \sqrt{| \gamma |} \rho u^0 \\
   * \tilde{S}_i  &= h \rho u_i + \alpha \sqrt{| \gamma |} \left( b^2 u^0 u_i - b^0 b_i \right) \\
   * \tilde{\tau} &= \alpha \sqrt{| \gamma |} \left[ \left( h \rho + b^2 \right) \left( u^0 \right)^2
   *                                              - \frac{P}{\alpha^2} - \frac{b^2}{2\alpha^2}
   *                                              - \left( b^0 \right)^2 \right]
   *               - \rho_* \\
   * \tilde{S}    &= \alpha \sqrt{| \gamma |} S u^0 \\
   * \tilde{Y_e}  &= \rho_* Y_e
   * \end{aligned}
   * \f]
   */
  cons->rho = alpha_sqrt_gamma * prims->rho * uU[0];
  cons->SD[0] = cons->rho*h_enthalpy*uD[1] + alpha_sqrt_gamma*(uU[0]*smallb2*uD[1] - smallb[0]*smallb_lower[1]);
  cons->SD[1] = cons->rho*h_enthalpy*uD[2] + alpha_sqrt_gamma*(uU[0]*smallb2*uD[2] - smallb[0]*smallb_lower[2]);
  cons->SD[2] = cons->rho*h_enthalpy*uD[3] + alpha_sqrt_gamma*(uU[0]*smallb2*uD[3] - smallb[0]*smallb_lower[3]);
  cons->tau = ADM_metric->lapse*alpha_sqrt_gamma*(rho0_h_plus_b2*SQR(uU[0]) - P_plus_half_b2*ADM_metric->lapseinv2 - SQR(smallb[0])) - cons->rho;
  cons->entropy = alpha_sqrt_gamma * prims->entropy * uU[0];
  cons->Y_e = cons->rho * prims->Y_e;

  /**
   * Finally, we compute the stress-energy tensor:
   *
   * \f[
   * T_{\mu\nu} = (h \rho + b^2) u_\mu u_\nu
   *            + (P + \frac{b^2}{2}) g_{\mu\nu}
   *            - b_\mu b_\nu
   * \f]
   *
   * We don't use the @grhayl-provided function for computing \f$ T_{\mu\nu} \f$
   * because we can reuse the precomputed quantities.
   */
  Tmunu->T4[0][0] = rho0_h_plus_b2*uD[0]*uD[0] + P_plus_half_b2*metric_aux->g4DD[0][0] - smallb_lower[0]*smallb_lower[0];
  Tmunu->T4[0][1] = rho0_h_plus_b2*uD[0]*uD[1] + P_plus_half_b2*metric_aux->g4DD[0][1] - smallb_lower[0]*smallb_lower[1];
  Tmunu->T4[0][2] = rho0_h_plus_b2*uD[0]*uD[2] + P_plus_half_b2*metric_aux->g4DD[0][2] - smallb_lower[0]*smallb_lower[2];
  Tmunu->T4[0][3] = rho0_h_plus_b2*uD[0]*uD[3] + P_plus_half_b2*metric_aux->g4DD[0][3] - smallb_lower[0]*smallb_lower[3];
  Tmunu->T4[1][1] = rho0_h_plus_b2*uD[1]*uD[1] + P_plus_half_b2*metric_aux->g4DD[1][1] - smallb_lower[1]*smallb_lower[1];
  Tmunu->T4[1][2] = rho0_h_plus_b2*uD[1]*uD[2] + P_plus_half_b2*metric_aux->g4DD[1][2] - smallb_lower[1]*smallb_lower[2];
  Tmunu->T4[1][3] = rho0_h_plus_b2*uD[1]*uD[3] + P_plus_half_b2*metric_aux->g4DD[1][3] - smallb_lower[1]*smallb_lower[3];
  Tmunu->T4[2][2] = rho0_h_plus_b2*uD[2]*uD[2] + P_plus_half_b2*metric_aux->g4DD[2][2] - smallb_lower[2]*smallb_lower[2];
  Tmunu->T4[2][3] = rho0_h_plus_b2*uD[2]*uD[3] + P_plus_half_b2*metric_aux->g4DD[2][3] - smallb_lower[2]*smallb_lower[3];
  Tmunu->T4[3][3] = rho0_h_plus_b2*uD[3]*uD[3] + P_plus_half_b2*metric_aux->g4DD[3][3] - smallb_lower[3]*smallb_lower[3];
}
