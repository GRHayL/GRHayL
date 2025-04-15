#include "ghl_con2prim.h"

/**
 * @ingroup GRHayL_Core
 * @brief Applies speed limit to \f$ v^i \f$ and computes \f$ u^0 \f$.
 *
 * @details
 * This function applies a speed limit to \f$ v^i \f$ based on
 * ghl_parameters::max_lorenz_factor and computes ghl_primitive_quantities::u0.
 *
 * Let the maximum Lorenz factor be \f$ W_\mathrm{max} \equiv 1 - \left( 1 - u^\mu u_\mu \right)^{-1/2} \f$.
 * Then, we want to compute \f$ u0 \f$ and enforce the speed limit.
 * However, @grhayl uses the velocity used by `IllinoisGRMHD`, which is
 * defined as \f$ v^i = u^i/u^0 \f$ where \f$ u^\mu \f$ is the 4-velocity.
 * Note that this is different from another common choice of velocity,
 * the Valencia 3-velocity \f$ \tilde{u}^i = (v_V^i + \beta^i)/\alpha \f$.
 * We then calculate
 * \f[
 * \begin{align}
 * \Omega &= \tilde{u}^i \tilde{u}_i = \gamma\_{i j} \alpha^{-2} (v^i + \beta^i)(v^j + \beta^j) \\
 *        &= \gamma_{ij} \left( u^0 \alpha \right)^{-2} \gamma^{ik} u_k \gamma^{jl} u_l \\
 *        &= \gamma_{j l} \left( u^0 \alpha \right)^{-2} u^j u^l \\
 *        &= \left( u^0 \alpha \right)^{-2} (u^0 \alpha - 1) \\
 *        &= 1 - \left( u^0 \alpha \right)^{-2}
 * \end{align}
 * \f]
 *
 * where step 2 comes from Eq. 53 of [Duez _et al._](https://arxiv.org/abs/astro-ph/0503420)
 * and step 4 comes from Eq. 56 of the same paper. The quantity \f$ \Omega \f$
 * is conveniently related to both the Lorenz factor and \f$ u^0 \f$, giving both
 * of the quantities we want. We first check that
 * \f[
 * \Omega < 1 - W_\mathrm{max}^{-2}
 * \f]
 *
 * If it is not, then we trigger a speed limit by rescaling \f$ \tilde{u}^i \f$ by
 * \f[
 * \frac{1 - W_\mathrm{max}^{-2}}{\Omega}
 * \f]
 *
 * which is guaranteed to give a velocity that has a sufficiently small Lorenz
 * factor. We recompute \f$ v^i \f$ from \f$ \tilde{u}^i \f$ and then compute
 * \f[
 * u^0 = \alpha^{-1} \left( 1 - \Omega \right)^{-1/2}
 * \f]
 *
 * @param[in] params:         pointer to ghl_parameters struct
 *
 * @param[in] ADM_metric:     pointer to ghl_metric_quantities struct with ADM metric data
 *
 * @param[in] prims:          pointer to ghl_primitive_quantities struct
 *
 * @param[out] speed_limited: whether speed limiter was triggered (True) or not (False)
 *
 * @returns @ref ghl_error_codes_t: if \f$ u^0 \f$ is singular, this error code will be non-zero
 */

ghl_error_codes_t ghl_limit_v_and_compute_u0(
      const ghl_parameters *restrict params,
      const ghl_metric_quantities *restrict ADM_metric,
      ghl_primitive_quantities *restrict prims,
      bool *restrict speed_limited) {

  // Derivation of first equation:
  // \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2
  //   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using Eq. 53 of arXiv:astro-ph/0503420
  //   = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <- Since \gamma_{ij} \gamma^{ik} = \delta^k_j
  //   = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2 - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420
  //   = 1 - 1/(u^0 \alpha)^2 <= 1
  const double utU[3] = {prims->vU[0] + ADM_metric->betaU[0], prims->vU[1] + ADM_metric->betaU[1], prims->vU[2] + ADM_metric->betaU[2]};
  double one_minus_one_over_alpha_u0_squared = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, utU)*ADM_metric->lapseinv2;

  /*** Limit velocity to GAMMA_SPEED_LIMIT ***/
  const double one_minus_one_over_W_max_squared = 1.0 - params->inv_sq_max_Lorentz_factor; // 1 - W_max^{-2}
  if(one_minus_one_over_alpha_u0_squared > one_minus_one_over_W_max_squared) {
    const double correction_fac = sqrt(one_minus_one_over_W_max_squared/one_minus_one_over_alpha_u0_squared);
    prims->vU[0] = utU[0]*correction_fac - ADM_metric->betaU[0];
    prims->vU[1] = utU[1]*correction_fac - ADM_metric->betaU[1];
    prims->vU[2] = utU[2]*correction_fac - ADM_metric->betaU[2];
    one_minus_one_over_alpha_u0_squared = one_minus_one_over_W_max_squared;
    *speed_limited |= true;
  }

  // A = 1.0-one_minus_one_over_alpha_u0_squared = 1-(1-1/(al u0)^2) = 1/(al u0)^2
  // 1/sqrt(A) = al u0
  //double alpha_u0_minus_one = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared)-1.0;
  //u0_out          = (alpha_u0_minus_one + 1.0)*lapseinv;
  const double alpha_u0 = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared);
  prims->u0 = alpha_u0*ADM_metric->lapseinv;
  if(isnan(prims->u0))
    return ghl_error_u0_singular;

  return ghl_success;
}
