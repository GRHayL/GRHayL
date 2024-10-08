#include "ghl_con2prim.h"

/*
 * Function     : ghl_limit_v_and_compute_u0()
 * Description  : Applies speed limit to v^i and computes u^0
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_limit_v_and_compute_u0
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
