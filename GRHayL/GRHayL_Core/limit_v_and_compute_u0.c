#include "con2prim.h"

/*
 * Function     : ghl_limit_v_and_compute_u0()
 * Description  : Applies speed limit to v^i and computes u^0
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_limit_v_and_compute_u0
*/

int ghl_limit_v_and_compute_u0(
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      ghl_primitive_quantities *restrict prims) {

  int speed_limited = 0;
  // Begin by computing vel^2, where vel^{i} = (v^{i} + beta^{i})/alpha
  const double velU[3] = {(prims->vU[0] + ADM_metric->betaU[0])*ADM_metric->lapseinv,
                          (prims->vU[1] + ADM_metric->betaU[1])*ADM_metric->lapseinv,
                          (prims->vU[2] + ADM_metric->betaU[2])*ADM_metric->lapseinv};
  double vsq           = ghl_compute_vec2_from_vec3D(ADM_metric->gammaDD, velU);
  const double vsq_max = 1 - 1/(eos->W_max*eos->W_max);

  /*** Limit velocity to GAMMA_SPEED_LIMIT ***/
  if(vsq > vsq_max) {
    const double fac     = ADM_metric->lapse * sqrt(vsq_max / vsq);
    // Recompute v^{i} = vel^{i}_{rescaled} * alpha - beta^{i}
    prims->vU[0]  = velU[0]*fac - ADM_metric->betaU[0];
    prims->vU[1]  = velU[1]*fac - ADM_metric->betaU[1];
    prims->vU[2]  = velU[2]*fac - ADM_metric->betaU[2];
    vsq           = vsq_max;
    speed_limited = 1;
  }

  // A = 1.0-one_minus_one_over_alpha_u0_squared = 1-(1-1/(al u0)^2) = 1/(al u0)^2
  // 1/sqrt(A) = al u0
  //double alpha_u0_minus_one = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared)-1.0;
  //u0_out          = (alpha_u0_minus_one + 1.0)*lapseinv;
  // const double alpha_u0 = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared);
  const double W = 1/sqrt(1-vsq);
  prims->u0 = W*ADM_metric->lapseinv;
  if(isnan(prims->u0)) {
    ghl_error("*********************************************\n"
                 "Found nan while computing u^{0}\nMetric: %e %e %e %e %e %e\n"
                 "Lapse/shift: %e (=1/%e) / %e %e %e\nVelocities : %e %e %e\n"
                 "*********************************************\n",
                 ADM_metric->gammaDD[0][0], ADM_metric->gammaDD[0][1], ADM_metric->gammaDD[0][2],
                 ADM_metric->gammaDD[1][1], ADM_metric->gammaDD[1][2], ADM_metric->gammaDD[2][2],
                 ADM_metric->lapse, ADM_metric->lapseinv,
                 ADM_metric->betaU[0], ADM_metric->betaU[1], ADM_metric->betaU[2],
                 prims->vU[0], prims->vU[1], prims->vU[2]);
  }
  return speed_limited;
}
