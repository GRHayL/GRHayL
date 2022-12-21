#include "con2prim.h"
#include <stdio.h>

void limit_v_and_output_u0(const eos_parameters *restrict eos,
                           const metric_quantities *restrict metric,
                           primitive_quantities *restrict prims,
                           double *restrict u0_out,
                           con2prim_diagnostics *restrict diagnostics) {

  // Derivation of first equation:
  // \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2
  //   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using Eq. 53 of arXiv:astro-ph/0503420
  //   = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <- Since \gamma_{ij} \gamma^{ik} = \delta^k_j
  //   = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2 - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420
  //   = 1 - 1/(u^0 \alpha)^2 <= 1
  double one_minus_one_over_alpha_u0_squared = (metric->adm_gxx * SQR(prims->vx + metric->betax) +
                                                2.0*metric->adm_gxy*(prims->vx + metric->betax)*(prims->vy + metric->betay) +
                                                2.0*metric->adm_gxz*(prims->vx + metric->betax)*(prims->vz + metric->betaz) +
                                                metric->adm_gyy * SQR(prims->vy + metric->betay) +
                                                2.0*metric->adm_gyz*(prims->vy + metric->betay)*(prims->vz + metric->betaz) +
                                                metric->adm_gzz * SQR(prims->vz + metric->betaz) )*metric->lapseinv2;

  /*** Limit velocity to GAMMA_SPEED_LIMIT ***/
  const double one_minus_one_over_W_max_squared = 1.0-1.0/SQR(eos->W_max); // 1 - W_max^{-2}
  if(one_minus_one_over_alpha_u0_squared > one_minus_one_over_W_max_squared) {
    double correction_fac = sqrt(one_minus_one_over_W_max_squared/one_minus_one_over_alpha_u0_squared);
    prims->vx = (prims->vx + metric->betax)*correction_fac - metric->betax;
    prims->vy = (prims->vy + metric->betay)*correction_fac - metric->betay;
    prims->vz = (prims->vz + metric->betaz)*correction_fac - metric->betaz;
    one_minus_one_over_alpha_u0_squared = one_minus_one_over_W_max_squared;
    diagnostics->failure_checker+=1000;
//printf("speed limited 2\n");
  }

  // A = 1.0-one_minus_one_over_alpha_u0_squared = 1-(1-1/(al u0)^2) = 1/(al u0)^2
  // 1/sqrt(A) = al u0
  //double alpha_u0_minus_one = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared)-1.0;
  //u0_out          = (alpha_u0_minus_one + 1.0)*metric.lapseinv;
  double alpha_u0 = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared);
  *u0_out = alpha_u0*metric->lapseinv;
  if(isnan(*u0_out)) {
    printf("*********************************************\n");
    printf("Metric/psi4: %e %e %e %e %e %e / %e\n", metric->adm_gxx, metric->adm_gxy, metric->adm_gxz, metric->adm_gyy, metric->adm_gyz, metric->adm_gzz, metric->psi4);
    printf("Lapse/shift: %e (=1/%e) / %e %e %e\n",metric->lapse, metric->lapseinv, metric->betax, metric->betay, metric->betaz);
    printf("Velocities : %e %e %e\n", prims->vx, prims->vx, prims->vz);
    printf("Found nan while computing u^{0} in function %s (file: %s)\n",__func__,__FILE__);
    printf("*********************************************\n");
    diagnostics->nan_found++;
  }
}
