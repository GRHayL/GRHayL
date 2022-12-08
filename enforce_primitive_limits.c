#include <stdio.h>
#include "con2prim_header.h"
// These were taken from Leo Werneck's IllinoisGRMHD code for con2prim to work,
// but they should actually be provided by the EOS code.

#include "con2prim_header.h"
#include "EOS/Hybrid/EOS_hybrid_header.h"
#include <stdio.h>

static inline void impose_speed_limit_output_u0(const eos_parameters *restrict eos, const metric_quantities *restrict metric,
                                         primitive_quantities *restrict prims, con2prim_diagnostics *restrict diagnostics,
                                         double *restrict u0_out) {

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
  if(isnan(alpha_u0*metric->lapseinv)) {
    printf("*********************************************\n");
    printf("Metric/psi4: %e %e %e %e %e %e / %e\n", metric->adm_gxx, metric->adm_gxy, metric->adm_gxz, metric->adm_gyy, metric->adm_gyz, metric->adm_gzz, metric->psi4);
    printf("Lapse/shift: %e (=1/%e) / %e %e %e\n",metric->lapse, metric->lapseinv, metric->betax, metric->betay, metric->betaz);
    printf("Velocities : %e %e %e\n", prims->vx, prims->vx, prims->vz);
    printf("Found nan while computing u^{0} in function %s (file: %s)\n",__func__,__FILE__);
    printf("*********************************************\n");
    diagnostics->nan_found++;
  }
  *u0_out = alpha_u0*metric->lapseinv;
}

void enforce_primitive_limits_and_output_u0(const GRMHD_parameters *restrict params, const eos_parameters *restrict eos,
                                            const metric_quantities *restrict metric, primitive_quantities *restrict prims,
                                            double *restrict u0, con2prim_diagnostics *restrict diagnostics) {

  // The density floor and ceiling is always applied
  prims->rho = MIN(MAX(prims->rho,eos->rho_atm),eos->rho_max);

  // Hybrid EOS specific floors and ceilings
  if( eos->eos_type == 0 ) {
    // Pressure and epsilon must be recomputed
    // Compute P and eps
    double prs_cold = 0.0;
    double eps_cold = 0.0;
    compute_P_cold_and_eps_cold(eos, prims->rho, &prs_cold, &eps_cold);
    // Set P_min and P_max
    double P_min = 0.9*prs_cold;
    double P_max = 100.0*prs_cold;
    // Set Psi6
    // Adjust P_max based on Psi6
    if(metric->psi6 > params->psi6threshold) P_max = 1e5*prs_cold; // <-- better than 10.
    // Now apply floors and ceilings to P
    if(prims->press<P_min) prims->press = P_min;
    // Finally, perform the last check
    if((prims->rho < 100.0*eos->rho_atm || metric->psi6 > params->psi6threshold) && prims->press>P_max) {
      prims->press = P_max;
    }
    // Now compute eps
    prims->eps = eps_cold + (prims->press-prs_cold)/(eos->Gamma_th-1.0)/prims->rho;
    // If needed, recompute the entropy function
    if( params->evolve_entropy ) compute_entropy_function(eos, prims->rho, prims->press, &prims->entropy);

  // Tabulated EOS specific floors and ceilings
  } else if( eos->eos_type==1 ) {
    printf("No tabulated EOS support yet! Sorry!");
  //  // Apply floors and ceilings to Y_e and T
  //  const double xye   = MIN(MAX(prims->Y_e,eos->Ye_min),eos->Ye_atm);
  //  const double xtemp = MIN(MAX(prims->temp,eos->temp_atm ),eos->temp_max );

  //  // Additional variables used for the EOS call
  //  const double xrho  = prims->rho;
  //  double xprs        = 0.0;
  //  double xeps        = 0.0;
  //  double xeps        = 0.0;
  //  double xent        = 0.0;
  //  WVU_EOS_P_eps_and_S_from_rho_Ye_T(xrho,xye,xtemp, &xprs,&xeps,&xent);
  //  // Now update the primitives (rho has already been set)
  //  prims->Y_e = xye;
  //  prims->temp = xtemp;
  //  prims->press = xprs;
  //  prims->eps = xeps;
  //  prims->entropy = xent;
  }

  impose_speed_limit_output_u0(eos, metric, prims, diagnostics, u0);
}
