#include <stdio.h>
#include "con2prim_header.h"

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

static inline void compute_smallba_b2_and_u_i_over_u0_psi4(const metric_quantities *restrict metric, const primitive_quantities *restrict prims,
                                             const double u0L, const double ONE_OVER_LAPSE_SQRT_4PI, double *restrict u_x_over_u0_psi4,
                                             double *restrict u_y_over_u0_psi4, double *restrict u_z_over_u0_psi4, double *restrict smallb) {

  double u_over_u0_psi4[3];
   u_over_u0_psi4[0] = *u_x_over_u0_psi4;
   u_over_u0_psi4[1] = *u_y_over_u0_psi4;
   u_over_u0_psi4[2] = *u_z_over_u0_psi4;

  // NOW COMPUTE b^{\mu} and b^2 = b^{\mu} b^{\nu} g_{\mu \nu}
  double ONE_OVER_U0 = 1.0/u0L;
  double shiftx_plus_vx = (metric->betax+prims->vx);
  double shifty_plus_vy = (metric->betay+prims->vy);
  double shiftz_plus_vz = (metric->betaz+prims->vz);

  // Eq. 56 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  //  u_i = gamma_{ij} u^0 (v^j + beta^j), gamma_{ij} is the physical metric, and gamma_{ij} = Psi4 * METRIC[Gij], since METRIC[Gij] is the conformal metric.
  u_over_u0_psi4[0] =  (metric->adm_gxx*shiftx_plus_vx + metric->adm_gxy*shifty_plus_vy + metric->adm_gxz*shiftz_plus_vz)/metric->psi4;
  u_over_u0_psi4[1] =  (metric->adm_gxy*shiftx_plus_vx + metric->adm_gyy*shifty_plus_vy + metric->adm_gyz*shiftz_plus_vz)/metric->psi4;
  u_over_u0_psi4[2] =  (metric->adm_gxz*shiftx_plus_vx + metric->adm_gyz*shifty_plus_vy + metric->adm_gzz*shiftz_plus_vz)/metric->psi4;

  // Eqs. 23 and 31 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  //   Compute alpha sqrt(4 pi) b^t = u_i B^i
  double alpha_sqrt_4pi_bt = ( u_over_u0_psi4[0]*prims->Bx + u_over_u0_psi4[1]*prims->By + u_over_u0_psi4[2]*prims->Bz ) * metric->psi4*u0L;

  // Eq. 24 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  // b^i = B^i_u / sqrt(4 pi)
  // b^i = ( B^i/alpha + B^0_u u^i ) / ( u^0 sqrt(4 pi) )
  // b^i = ( B^i/alpha +  sqrt(4 pi) b^t u^i ) / ( u^0 sqrt(4 pi) )
  // b^i = ( B^i +  alpha sqrt(4 pi) b^t u^i ) / ( alpha u^0 sqrt(4 pi) )
  // b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t u^i/u^0 ) / ( alpha sqrt(4 pi) )
  // b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t v^i ) / ( alpha sqrt(4 pi) )
  smallb[SMALLBX] = (prims->Bx*ONE_OVER_U0 + prims->vx*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  smallb[SMALLBY] = (prims->By*ONE_OVER_U0 + prims->vy*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  smallb[SMALLBZ] = (prims->Bz*ONE_OVER_U0 + prims->vz*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  // Eq. 23 in http://arxiv.org/pdf/astro-ph/0503420.pdf, with alpha sqrt (4 pi) b^2 = u_i B^i already computed above
  smallb[SMALLBT] = alpha_sqrt_4pi_bt * ONE_OVER_LAPSE_SQRT_4PI;

  // b^2 = g_{\mu \nu} b^{\mu} b^{\nu}
  //     = gtt bt^2 + gxx bx^2 + gyy by^2 + gzz bz^2 + 2 (gtx bt bx + gty bt by + gtz bt bz + gxy bx by + gxz bx bz + gyz by bz)
  //     = (-al^2 + gamma_{ij} betai betaj) bt^2 + b^i b^j gamma_{ij} + 2 g_{t i} b^t b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t g_{t i} b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t (gamma_{ij} beta^j) b^i
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + b^i b^j + 2 b^t beta^j b^i)
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + 2 b^t beta^j b^i + b^i b^j)
  //     = - (alpha b^t)^2 + gamma_{ij} (b^i + b^t beta^i) (b^j + b^t beta^j)
  double bx_plus_shiftx_bt = smallb[SMALLBX]+metric->betax*smallb[SMALLBT];
  double by_plus_shifty_bt = smallb[SMALLBY]+metric->betay*smallb[SMALLBT];
  double bz_plus_shiftz_bt = smallb[SMALLBZ]+metric->betaz*smallb[SMALLBT];
  smallb[SMALLB2] = -SQR(metric->lapse*smallb[SMALLBT]) +
    metric->adm_gxx*SQR(bx_plus_shiftx_bt) + metric->adm_gyy*SQR(by_plus_shifty_bt) + metric->adm_gzz*SQR(bz_plus_shiftz_bt) +
       2.0*( metric->adm_gxy*(bx_plus_shiftx_bt)*(by_plus_shifty_bt) +
             metric->adm_gxz*(bx_plus_shiftx_bt)*(bz_plus_shiftz_bt) +
             metric->adm_gyz*(by_plus_shifty_bt)*(bz_plus_shiftz_bt) );
  *u_x_over_u0_psi4 = u_over_u0_psi4[0];
  *u_y_over_u0_psi4 = u_over_u0_psi4[1];
  *u_z_over_u0_psi4 = u_over_u0_psi4[2];
  /***********************************************************/
}

void enforce_limits_on_primitives_and_recompute_conservs(const GRMHD_parameters *restrict params, const eos_parameters *restrict eos,
                                                         const metric_quantities *restrict metric, primitive_quantities *restrict prims,
                                                         conservative_quantities *restrict cons, double *TUPMUNU,double *TDNMUNU,
                                                         stress_energy *restrict Tmunu,
                                                         con2prim_diagnostics *restrict diagnostics) {

  // The goal here is to apply floors and ceilings to the primitives
  // and compute the enthalpy. Note that the derivative of the
  // pressure with respect to rho was unecessary in this function.
  //
  // This function will apply floors and ceilings and recompute:
  //
  // rho_b
  // P
  // eps
  // S  (if evolving the entropy)
  // Ye (if tabulated EOS is enabled)
  // T  (if tabulated EOS is enabled)
  prims_enforce_extrema_and_recompute(params,eos,metric,prims);

  // Now compute the enthalpy
  const double h_enthalpy = 1.0 + prims->eps + prims->press/prims->rho;

  double uUP[4];
  impose_speed_limit_output_u0(eos, metric, prims, diagnostics, &uUP[0]);

  // Compute u^i. We've already set uUP[0] in the lines above.
  uUP[1] = uUP[0]*prims->vx;
  uUP[2] = uUP[0]*prims->vy;
  uUP[3] = uUP[0]*prims->vz;

  /***************************************************************/
  // COMPUTE TUPMUNU, TDNMUNU, AND  CONSERVATIVES FROM PRIMITIVES
  /***************************************************************/
  // Compute b^{\mu}, b^2, and u_i/(u^0 Psi4)
  double ONE_OVER_LAPSE_SQRT_4PI = metric->lapseinv*ONE_OVER_SQRT_4PI;
  double u_x_over_u0_psi4,u_y_over_u0_psi4,u_z_over_u0_psi4;
  double smallb[NUMVARS_SMALLB];
  compute_smallba_b2_and_u_i_over_u0_psi4(metric, prims, uUP[0], ONE_OVER_LAPSE_SQRT_4PI, &u_x_over_u0_psi4,
                                          &u_y_over_u0_psi4, &u_z_over_u0_psi4, smallb);

  // Compute u_i; we compute u_0 below.
  double uDN[4] = { 1e200, u_x_over_u0_psi4*uUP[0]*metric->psi4,u_y_over_u0_psi4*uUP[0]*metric->psi4,u_z_over_u0_psi4*uUP[0]*metric->psi4 };

  // Precompute some useful quantities, for later:
  double alpha_sqrt_gamma=metric->lapse*metric->psi6;
  double rho0_h_plus_b2 = (prims->rho*h_enthalpy + smallb[SMALLB2]);
  double P_plus_half_b2 = (prims->press+0.5*smallb[SMALLB2]);

  int count;
  // Next compute T^{\mu \nu}. This is very useful when computing fluxes and source terms in the GRMHD evolution equations.
  // (Eq. 33 in http://arxiv.org/pdf/astro-ph/0503420.pdf):
  // T^{mn} = (rho_0 h + b^2) u^m u^n + (P + 0.5 b^2) g^{mn} - b^m b^n, where m and n both run from 0 to 3.
  count=0; for(int ii=0;ii<4;ii++) for(int jj=ii;jj<4;jj++) { TUPMUNU[count] = rho0_h_plus_b2*uUP[ii]*uUP[jj] + P_plus_half_b2*metric->g4up[ii][jj] - smallb[SMALLBT+ii]*smallb[SMALLBT+jj]; count++; }


  // Next compute T_{\mu \nu}
  // T_{mn} = (rho_0 h + b^2) u_m u_n + (P + 0.5 b^2) g_{mn} - b_m b_n, where m and n both run from 0 to 3.

  double smallb_lower[NUMVARS_SMALLB];
  // FIXME: This could be replaced by the function call
  //           lower_4vector_output_spatial_part(psi4,smallb,smallb_lower);
  // b_a = b^c g_{ac}
  for(int ii=0;ii<4;ii++) { smallb_lower[SMALLBT+ii]=0; for(int jj=0;jj<4;jj++) smallb_lower[SMALLBT+ii] += smallb[SMALLBT+jj]*metric->g4dn[ii][jj]; }

  // Compute u_0, as we've already computed u_i above.
  uDN[0]=0.0; for(int jj=0;jj<4;jj++) uDN[0] += uUP[jj]*metric->g4dn[0][jj];

  // Compute T_{\mu \nu}
  if(params->update_Tmunu) {
    count=0; for(int ii=0;ii<4;ii++) for(int jj=ii;jj<4;jj++) { TDNMUNU[count] = rho0_h_plus_b2*uDN[ii]*uDN[jj] + P_plus_half_b2*metric->g4dn[ii][jj] - smallb_lower[SMALLBT+ii]*smallb_lower[SMALLBT+jj]; count++; }
  }

  // Finally, compute conservatives:
  cons->rho = alpha_sqrt_gamma * prims->rho * uUP[0];
  cons->S_x = cons->rho*h_enthalpy*uDN[1] + alpha_sqrt_gamma*(uUP[0]*smallb[SMALLB2]*uDN[1] - smallb[SMALLBT]*smallb_lower[SMALLBX]);
  cons->S_y = cons->rho*h_enthalpy*uDN[2] + alpha_sqrt_gamma*(uUP[0]*smallb[SMALLB2]*uDN[2] - smallb[SMALLBT]*smallb_lower[SMALLBY]);
  cons->S_z = cons->rho*h_enthalpy*uDN[3] + alpha_sqrt_gamma*(uUP[0]*smallb[SMALLB2]*uDN[3] - smallb[SMALLBT]*smallb_lower[SMALLBZ]);
  // tauL = alpha^2 sqrt(gamma) T^{00} - CONSERVS[RHOSTAR]
  cons->tau =  metric->lapse*alpha_sqrt_gamma*(rho0_h_plus_b2*SQR(uUP[0]) - P_plus_half_b2*metric->lapseinv2 - SQR(smallb[SMALLBT])) - cons->rho;

  if( params->evolve_entropy ) {
    // Entropy equation evolves S_star = alpha * sqrt(gamma) * S * u^{0}
    cons->entropy = alpha_sqrt_gamma * prims->entropy * uUP[0];
  }
  if( eos->eos_type == 1 ) {
    // Tabulated EOS evolves Y_e_star = alpha * sqrt(gamma) * rho_b * Y_e * u^{0} = rho_star * Y_e
    cons->Y_e = cons->rho * prims->Y_e;
  }
    if(params->update_Tmunu) {
      int ww=0;
      Tmunu->Ttt = TDNMUNU[ww++];
      Tmunu->Ttx = TDNMUNU[ww++];
      Tmunu->Tty = TDNMUNU[ww++];
      Tmunu->Ttz = TDNMUNU[ww++];
      Tmunu->Txx = TDNMUNU[ww++];
      Tmunu->Txy = TDNMUNU[ww++];
      Tmunu->Txz = TDNMUNU[ww++];
      Tmunu->Tyy = TDNMUNU[ww++];
      Tmunu->Tyz = TDNMUNU[ww++];
      Tmunu->Tzz = TDNMUNU[ww  ];
    }
}
