#include "con2prim.h"

static const int SMALLBT=0,SMALLBX=1,SMALLBY=2,SMALLBZ=3,SMALLB2=4,NUMVARS_SMALLB=5;

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

// TODO: really shouldn't set prims->eps in here
void compute_conservs_and_Tmunu(const GRHayL_parameters *restrict params,
                                const eos_parameters *restrict eos,
                                const metric_quantities *restrict metric,
                                primitive_quantities *restrict prims,
                                const double u0,
                                conservative_quantities *restrict cons,
                                stress_energy *restrict Tmunu) {


  double prs_cold = 0.0;
  double eps_cold = 0.0;
  eos->hybrid_compute_P_cold_and_eps_cold(eos, prims->rho, &prs_cold, &eps_cold);
  prims->eps = eps_cold + (prims->press-prs_cold)/(eos->Gamma_th-1.0)/prims->rho;

  // Now compute the enthalpy
  const double h_enthalpy = 1.0 + prims->eps + prims->press/prims->rho;

  double uUP[4];

  // Compute u^i. u^0 is provided to the function.
  uUP[0] = u0;
  uUP[1] = uUP[0]*prims->vx;
  uUP[2] = uUP[0]*prims->vy;
  uUP[3] = uUP[0]*prims->vz;

  /***************************************************************/
  //     COMPUTE TDNMUNU AND  CONSERVATIVES FROM PRIMITIVES      //
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


  double smallb_lower[NUMVARS_SMALLB];
  // FIXME: This could be replaced by the function call
  //           lower_4vector_output_spatial_part(psi4,smallb,smallb_lower);
  // b_a = b^c g_{ac}
  for(int ii=0;ii<4;ii++) { smallb_lower[SMALLBT+ii]=0; for(int jj=0;jj<4;jj++) smallb_lower[SMALLBT+ii] += smallb[SMALLBT+jj]*metric->g4dn[ii][jj]; }

  // Compute u_0, as we've already computed u_i above.
  uDN[0]=0.0; for(int jj=0;jj<4;jj++) uDN[0] += uUP[jj]*metric->g4dn[0][jj];

  // Compute conservatives:
  cons->rho = alpha_sqrt_gamma * prims->rho * uUP[0];
  cons->S_x = cons->rho*h_enthalpy*uDN[1] + alpha_sqrt_gamma*(uUP[0]*smallb[SMALLB2]*uDN[1] - smallb[SMALLBT]*smallb_lower[SMALLBX]);
  cons->S_y = cons->rho*h_enthalpy*uDN[2] + alpha_sqrt_gamma*(uUP[0]*smallb[SMALLB2]*uDN[2] - smallb[SMALLBT]*smallb_lower[SMALLBY]);
  cons->S_z = cons->rho*h_enthalpy*uDN[3] + alpha_sqrt_gamma*(uUP[0]*smallb[SMALLB2]*uDN[3] - smallb[SMALLBT]*smallb_lower[SMALLBZ]);
  // tauL = alpha^2 sqrt(gamma) T^{00} - CONSERVS[RHOSTAR]
  cons->tau =  metric->lapse*alpha_sqrt_gamma*(rho0_h_plus_b2*SQR(uUP[0]) - P_plus_half_b2*metric->lapseinv2 - SQR(smallb[SMALLBT])) - cons->rho;
  // Entropy equation evolves S_star = alpha * sqrt(gamma) * S * u^{0}
  cons->entropy = alpha_sqrt_gamma * prims->entropy * uUP[0];
  // Tabulated EOS evolves Y_e_star = alpha * sqrt(gamma) * rho_b * Y_e * u^{0} = rho_star * Y_e
  cons->Y_e = cons->rho * prims->Y_e;

  // Finally, compute T_{\mu \nu}
  // T_{mn} = (rho_0 h + b^2) u_m u_n + (P + 0.5 b^2) g_{mn} - b_m b_n, where m and n both run from 0 to 3.
  if(params->update_Tmunu) {
    Tmunu->Ttt = rho0_h_plus_b2*uDN[0]*uDN[0] + P_plus_half_b2*metric->g4dn[0][0] - smallb_lower[SMALLBT+0]*smallb_lower[SMALLBT+0];;
    Tmunu->Ttx = rho0_h_plus_b2*uDN[0]*uDN[1] + P_plus_half_b2*metric->g4dn[0][1] - smallb_lower[SMALLBT+0]*smallb_lower[SMALLBT+1];;
    Tmunu->Tty = rho0_h_plus_b2*uDN[0]*uDN[2] + P_plus_half_b2*metric->g4dn[0][2] - smallb_lower[SMALLBT+0]*smallb_lower[SMALLBT+2];;
    Tmunu->Ttz = rho0_h_plus_b2*uDN[0]*uDN[3] + P_plus_half_b2*metric->g4dn[0][3] - smallb_lower[SMALLBT+0]*smallb_lower[SMALLBT+3];;
    Tmunu->Txx = rho0_h_plus_b2*uDN[1]*uDN[1] + P_plus_half_b2*metric->g4dn[1][1] - smallb_lower[SMALLBT+1]*smallb_lower[SMALLBT+1];;
    Tmunu->Txy = rho0_h_plus_b2*uDN[1]*uDN[2] + P_plus_half_b2*metric->g4dn[1][2] - smallb_lower[SMALLBT+1]*smallb_lower[SMALLBT+2];;
    Tmunu->Txz = rho0_h_plus_b2*uDN[1]*uDN[3] + P_plus_half_b2*metric->g4dn[1][3] - smallb_lower[SMALLBT+1]*smallb_lower[SMALLBT+3];;
    Tmunu->Tyy = rho0_h_plus_b2*uDN[2]*uDN[2] + P_plus_half_b2*metric->g4dn[2][2] - smallb_lower[SMALLBT+2]*smallb_lower[SMALLBT+2];;
    Tmunu->Tyz = rho0_h_plus_b2*uDN[2]*uDN[3] + P_plus_half_b2*metric->g4dn[2][3] - smallb_lower[SMALLBT+2]*smallb_lower[SMALLBT+3];;
    Tmunu->Tzz = rho0_h_plus_b2*uDN[3]*uDN[3] + P_plus_half_b2*metric->g4dn[3][3] - smallb_lower[SMALLBT+3]*smallb_lower[SMALLBT+3];;
//    double TDNMUNU[10];
//    int ww=0;
//    for(int ii=0;ii<4;ii++) for(int jj=ii;jj<4;jj++) { TDNMUNU[ww] = rho0_h_plus_b2*uDN[ii]*uDN[jj] + P_plus_half_b2*metric->g4dn[ii][jj] - smallb_lower[SMALLBT+ii]*smallb_lower[SMALLBT+jj]; ww++; }
//    ww=0;
//    Tmunu->Ttt = TDNMUNU[ww++];
//    Tmunu->Ttx = TDNMUNU[ww++];
//    Tmunu->Tty = TDNMUNU[ww++];
//    Tmunu->Ttz = TDNMUNU[ww++];
//    Tmunu->Txx = TDNMUNU[ww++];
//    Tmunu->Txy = TDNMUNU[ww++];
//    Tmunu->Txz = TDNMUNU[ww++];
//    Tmunu->Tyy = TDNMUNU[ww++];
//    Tmunu->Tyz = TDNMUNU[ww++];
//    Tmunu->Tzz = TDNMUNU[ww  ];
  }
}
