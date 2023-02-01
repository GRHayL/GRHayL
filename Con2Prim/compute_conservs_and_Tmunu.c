#include "con2prim.h"

void compute_u_i_over_u0_psi4(const metric_quantities *restrict metric, const primitive_quantities *restrict prims,
                              double *restrict u_x_over_u0_psi4, double *restrict u_y_over_u0_psi4,
                              double *restrict u_z_over_u0_psi4);
void compute_smallb_and_b2(const metric_quantities *restrict metric, const primitive_quantities *restrict prims,
                           const double u0L, const double u_x_over_u0_psi4, const double u_y_over_u0_psi4,
                           const double u_z_over_u0_psi4, double *restrict smallb, double *restrict smallb2);

/* Function    : compute_conservs_and_Tmunu()
 * Description : This function computes the conservatives from the given primitives.
 *               It also computes Tmunu if params->update_Tmunu is true.
 *
 * Inputs      : params         - GRHayL_parameters struct with parameters
 *                                for the simulation
 *             : eos            - eos_parameters struct with data for the
 *                                EOS of the simulation
 *             : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : prims          - primitive_quantities struct with data
 *                                for the gridpoint of interest
 *             : u0             - t component of 4-velocity
 *
 * Outputs     : cons           - computed conservative values
 *             : Tmunu          - computed stress-energy tensor
 *
 */

void compute_conservs_and_Tmunu(const GRHayL_parameters *restrict params,
                                const eos_parameters *restrict eos,
                                const metric_quantities *restrict metric,
                                const primitive_quantities *restrict prims,
                                const double u0,
                                conservative_quantities *restrict cons,
                                stress_energy *restrict Tmunu) {

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
  double u_x_over_u0_psi4,u_y_over_u0_psi4,u_z_over_u0_psi4;
  double smallb[4], smallb2;
  compute_u_i_over_u0_psi4(metric, prims, &u_x_over_u0_psi4, &u_y_over_u0_psi4, &u_z_over_u0_psi4);
  compute_smallb_and_b2(metric, prims, uUP[0], u_x_over_u0_psi4, u_y_over_u0_psi4, u_z_over_u0_psi4, smallb, &smallb2);

  // Compute u_0; we compute u_i below.
  double udn0=0.0;
  for(int jj=0; jj<4; jj++)
    udn0 += uUP[jj]*metric->g4dn[0][jj];

  // Compute u_i; we computed u_0 above.
  const double uDN[4] = { udn0, u_x_over_u0_psi4*uUP[0]*metric->psi4,u_y_over_u0_psi4*uUP[0]*metric->psi4,u_z_over_u0_psi4*uUP[0]*metric->psi4 };

  // Precompute some useful quantities, for later:
  const double alpha_sqrt_gamma=metric->lapse*metric->psi6;
  const double rho0_h_plus_b2 = (prims->rho*h_enthalpy + smallb2);
  const double P_plus_half_b2 = (prims->press+0.5*smallb2);

  double smallb_lower[4];
  // FIXME: This could be replaced by the function call
  //           lower_4vector_output_spatial_part(psi4,smallb,smallb_lower);
  // b_a = b^c g_{ac}
  for(int ii=0;ii<4;ii++) {
    smallb_lower[0+ii]=0;
    for(int jj=0;jj<4;jj++)
      smallb_lower[0+ii] += smallb[0+jj]*metric->g4dn[ii][jj];
  }

  // Compute conservatives:
  cons->rho = alpha_sqrt_gamma * prims->rho * uUP[0];
  cons->S_x = cons->rho*h_enthalpy*uDN[1] + alpha_sqrt_gamma*(uUP[0]*smallb2*uDN[1] - smallb[0]*smallb_lower[1]);
  cons->S_y = cons->rho*h_enthalpy*uDN[2] + alpha_sqrt_gamma*(uUP[0]*smallb2*uDN[2] - smallb[0]*smallb_lower[2]);
  cons->S_z = cons->rho*h_enthalpy*uDN[3] + alpha_sqrt_gamma*(uUP[0]*smallb2*uDN[3] - smallb[0]*smallb_lower[3]);
  // tauL = alpha^2 sqrt(gamma) T^{00} - CONSERVS[RHOSTAR]
  cons->tau =  metric->lapse*alpha_sqrt_gamma*(rho0_h_plus_b2*SQR(uUP[0]) - P_plus_half_b2*metric->lapseinv2 - SQR(smallb[0])) - cons->rho;
  // Entropy equation evolves S_star = alpha * sqrt(gamma) * S * u^{0}
  cons->entropy = alpha_sqrt_gamma * prims->entropy * uUP[0];
  // Tabulated EOS evolves Y_e_star = alpha * sqrt(gamma) * rho_b * Y_e * u^{0} = rho_star * Y_e
  cons->Y_e = cons->rho * prims->Y_e;

  // Finally, compute T_{\mu \nu}
  // T_{mn} = (rho_0 h + b^2) u_m u_n + (P + 0.5 b^2) g_{mn} - b_m b_n, where m and n both run from 0 to 3.
  if(params->update_Tmunu) {
    Tmunu->Ttt = rho0_h_plus_b2*uDN[0]*uDN[0] + P_plus_half_b2*metric->g4dn[0][0] - smallb_lower[0]*smallb_lower[0];;
    Tmunu->Ttx = rho0_h_plus_b2*uDN[0]*uDN[1] + P_plus_half_b2*metric->g4dn[0][1] - smallb_lower[0]*smallb_lower[1];;
    Tmunu->Tty = rho0_h_plus_b2*uDN[0]*uDN[2] + P_plus_half_b2*metric->g4dn[0][2] - smallb_lower[0]*smallb_lower[2];;
    Tmunu->Ttz = rho0_h_plus_b2*uDN[0]*uDN[3] + P_plus_half_b2*metric->g4dn[0][3] - smallb_lower[0]*smallb_lower[3];;
    Tmunu->Txx = rho0_h_plus_b2*uDN[1]*uDN[1] + P_plus_half_b2*metric->g4dn[1][1] - smallb_lower[1]*smallb_lower[1];;
    Tmunu->Txy = rho0_h_plus_b2*uDN[1]*uDN[2] + P_plus_half_b2*metric->g4dn[1][2] - smallb_lower[1]*smallb_lower[2];;
    Tmunu->Txz = rho0_h_plus_b2*uDN[1]*uDN[3] + P_plus_half_b2*metric->g4dn[1][3] - smallb_lower[1]*smallb_lower[3];;
    Tmunu->Tyy = rho0_h_plus_b2*uDN[2]*uDN[2] + P_plus_half_b2*metric->g4dn[2][2] - smallb_lower[2]*smallb_lower[2];;
    Tmunu->Tyz = rho0_h_plus_b2*uDN[2]*uDN[3] + P_plus_half_b2*metric->g4dn[2][3] - smallb_lower[2]*smallb_lower[3];;
    Tmunu->Tzz = rho0_h_plus_b2*uDN[3]*uDN[3] + P_plus_half_b2*metric->g4dn[3][3] - smallb_lower[3]*smallb_lower[3];;
  }
}

void compute_u_i_over_u0_psi4(
      const metric_quantities *restrict metric,
      const primitive_quantities *restrict prims,
      double *restrict u_x_over_u0_psi4,
      double *restrict u_y_over_u0_psi4,
      double *restrict u_z_over_u0_psi4) {

  double shiftx_plus_vx = (metric->betax+prims->vx);
  double shifty_plus_vy = (metric->betay+prims->vy);
  double shiftz_plus_vz = (metric->betaz+prims->vz);

  // Eq. 56 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  //  u_i = gamma_{ij} u^0 (v^j + beta^j), gamma_{ij} is the physical metric, and gamma_{ij} = Psi4 * METRIC[Gij], since METRIC[Gij] is the conformal metric.
  *u_x_over_u0_psi4 =  (metric->adm_gxx*shiftx_plus_vx + metric->adm_gxy*shifty_plus_vy + metric->adm_gxz*shiftz_plus_vz)/metric->psi4;
  *u_y_over_u0_psi4 =  (metric->adm_gxy*shiftx_plus_vx + metric->adm_gyy*shifty_plus_vy + metric->adm_gyz*shiftz_plus_vz)/metric->psi4;
  *u_z_over_u0_psi4 =  (metric->adm_gxz*shiftx_plus_vx + metric->adm_gyz*shifty_plus_vy + metric->adm_gzz*shiftz_plus_vz)/metric->psi4;
}

// Computes b^{\mu} and b^2 = b^{\mu} b^{\nu} g_{\mu \nu}
void compute_smallb_and_b2(
      const metric_quantities *restrict metric,
      const primitive_quantities *restrict prims,
      const double u0L,
      const double u_x_over_u0_psi4,
      const double u_y_over_u0_psi4,
      const double u_z_over_u0_psi4,
      double *restrict smallb,
      double *restrict smallb2) {

  double ONE_OVER_LAPSE_SQRT_4PI = metric->lapseinv*ONE_OVER_SQRT_4PI;
  double ONE_OVER_U0 = 1.0/u0L;

  // Eqs. 23 and 31 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  //   Compute alpha sqrt(4 pi) b^t = u_i B^i
  double alpha_sqrt_4pi_bt = ( u_x_over_u0_psi4*prims->Bx + u_y_over_u0_psi4*prims->By + u_z_over_u0_psi4*prims->Bz ) * metric->psi4*u0L;

  // Eq. 24 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
  // b^i = B^i_u / sqrt(4 pi)
  // b^i = ( B^i/alpha + B^0_u u^i ) / ( u^0 sqrt(4 pi) )
  // b^i = ( B^i/alpha +  sqrt(4 pi) b^t u^i ) / ( u^0 sqrt(4 pi) )
  // b^i = ( B^i +  alpha sqrt(4 pi) b^t u^i ) / ( alpha u^0 sqrt(4 pi) )
  // b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t u^i/u^0 ) / ( alpha sqrt(4 pi) )
  // b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t v^i ) / ( alpha sqrt(4 pi) )
  smallb[1] = (prims->Bx*ONE_OVER_U0 + prims->vx*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  smallb[2] = (prims->By*ONE_OVER_U0 + prims->vy*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  smallb[3] = (prims->Bz*ONE_OVER_U0 + prims->vz*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
  // Eq. 23 in http://arxiv.org/pdf/astro-ph/0503420.pdf, with alpha sqrt (4 pi) b^2 = u_i B^i already computed above
  smallb[0] = alpha_sqrt_4pi_bt * ONE_OVER_LAPSE_SQRT_4PI;

  // b^2 = g_{\mu \nu} b^{\mu} b^{\nu}
  //     = gtt bt^2 + gxx bx^2 + gyy by^2 + gzz bz^2 + 2 (gtx bt bx + gty bt by + gtz bt bz + gxy bx by + gxz bx bz + gyz by bz)
  //     = (-al^2 + gamma_{ij} betai betaj) bt^2 + b^i b^j gamma_{ij} + 2 g_{t i} b^t b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t g_{t i} b^i
  //     = - (alpha b^t)^2 + (b^t)^2 gamma_{ij} beta^i beta^j + b^i b^j gamma_{ij} + 2 b^t (gamma_{ij} beta^j) b^i
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + b^i b^j + 2 b^t beta^j b^i)
  //     = - (alpha b^t)^2 + gamma_{ij} ((b^t)^2 beta^i beta^j + 2 b^t beta^j b^i + b^i b^j)
  //     = - (alpha b^t)^2 + gamma_{ij} (b^i + b^t beta^i) (b^j + b^t beta^j)
  double bx_plus_shiftx_bt = smallb[1]+metric->betax*smallb[0];
  double by_plus_shifty_bt = smallb[2]+metric->betay*smallb[0];
  double bz_plus_shiftz_bt = smallb[3]+metric->betaz*smallb[0];
  *smallb2 = -SQR(metric->lapse*smallb[0]) +
    metric->adm_gxx*SQR(bx_plus_shiftx_bt) + metric->adm_gyy*SQR(by_plus_shifty_bt) + metric->adm_gzz*SQR(bz_plus_shiftz_bt) +
       2.0*( metric->adm_gxy*(bx_plus_shiftx_bt)*(by_plus_shifty_bt) +
             metric->adm_gxz*(bx_plus_shiftx_bt)*(bz_plus_shiftz_bt) +
             metric->adm_gyz*(by_plus_shifty_bt)*(bz_plus_shiftz_bt) );
  /***********************************************************/
}

