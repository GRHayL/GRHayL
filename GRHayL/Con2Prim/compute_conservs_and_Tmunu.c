#include "ghl_con2prim.h"

/*
 * Function     : ghl_compute_conservs_and_Tmunu()
 * Description  : Computes the conservatives and T_munu from the given primitives
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_compute_conservs_and_Tmunu
*/
GHL_DEVICE
void ghl_compute_conservs_and_Tmunu(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons,
      ghl_stress_energy *restrict Tmunu) {

  // First compute the enthalpy
  const double h_enthalpy = 1.0 + prims->eps + prims->press/prims->rho;

  // Compute u^i. u^0 is provided to the function.
  const double uU[4] = {prims->u0,
                        prims->u0*prims->vU[0],
                        prims->u0*prims->vU[1],
                        prims->u0*prims->vU[2]};

  // Compute u_\alpha
  double uD[4];
  ghl_raise_lower_vector_4D(metric_aux->g4DD, uU, uD);

  /***************************************************************/
  //     COMPUTE TDNMUNU AND CONSERVATIVES FROM PRIMITIVES       //
  /***************************************************************/
  // Compute b^{\mu} and b^2
  double smallb[4], smallb2;
  ghl_compute_smallb_and_b2(ADM_metric, prims, uD, smallb, &smallb2);

  // Precompute some useful quantities, for later:
  const double alpha_sqrt_gamma = ADM_metric->lapse*ADM_metric->sqrt_detgamma;
  const double rho0_h_plus_b2 = (prims->rho*h_enthalpy + smallb2);
  const double P_plus_half_b2 = (prims->press+0.5*smallb2);

  double smallb_lower[4];
  ghl_raise_lower_vector_4D(metric_aux->g4DD, smallb, smallb_lower);

  // Compute conservatives:
  cons->rho = alpha_sqrt_gamma * prims->rho * uU[0];
  cons->SD[0] = cons->rho*h_enthalpy*uD[1] + alpha_sqrt_gamma*(uU[0]*smallb2*uD[1] - smallb[0]*smallb_lower[1]);
  cons->SD[1] = cons->rho*h_enthalpy*uD[2] + alpha_sqrt_gamma*(uU[0]*smallb2*uD[2] - smallb[0]*smallb_lower[2]);
  cons->SD[2] = cons->rho*h_enthalpy*uD[3] + alpha_sqrt_gamma*(uU[0]*smallb2*uD[3] - smallb[0]*smallb_lower[3]);
  // tau = alpha^2 sqrt(gamma) T^{00} - rho_star
  cons->tau = ADM_metric->lapse*alpha_sqrt_gamma*(rho0_h_plus_b2*SQR(uU[0]) - P_plus_half_b2*ADM_metric->lapseinv2 - SQR(smallb[0])) - cons->rho;
  // Entropy equation evolves S_star = alpha * sqrt(gamma) * S * u^{0}
  cons->entropy = alpha_sqrt_gamma * prims->entropy * uU[0];
  // Tabulated EOS evolves Y_e_star = alpha * sqrt(gamma) * rho_b * Y_e * u^{0} = rho_star * Y_e
  cons->Y_e = cons->rho * prims->Y_e;

  // Finally, compute T_{\mu \nu}
  // T_{mn} = (rho_0 h + b^2) u_m u_n + (P + 0.5 b^2) g_{mn} - b_m b_n, where m and n both run from 0 to 3.
  // We don't use the GRHayL-provided function for computing T_{\mu \nu} because we can reuse a lot of precomputed quantities
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
