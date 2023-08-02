#include "con2prim.h"

/*
 * Function     : ghl_compute_conservs()
 * Description  : Computes the conservatives from the given primitives
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_compute_conservs
*/

void ghl_compute_conservs(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_conservative_quantities *restrict cons) {

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
  //           COMPUTE CONSERVATIVES FROM PRIMITIVES             //
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
  cons->tau =  ADM_metric->lapse*alpha_sqrt_gamma*(rho0_h_plus_b2*SQR(uU[0]) - P_plus_half_b2*ADM_metric->lapseinv2 - SQR(smallb[0])) - cons->rho;
  // Entropy equation evolves S_star = alpha * sqrt(gamma) * S * u^{0}
  cons->entropy = alpha_sqrt_gamma * prims->entropy * uU[0];
  // Tabulated EOS evolves Y_e_star = alpha * sqrt(gamma) * rho_b * Y_e * u^{0} = rho_star * Y_e
  cons->Y_e = cons->rho * prims->Y_e;
}
