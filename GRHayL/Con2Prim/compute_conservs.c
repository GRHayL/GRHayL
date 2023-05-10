#include "con2prim.h"

/* Function    : compute_conservs()
 * Description : Computes the conservatives from the given primitives
 *
 * Inputs      : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : prims          - primitive_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : cons           - returns computed conservative values
 *
 */

void compute_conservs(const metric_quantities *restrict metric,
                      const primitive_quantities *restrict prims,
                      conservative_quantities *restrict cons) {

  // First compute the enthalpy
  const double h_enthalpy = 1.0 + prims->eps + prims->press/prims->rho;

  double uUP[4], uDN[4];

  // Compute u^i. u^0 is provided to the function.
  uUP[0] = prims->u0;
  uUP[1] = uUP[0]*prims->vx;
  uUP[2] = uUP[0]*prims->vy;
  uUP[3] = uUP[0]*prims->vz;

  // Compute u_\alpha
  lower_vector(metric, uUP, uDN);

  /***************************************************************/
  //     COMPUTE TDNMUNU AND  CONSERVATIVES FROM PRIMITIVES      //
  /***************************************************************/
  // Compute b^{\mu} and b^2
  double smallb[4], smallb2;
  compute_smallb_and_b2(metric, prims, uDN, smallb, &smallb2);

  // Precompute some useful quantities, for later:
  const double alpha_sqrt_gamma = metric->lapse*metric->psi6;
  const double rho0_h_plus_b2 = (prims->rho*h_enthalpy + smallb2);
  const double P_plus_half_b2 = (prims->press+0.5*smallb2);

  double smallb_lower[4];
  lower_vector(metric, smallb, smallb_lower);

  // Compute conservatives:
  cons->rho = alpha_sqrt_gamma * prims->rho * uUP[0];
  cons->S_x = cons->rho*h_enthalpy*uDN[1] + alpha_sqrt_gamma*(uUP[0]*smallb2*uDN[1] - smallb[0]*smallb_lower[1]);
  cons->S_y = cons->rho*h_enthalpy*uDN[2] + alpha_sqrt_gamma*(uUP[0]*smallb2*uDN[2] - smallb[0]*smallb_lower[2]);
  cons->S_z = cons->rho*h_enthalpy*uDN[3] + alpha_sqrt_gamma*(uUP[0]*smallb2*uDN[3] - smallb[0]*smallb_lower[3]);
  // tau = alpha^2 sqrt(gamma) T^{00} - rho_star
  cons->tau =  metric->lapse*alpha_sqrt_gamma*(rho0_h_plus_b2*SQR(uUP[0]) - P_plus_half_b2*metric->lapseinv2 - SQR(smallb[0])) - cons->rho;
  // Entropy equation evolves S_star = alpha * sqrt(gamma) * S * u^{0}
  cons->entropy = alpha_sqrt_gamma * prims->entropy * uUP[0];
  // Tabulated EOS evolves Y_e_star = alpha * sqrt(gamma) * rho_b * Y_e * u^{0} = rho_star * Y_e
  cons->Y_e = cons->rho * prims->Y_e;
}
