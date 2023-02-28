#include "con2prim.h"

/* Function    : compute_conservs_and_Tmunu()
 * Description : Computes the conservatives from the given primitives and
 *               computes Tmunu if params->update_Tmunu is true.
 *
 * Inputs      : params         - GRHayL_parameters struct with parameters
 *                                for the simulation
 *             : metric         - metric_quantities struct with data for
 *                                the gridpoint of interest
 *             : prims          - primitive_quantities struct with data
 *                                for the gridpoint of interest
 *
 * Outputs     : cons           - returns computed conservative values
 *             : Tmunu          - returns computed stress-energy tensor
 *
 */

void compute_conservs_and_Tmunu(const GRHayL_parameters *restrict params,
                                const metric_quantities *restrict metric,
                                const primitive_quantities *restrict prims,
                                conservative_quantities *restrict cons,
                                stress_energy *restrict Tmunu) {

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

  printf("----------------------------\n");
  printf("rho = %.15e\n", prims->rho);
  printf("Y_e = %.15e\n", prims->Y_e);
  printf(" T  = %.15e\n", prims->temperature);
  printf(" P  = %.15e\n", prims->press);
  printf("eps = %.15e\n", prims->eps);
  printf("v^x = %.15e\n", prims->vx);
  printf("v^y = %.15e\n", prims->vy);
  printf("v^z = %.15e\n", prims->vz);
  printf("B^x = %.15e\n", prims->Bx);
  printf("B^y = %.15e\n", prims->By);
  printf("B^z = %.15e\n", prims->Bz);
  printf(" h  = %.15e\n", h_enthalpy);
  printf(" ~D = %.15e\n", cons->rho);
  printf("~Sx = %.15e\n", cons->S_x);
  printf("~Sy = %.15e\n", cons->S_y);
  printf("~Sz = %.15e\n", cons->S_z);
  printf("~tau= %.15e\n", cons->tau);
  printf("~DYe= %.15e\n", cons->Y_e);
  printf("----------------------------\n");

  // Finally, compute T_{\mu \nu}
  // T_{mn} = (rho_0 h + b^2) u_m u_n + (P + 0.5 b^2) g_{mn} - b_m b_n, where m and n both run from 0 to 3.
  // We don't use the GRHayL-provided function for computing T_{\mu \nu} because we can reuse a lot of precomputed quantities
  if(params->update_Tmunu) {
    Tmunu->Ttt = rho0_h_plus_b2*uDN[0]*uDN[0] + P_plus_half_b2*metric->g4dn[0][0] - smallb_lower[0]*smallb_lower[0];
    Tmunu->Ttx = rho0_h_plus_b2*uDN[0]*uDN[1] + P_plus_half_b2*metric->g4dn[0][1] - smallb_lower[0]*smallb_lower[1];
    Tmunu->Tty = rho0_h_plus_b2*uDN[0]*uDN[2] + P_plus_half_b2*metric->g4dn[0][2] - smallb_lower[0]*smallb_lower[2];
    Tmunu->Ttz = rho0_h_plus_b2*uDN[0]*uDN[3] + P_plus_half_b2*metric->g4dn[0][3] - smallb_lower[0]*smallb_lower[3];
    Tmunu->Txx = rho0_h_plus_b2*uDN[1]*uDN[1] + P_plus_half_b2*metric->g4dn[1][1] - smallb_lower[1]*smallb_lower[1];
    Tmunu->Txy = rho0_h_plus_b2*uDN[1]*uDN[2] + P_plus_half_b2*metric->g4dn[1][2] - smallb_lower[1]*smallb_lower[2];
    Tmunu->Txz = rho0_h_plus_b2*uDN[1]*uDN[3] + P_plus_half_b2*metric->g4dn[1][3] - smallb_lower[1]*smallb_lower[3];
    Tmunu->Tyy = rho0_h_plus_b2*uDN[2]*uDN[2] + P_plus_half_b2*metric->g4dn[2][2] - smallb_lower[2]*smallb_lower[2];
    Tmunu->Tyz = rho0_h_plus_b2*uDN[2]*uDN[3] + P_plus_half_b2*metric->g4dn[2][3] - smallb_lower[2]*smallb_lower[3];
    Tmunu->Tzz = rho0_h_plus_b2*uDN[3]*uDN[3] + P_plus_half_b2*metric->g4dn[3][3] - smallb_lower[3]*smallb_lower[3];
  }
}
