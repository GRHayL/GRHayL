#include "GRHayL.h"

void compute_TUPmunu(
const eos_parameters *restrict eos,
const metric_quantities *restrict metric,
const primitive_quantities *restrict prims,
stress_energy *restrict Tmunu) {

  //   First set h, the enthalpy:
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
  //                     COMPUTE TDNMUNU                         //
  /***************************************************************/
  // Compute b^{\mu} and b^2
  double smallb[4], smallb2;
  compute_smallb_and_b2(metric, prims, uDN, smallb, &smallb2);

  // Precompute some useful quantities, for later:
  const double rho0_h_plus_b2 = (prims->rho*h_enthalpy + smallb2);
  const double P_plus_half_b2 = (prims->press+0.5*smallb2);

  // Next compute T^{\mu \nu}:
  // (Eq. 33 in http://arxiv.org/pdf/astro-ph/0503420.pdf):
  // T^{mn} = (rho_0 h + b^2) u^m u^n + (P + 0.5 b^2) g^{mn} - b^m b^n, where m and n both run from 0 to 3.
  Tmunu->Ttt = rho0_h_plus_b2*uUP[0]*uUP[0] + P_plus_half_b2*metric->g4up[0][0] - smallb[0]*smallb[0];
  Tmunu->Ttx = rho0_h_plus_b2*uUP[0]*uUP[1] + P_plus_half_b2*metric->g4up[0][1] - smallb[0]*smallb[1];
  Tmunu->Tty = rho0_h_plus_b2*uUP[0]*uUP[2] + P_plus_half_b2*metric->g4up[0][2] - smallb[0]*smallb[2];
  Tmunu->Ttz = rho0_h_plus_b2*uUP[0]*uUP[3] + P_plus_half_b2*metric->g4up[0][3] - smallb[0]*smallb[3];
  Tmunu->Txx = rho0_h_plus_b2*uUP[1]*uUP[1] + P_plus_half_b2*metric->g4up[1][1] - smallb[1]*smallb[1];
  Tmunu->Txy = rho0_h_plus_b2*uUP[1]*uUP[2] + P_plus_half_b2*metric->g4up[1][2] - smallb[1]*smallb[2];
  Tmunu->Txz = rho0_h_plus_b2*uUP[1]*uUP[3] + P_plus_half_b2*metric->g4up[1][3] - smallb[1]*smallb[3];
  Tmunu->Tyy = rho0_h_plus_b2*uUP[2]*uUP[2] + P_plus_half_b2*metric->g4up[2][2] - smallb[2]*smallb[2];
  Tmunu->Tyz = rho0_h_plus_b2*uUP[2]*uUP[3] + P_plus_half_b2*metric->g4up[2][3] - smallb[2]*smallb[3];
  Tmunu->Tzz = rho0_h_plus_b2*uUP[3]*uUP[3] + P_plus_half_b2*metric->g4up[3][3] - smallb[3]*smallb[3];
}