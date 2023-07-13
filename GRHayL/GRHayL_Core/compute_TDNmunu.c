#include "ghl.h"

void ghl_compute_TDNmunu(
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_primitive_quantities *restrict prims,
      ghl_stress_energy *restrict Tmunu) {

  //   First set h, the enthalpy:
  const double h_enthalpy = 1.0 + prims->eps + prims->press/prims->rho;

  // Compute u^i. u^0 is provided to the function.
  const double uU[4] = {prims->u0,
                        prims->u0*prims->vU[0],
                        prims->u0*prims->vU[1],
                        prims->u0*prims->vU[2]};

  // Compute u_\alpha
  double uD[4];
  ghl_lower_vector_4D(metric_aux->g4DD, uU, uD);

  /***************************************************************/
  //                     COMPUTE TDNMUNU                         //
  /***************************************************************/
  // Compute b^{\mu} and b^2
  double smallb[4], smallb2;
  ghl_compute_smallb_and_b2(ADM_metric, prims, uD, smallb, &smallb2);

  // Precompute some useful quantities, for later:
  const double rho0_h_plus_b2 = (prims->rho*h_enthalpy + smallb2);
  const double P_plus_half_b2 = (prims->press+0.5*smallb2);

  double smallb_lower[4];
  ghl_lower_vector_4D(metric_aux->g4DD, smallb, smallb_lower);

  // Next compute T^{\mu \nu}:
  // (Eq. 33 in http://arxiv.org/pdf/astro-ph/0503420.pdf):
  // T^{mn} = (rho_0 h + b^2) u^m u^n + (P + 0.5 b^2) g^{mn} - b^m b^n, where m and n both run from 0 to 3.
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
