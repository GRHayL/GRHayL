#include "ghl_flux_source.h"
/*
 * Add source terms for Stilde and tau_tilde
 */
void ghl_calculate_source_terms(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      const ghl_metric_quantities *restrict metric,
      const ghl_metric_quantities *restrict metric_d0,
      const ghl_metric_quantities *restrict metric_d1,
      const ghl_metric_quantities *restrict metric_d2,
      const ghl_extrinsic_curvature *restrict curv,
      ghl_conservative_quantities *restrict cons) {

  double h, cs2;
  ghl_compute_h_and_cs2(eos, prims, &h, &cs2);

  const double sqrt_detg = metric->lapse*metric->sqrt_detgamma;
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(metric, &metric_aux);
  const double betaD[3] = {metric_aux.g4DD[0][1],
                           metric_aux.g4DD[0][2],
                           metric_aux.g4DD[0][3]};

  // All of these vectors compute derivatives D^j(beta_i) using D^j(beta^k gamma_ik)
  const double d0_betaD[3] = {metric->betaU[0]*metric_d0->gammaDD[0][0] + metric->betaU[1]*metric_d0->gammaDD[0][1] + metric->betaU[2]*metric_d0->gammaDD[0][2]
                            + metric_d0->betaU[0]*metric->gammaDD[0][0] + metric_d0->betaU[1]*metric->gammaDD[0][1] + metric_d0->betaU[2]*metric->gammaDD[0][2],
                              metric->betaU[0]*metric_d0->gammaDD[0][1] + metric->betaU[1]*metric_d0->gammaDD[1][1] + metric->betaU[2]*metric_d0->gammaDD[1][2]
                            + metric_d0->betaU[0]*metric->gammaDD[0][1] + metric_d0->betaU[1]*metric->gammaDD[1][1] + metric_d0->betaU[2]*metric->gammaDD[1][2],
                              metric->betaU[0]*metric_d0->gammaDD[0][2] + metric->betaU[1]*metric_d0->gammaDD[1][2] + metric->betaU[2]*metric_d0->gammaDD[2][2]
                            + metric_d0->betaU[0]*metric->gammaDD[0][2] + metric_d0->betaU[1]*metric->gammaDD[1][2] + metric_d0->betaU[2]*metric->gammaDD[2][2]};

  const double d1_betaD[3] = {metric->betaU[0]*metric_d1->gammaDD[0][0] + metric->betaU[1]*metric_d1->gammaDD[0][1] + metric->betaU[2]*metric_d1->gammaDD[0][2]
                            + metric_d1->betaU[0]*metric->gammaDD[0][0] + metric_d1->betaU[1]*metric->gammaDD[0][1] + metric_d1->betaU[2]*metric->gammaDD[0][2],
                              metric->betaU[0]*metric_d1->gammaDD[0][1] + metric->betaU[1]*metric_d1->gammaDD[1][1] + metric->betaU[2]*metric_d1->gammaDD[1][2]
                            + metric_d1->betaU[0]*metric->gammaDD[0][1] + metric_d1->betaU[1]*metric->gammaDD[1][1] + metric_d1->betaU[2]*metric->gammaDD[1][2],
                              metric->betaU[0]*metric_d1->gammaDD[0][2] + metric->betaU[1]*metric_d1->gammaDD[1][2] + metric->betaU[2]*metric_d1->gammaDD[2][2]
                            + metric_d1->betaU[0]*metric->gammaDD[0][2] + metric_d1->betaU[1]*metric->gammaDD[1][2] + metric_d1->betaU[2]*metric->gammaDD[2][2]};

  const double d2_betaD[3] = {metric->betaU[0]*metric_d2->gammaDD[0][0] + metric->betaU[1]*metric_d2->gammaDD[0][1] + metric->betaU[2]*metric_d2->gammaDD[0][2]
                            + metric_d2->betaU[0]*metric->gammaDD[0][0] + metric_d2->betaU[1]*metric->gammaDD[0][1] + metric_d2->betaU[2]*metric->gammaDD[0][2],
                              metric->betaU[0]*metric_d2->gammaDD[0][1] + metric->betaU[1]*metric_d2->gammaDD[1][1] + metric->betaU[2]*metric_d2->gammaDD[1][2]
                            + metric_d2->betaU[0]*metric->gammaDD[0][1] + metric_d2->betaU[1]*metric->gammaDD[1][1] + metric_d2->betaU[2]*metric->gammaDD[1][2],
                              metric->betaU[0]*metric_d2->gammaDD[0][2] + metric->betaU[1]*metric_d2->gammaDD[1][2] + metric->betaU[2]*metric_d2->gammaDD[2][2]
                            + metric_d2->betaU[0]*metric->gammaDD[0][2] + metric_d2->betaU[1]*metric->gammaDD[1][2] + metric_d2->betaU[2]*metric->gammaDD[2][2]};

  // This contains the derivatives of g4_00: D(beta^i beta_i - alpha^2)
  const double d_g4DD00[3] = {metric->betaU[0]*d0_betaD[0] + metric->betaU[1]*d0_betaD[1] + metric->betaU[2]*d0_betaD[2]
                                + metric_d0->betaU[0]*betaD[0] + metric_d0->betaU[1]*betaD[1] + metric_d0->betaU[2]*betaD[2]
                                - 2.0*metric->lapse*metric_d0->lapse,
                              metric->betaU[0]*d1_betaD[0] + metric->betaU[1]*d1_betaD[1] + metric->betaU[2]*d1_betaD[2]
                                + metric_d1->betaU[0]*betaD[0] + metric_d1->betaU[1]*betaD[1] + metric_d1->betaU[2]*betaD[2]
                                - 2.0*metric->lapse*metric_d1->lapse,
                              metric->betaU[0]*d2_betaD[0] + metric->betaU[1]*d2_betaD[1] + metric->betaU[2]*d2_betaD[2]
                                + metric_d2->betaU[0]*betaD[0] + metric_d2->betaU[1]*betaD[1] + metric_d2->betaU[2]*betaD[2]
                                - 2.0*metric->lapse*metric_d2->lapse};

  double vD[3];
  ghl_raise_lower_vector_3D(metric->gammaDD, prims->vU, vD);
  const double Bvtilde = prims->u0*(prims->BU[0]*(vD[0] + betaD[0])
                                  + prims->BU[1]*(vD[1] + betaD[1])
                                  + prims->BU[2]*(vD[2] + betaD[2]));

  const double tmp_4vecU[4] = {metric->lapseinv*Bvtilde,
                                metric->lapseinv*(prims->BU[0]/prims->u0 + Bvtilde*prims->vU[0]),
                                metric->lapseinv*(prims->BU[1]/prims->u0 + Bvtilde*prims->vU[1]),
                                metric->lapseinv*(prims->BU[2]/prims->u0 + Bvtilde*prims->vU[2])};
  const double tmp_vec2 = ghl_compute_vec2_from_vec4D( metric_aux.g4DD, tmp_4vecU);

  const double tmp_r = (tmp_vec2 + h*prims->rho)*prims->u0*prims->u0;
  const double tmp_p = prims->press + tmp_vec2/2.0;

  const double soln_scal   =  -tmp_4vecU[0]*tmp_4vecU[0] + tmp_p*metric_aux.g4UU[0][0] + tmp_r  /* u^0 */    /* u^0 */  ;
  const double soln_vec[3] = {-tmp_4vecU[0]*tmp_4vecU[1] + tmp_p*metric_aux.g4UU[0][1] + tmp_r*prims->vU[0]  /* u^0 */  ,
                              -tmp_4vecU[0]*tmp_4vecU[2] + tmp_p*metric_aux.g4UU[0][2] + tmp_r*prims->vU[1]  /* u^0 */  ,
                              -tmp_4vecU[0]*tmp_4vecU[3] + tmp_p*metric_aux.g4UU[0][3] + tmp_r*prims->vU[2]  /* u^0 */  };
  const double soln_tsr[6] = {-tmp_4vecU[1]*tmp_4vecU[1] + tmp_p*metric_aux.g4UU[1][1] + tmp_r*prims->vU[0]*prims->vU[0],
                              -tmp_4vecU[1]*tmp_4vecU[2] + tmp_p*metric_aux.g4UU[1][2] + tmp_r*prims->vU[0]*prims->vU[1],
                              -tmp_4vecU[1]*tmp_4vecU[3] + tmp_p*metric_aux.g4UU[1][3] + tmp_r*prims->vU[0]*prims->vU[2],
                              -tmp_4vecU[2]*tmp_4vecU[2] + tmp_p*metric_aux.g4UU[2][2] + tmp_r*prims->vU[1]*prims->vU[1],
                              -tmp_4vecU[2]*tmp_4vecU[3] + tmp_p*metric_aux.g4UU[2][3] + tmp_r*prims->vU[1]*prims->vU[2],
                              -tmp_4vecU[3]*tmp_4vecU[3] + tmp_p*metric_aux.g4UU[3][3] + tmp_r*prims->vU[2]*prims->vU[2]};

  const double scal_betaU[3] = {metric->betaU[0]*soln_scal,
                                metric->betaU[1]*soln_scal,
                                metric->betaU[2]*soln_scal};

  cons->tau = sqrt_detg*(metric_d0->lapse*(-scal_betaU[0] - soln_vec[0])
                       + metric_d1->lapse*(-scal_betaU[1] - soln_vec[1])
                       + metric_d2->lapse*(-scal_betaU[2] - soln_vec[2])
                       + 2.0*curv->K[0][1]*(metric->betaU[0]*soln_vec[1] + soln_vec[0]*metric->betaU[1] + metric->betaU[0]*scal_betaU[1] + soln_tsr[1])
                       + 2.0*curv->K[0][2]*(metric->betaU[0]*soln_vec[2] + soln_vec[0]*metric->betaU[2] + metric->betaU[0]*scal_betaU[2] + soln_tsr[2])
                       + 2.0*curv->K[1][2]*(metric->betaU[1]*soln_vec[2] + soln_vec[1]*metric->betaU[2] + metric->betaU[1]*scal_betaU[2] + soln_tsr[4])
                       + curv->K[0][0]*(2.0*metric->betaU[0]*soln_vec[0] + metric->betaU[0]*scal_betaU[0] + soln_tsr[0])
                       + curv->K[1][1]*(2.0*metric->betaU[1]*soln_vec[1] + metric->betaU[1]*scal_betaU[1] + soln_tsr[3])
                       + curv->K[2][2]*(2.0*metric->betaU[2]*soln_vec[2] + metric->betaU[2]*scal_betaU[2] + soln_tsr[5]));

  cons->SD[0] = sqrt_detg*(
      soln_scal*d_g4DD00[0]
    + d0_betaD[0]*soln_vec[0] + d0_betaD[1]*soln_vec[1] + d0_betaD[2]*soln_vec[2]
    + metric_d0->gammaDD[0][1]*soln_tsr[1] + metric_d0->gammaDD[0][2]*soln_tsr[2] + metric_d0->gammaDD[1][2]*soln_tsr[4]
    + (metric_d0->gammaDD[0][0]*soln_tsr[0] + metric_d0->gammaDD[1][1]*soln_tsr[3] + metric_d0->gammaDD[2][2]*soln_tsr[5])/2.0);

  cons->SD[1] = sqrt_detg*(
      soln_scal*d_g4DD00[1]
    + d1_betaD[0]*soln_vec[0] + d1_betaD[1]*soln_vec[1] + d1_betaD[2]*soln_vec[2]
    + metric_d1->gammaDD[0][1]*soln_tsr[1] + metric_d1->gammaDD[0][2]*soln_tsr[2] + metric_d1->gammaDD[1][2]*soln_tsr[4]
    + (metric_d1->gammaDD[0][0]*soln_tsr[0] + metric_d1->gammaDD[1][1]*soln_tsr[3] + metric_d1->gammaDD[2][2]*soln_tsr[5])/2.0);

  cons->SD[2] = sqrt_detg*(
      soln_scal*d_g4DD00[2]
    + soln_vec[0]*d2_betaD[0] + soln_vec[1]*d2_betaD[1] + soln_vec[2]*d2_betaD[2]
    + metric_d2->gammaDD[0][1]*soln_tsr[1] + metric_d2->gammaDD[0][2]*soln_tsr[2] + metric_d2->gammaDD[1][2]*soln_tsr[4]
    + (metric_d2->gammaDD[0][0]*soln_tsr[0] + metric_d2->gammaDD[1][1]*soln_tsr[3] + metric_d2->gammaDD[2][2]*soln_tsr[5])/2.0);
}
