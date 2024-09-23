#include "ghl_flux_source.h"
/*
 * Compute the HLLE-derived fluxes on the left face in the 0direction for all components.
 */
void ghl_calculate_HLLE_fluxes_hybrid_entropy(
      const int direction,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      const double cmin,
      const double cmax,
      ghl_conservative_quantities *restrict cons) {

  double h_r, h_l, cs2_r, cs2_l;
  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  const double sqrt_detg = ADM_metric->lapse*ADM_metric->sqrt_detgamma;
  const double sqrt_detg_over_csum = sqrt_detg/(cmax + cmin);
  const double cmul = cmax*cmin;

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);
  const double betaD[3] = {metric_aux.g4DD[0][1],
                           metric_aux.g4DD[0][2],
                           metric_aux.g4DD[0][3]};

  double vlD[3];
  ghl_raise_lower_vector_3D(ADM_metric->gammaDD, prims_l->vU, vlD);
  const double Bvtilde_l = prims_l->u0*(prims_l->BU[0]*(vlD[0] + betaD[0])
                                      + prims_l->BU[1]*(vlD[1] + betaD[1])
                                      + prims_l->BU[2]*(vlD[2] + betaD[2]));

  const double tmp_4vec_lU[4] = {ADM_metric->lapseinv*Bvtilde_l,
                                ADM_metric->lapseinv*(prims_l->BU[0]/prims_l->u0 + Bvtilde_l*prims_l->vU[0]),
                                ADM_metric->lapseinv*(prims_l->BU[1]/prims_l->u0 + Bvtilde_l*prims_l->vU[1]),
                                ADM_metric->lapseinv*(prims_l->BU[2]/prims_l->u0 + Bvtilde_l*prims_l->vU[2])};
  const double tmp_vec2_l = ghl_compute_vec2_from_vec4D( metric_aux.g4DD, tmp_4vec_lU);

  const double tmp_1 = (tmp_vec2_l + h_l*prims_l->rho)*prims_l->u0*prims_l->u0;
  const double tmp_2 = prims_l->press + tmp_vec2_l/2.0;

  // These variables represent various elements of the solution, some of which
  // change w.r.t. the direction of the flux.
  const double soln_scal_l    =  -tmp_4vec_lU[0]*tmp_4vec_lU[0] + tmp_1   /* u^0 */      /* u^0 */    + tmp_2*metric_aux.g4UU[0][0];
  const double soln_vec1_l[3] = {-tmp_4vec_lU[0]*tmp_4vec_lU[1] + tmp_1*prims_l->vU[0]   /* u^0 */    + tmp_2*metric_aux.g4UU[0][1],
                                 -tmp_4vec_lU[0]*tmp_4vec_lU[2] + tmp_1*prims_l->vU[1]   /* u^0 */    + tmp_2*metric_aux.g4UU[0][2],
                                 -tmp_4vec_lU[0]*tmp_4vec_lU[3] + tmp_1*prims_l->vU[2]   /* u^0 */    + tmp_2*metric_aux.g4UU[0][3]};
  const double soln_vec2_l[3] = {-tmp_4vec_lU[direction+1]*tmp_4vec_lU[1] + tmp_1*prims_l->vU[direction]*prims_l->vU[0] + tmp_2*metric_aux.g4UU[direction+1][1],
                                 -tmp_4vec_lU[direction+1]*tmp_4vec_lU[2] + tmp_1*prims_l->vU[direction]*prims_l->vU[1] + tmp_2*metric_aux.g4UU[direction+1][2],
                                 -tmp_4vec_lU[direction+1]*tmp_4vec_lU[3] + tmp_1*prims_l->vU[direction]*prims_l->vU[2] + tmp_2*metric_aux.g4UU[direction+1][3]};

  double vrD[3];
  ghl_raise_lower_vector_3D(ADM_metric->gammaDD, prims_r->vU, vrD);
  const double Bvtilde_r = prims_r->u0*(prims_r->BU[0]*(vrD[0] + betaD[0])
                                      + prims_r->BU[1]*(vrD[1] + betaD[1])
                                      + prims_r->BU[2]*(vrD[2] + betaD[2]));

  const double tmp_4vec_rU[4] = {ADM_metric->lapseinv*Bvtilde_r,
                                ADM_metric->lapseinv*(prims_r->BU[0]/prims_r->u0 + Bvtilde_r*prims_r->vU[0]),
                                ADM_metric->lapseinv*(prims_r->BU[1]/prims_r->u0 + Bvtilde_r*prims_r->vU[1]),
                                ADM_metric->lapseinv*(prims_r->BU[2]/prims_r->u0 + Bvtilde_r*prims_r->vU[2])};
  const double tmp_vec2_r = ghl_compute_vec2_from_vec4D( metric_aux.g4DD, tmp_4vec_rU);

  const double tmp_3 = (tmp_vec2_r + h_r*prims_r->rho)*prims_r->u0*prims_r->u0;
  const double tmp_4 = prims_r->press + tmp_vec2_r/2.0;

  const double soln_scal_r    =  -tmp_4vec_rU[0]*tmp_4vec_rU[0] + tmp_3   /* u^0 */      /* u^0 */    + tmp_4*metric_aux.g4UU[0][0];
  const double soln_vec1_r[3] = {-tmp_4vec_rU[0]*tmp_4vec_rU[1] + tmp_3*prims_r->vU[0]   /* u^0 */    + tmp_4*metric_aux.g4UU[0][1],
                                 -tmp_4vec_rU[0]*tmp_4vec_rU[2] + tmp_3*prims_r->vU[1]   /* u^0 */    + tmp_4*metric_aux.g4UU[0][2],
                                 -tmp_4vec_rU[0]*tmp_4vec_rU[3] + tmp_3*prims_r->vU[2]   /* u^0 */    + tmp_4*metric_aux.g4UU[0][3]};
  const double soln_vec2_r[3] = {-tmp_4vec_rU[direction+1]*tmp_4vec_rU[1] + tmp_3*prims_r->vU[direction]*prims_r->vU[0] + tmp_4*metric_aux.g4UU[direction+1][1],
                                 -tmp_4vec_rU[direction+1]*tmp_4vec_rU[2] + tmp_3*prims_r->vU[direction]*prims_r->vU[1] + tmp_4*metric_aux.g4UU[direction+1][2],
                                 -tmp_4vec_rU[direction+1]*tmp_4vec_rU[3] + tmp_3*prims_r->vU[direction]*prims_r->vU[2] + tmp_4*metric_aux.g4UU[direction+1][3]};

  cons->rho = sqrt_detg_over_csum*(prims_l->rho*cmax*prims_l->vU[direction]*prims_l->u0
                                 + prims_r->rho*cmin*prims_r->vU[direction]*prims_r->u0
                                 + cmul*(prims_l->rho*prims_l->u0 - prims_r->rho*prims_r->u0));

  cons->tau = ADM_metric->lapse*sqrt_detg_over_csum*(cmax*soln_vec1_l[direction] + cmin*soln_vec1_r[direction] - cmul*(soln_scal_r - soln_scal_l)) - cons->rho;

  cons->SD[0] = sqrt_detg_over_csum*(metric_aux.g4DD[1][1]*(cmax*soln_vec2_l[0] + cmin*soln_vec2_r[0] + cmul*(soln_vec1_l[0] - soln_vec1_r[0]))
                                   + metric_aux.g4DD[1][2]*(cmax*soln_vec2_l[1] + cmin*soln_vec2_r[1] + cmul*(soln_vec1_l[1] - soln_vec1_r[1]))
                                   + metric_aux.g4DD[1][3]*(cmax*soln_vec2_l[2] + cmin*soln_vec2_r[2] + cmul*(soln_vec1_l[2] - soln_vec1_r[2]))
                                   + betaD[0]*(cmax*soln_vec1_l[direction] + cmin*soln_vec1_r[direction] + cmul*(soln_scal_l - soln_scal_r)));

  cons->SD[1] = sqrt_detg_over_csum*(metric_aux.g4DD[1][2]*(cmax*soln_vec2_l[0] + cmin*soln_vec2_r[0] + cmul*(soln_vec1_l[0] - soln_vec1_r[0]))
                                   + metric_aux.g4DD[2][2]*(cmax*soln_vec2_l[1] + cmin*soln_vec2_r[1] + cmul*(soln_vec1_l[1] - soln_vec1_r[1]))
                                   + metric_aux.g4DD[2][3]*(cmax*soln_vec2_l[2] + cmin*soln_vec2_r[2] + cmul*(soln_vec1_l[2] - soln_vec1_r[2]))
                                   + betaD[1]*(cmax*soln_vec1_l[direction] + cmin*soln_vec1_r[direction] + cmul*(soln_scal_l - soln_scal_r)));

  cons->SD[2] = sqrt_detg_over_csum*(metric_aux.g4DD[1][3]*(cmax*soln_vec2_l[0] + cmin*soln_vec2_r[0] + cmul*(soln_vec1_l[0] - soln_vec1_r[0]))
                                   + metric_aux.g4DD[2][3]*(cmax*soln_vec2_l[1] + cmin*soln_vec2_r[1] + cmul*(soln_vec1_l[1] - soln_vec1_r[1]))
                                   + metric_aux.g4DD[3][3]*(cmax*soln_vec2_l[2] + cmin*soln_vec2_r[2] + cmul*(soln_vec1_l[2] - soln_vec1_r[2]))
                                   + betaD[2]*(cmax*soln_vec1_l[direction] + cmin*soln_vec1_r[direction] + cmul*(soln_scal_l - soln_scal_r)));

  cons->entropy = sqrt_detg_over_csum*(prims_l->entropy*cmax*prims_l->vU[direction]*prims_l->u0
                                     + prims_r->entropy*cmin*prims_r->vU[direction]*prims_r->u0
                                     + cmul*(prims_l->entropy*prims_l->u0 - prims_r->entropy*prims_r->u0));
}
