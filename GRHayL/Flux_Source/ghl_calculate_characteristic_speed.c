#include "ghl_flux_source.h"
void ghl_calculate_characteristic_speed(
      const int direction,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      ghl_primitive_quantities *restrict prims_r,
      ghl_primitive_quantities *restrict prims_l,
      double *cmin,
      double *cmax) {

  double h_r, h_l, cs2_r, cs2_l;

  ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
  ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);

  double betaD[3];
  ghl_raise_lower_vector_3D(ADM_metric->gammaDD, ADM_metric->betaU, betaD);
  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(ADM_metric, &metric_aux);

  double vlD[3];
  ghl_raise_lower_vector_3D(ADM_metric->gammaDD, prims_l->vU, vlD);
  const double Bvtilde_l = prims_l->u0*prims_l->u0*(prims_l->BU[0]*(vlD[0] + betaD[0]) + prims_l->BU[1]*(vlD[1] + betaD[1]) + prims_l->BU[2]*(vlD[2] + betaD[2]));
  const double tmp_vec_l[4] = {Bvtilde_l, prims_l->BU[0] + Bvtilde_l*prims_l->vU[0], prims_l->BU[1] + Bvtilde_l*prims_l->vU[1], prims_l->BU[2] + Bvtilde_l*prims_l->vU[2]};
  const double tmp_vec2_l = ghl_compute_vec2_from_vec4D(metric_aux.g4DD, tmp_vec_l)/SQR(ADM_metric->lapse*prims_l->u0);
  const double tmp_1 = tmp_vec2_l/(h_l*prims_l->rho + tmp_vec2_l);
  const double tmp_2 = cs2_l*(1.0 - tmp_1) + tmp_1;
  const double tmp_3 = (1.0 - tmp_2)*prims_l->u0*prims_l->u0;
  const double tmp_4 = -metric_aux.g4UU[0][direction+1]*tmp_2 + tmp_3*prims_l->vU[direction];
  const double tmp_5 = ADM_metric->lapseinv2*tmp_2 + tmp_3;
  const double tmp_6 = (1.0/(tmp_5));

  double vrD[3];
  ghl_raise_lower_vector_3D(ADM_metric->gammaDD, prims_r->vU, vrD);
  const double Bvtilde_r = prims_r->u0*prims_r->u0*(prims_r->BU[0]*(vrD[0] + betaD[0]) + prims_r->BU[1]*(vrD[1] + betaD[1]) + prims_r->BU[2]*(vrD[2] + betaD[2]));
  const double tmp_vec_r[4] = {Bvtilde_r, prims_r->BU[0] + Bvtilde_r*prims_r->vU[0], prims_r->BU[1] + Bvtilde_r*prims_r->vU[1], prims_r->BU[2] + Bvtilde_r*prims_r->vU[2]};
  const double tmp_vec2_r = ghl_compute_vec2_from_vec4D(metric_aux.g4DD, tmp_vec_r)/SQR(ADM_metric->lapse*prims_r->u0);
  const double tmp_7 = tmp_vec2_r/(h_r*prims_r->rho + tmp_vec2_r);
  const double tmp_8 = cs2_r*(1.0 - tmp_7) + tmp_7;
  const double tmp_9 = (1.0 - tmp_8)*prims_r->u0*prims_r->u0;
  const double tmp_10 = -metric_aux.g4UU[0][direction+1]*tmp_8 + tmp_9*prims_r->vU[direction];
  const double tmp_11 = ADM_metric->lapseinv2*tmp_8 + tmp_9;
  const double tmp_12 = (1.0/(tmp_11));

  const double tmp_15 = -2.0*(tmp_4*tmp_6 + tmp_10*tmp_12);
  const double tmp_16 = tmp_4*tmp_6 - tmp_10*tmp_12;

  const double tmp_17 = tmp_4*tmp_4 - tmp_5*(-tmp_2*metric_aux.g4UU[direction+1][direction+1] + tmp_3*prims_l->vU[direction]*prims_l->vU[direction]);
  const double tmp_18 = tmp_10*tmp_10 - tmp_11*(-tmp_8*metric_aux.g4UU[direction+1][direction+1] + tmp_9*prims_r->vU[direction]*prims_r->vU[direction]);
  const double tmp_19 = fabs(tmp_6*sqrt(2.0*(tmp_17 + fabs(tmp_17))));
  const double tmp_20 = fabs(tmp_12*sqrt(2.0*(tmp_18 + fabs(tmp_18))));
  const double tmp_21 = tmp_19 + tmp_20;
  const double tmp_22 = tmp_19 - tmp_20;
  const double tmp_24 = fabs(-tmp_22 + 2.0*tmp_16);
  const double tmp_25 = fabs( tmp_22 + 2.0*tmp_16);
  *cmin = ( tmp_15 + fabs(-tmp_24 - tmp_21 - tmp_15) + tmp_24 + tmp_21)/8.0;
  *cmax = (-tmp_15 + fabs( tmp_25 + tmp_21 - tmp_15) + tmp_25 + tmp_21)/8.0;
}
