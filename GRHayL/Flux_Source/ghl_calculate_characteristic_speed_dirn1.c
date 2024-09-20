#include "ghl_flux_source.h"
/*
 * Compute the characteristic speeds in direction 1
 */
void ghl_calculate_characteristic_speed_dirn1(ghl_primitive_quantities *restrict prims_r, ghl_primitive_quantities *restrict prims_l, const ghl_eos_parameters *restrict eos, const ghl_metric_quantities *restrict metric_face, double *cmin_dirn1, double *cmax_dirn1) {

{

double h_r, h_l, cs2_r, cs2_l;

ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);
const double u4rU0 = prims_r->u0;
const double u4lU0 = prims_l->u0;
const double u4rU1 = prims_r->vU[0]*u4rU0;
const double u4lU1 = prims_l->vU[0]*u4lU0;
const double u4rU2 = prims_r->vU[1]*u4rU0;
const double u4lU2 = prims_l->vU[1]*u4lU0;
const double u4rU3 = prims_r->vU[2]*u4rU0;
const double u4lU3 = prims_l->vU[2]*u4lU0;
const double BrU0 = prims_r->BU[0];
const double BlU0 = prims_l->BU[0];
const double BrU1 = prims_r->BU[1];
const double BlU1 = prims_l->BU[1];
const double BrU2 = prims_r->BU[2];
const double BlU2 = prims_l->BU[2];
const double rhob_r = prims_r->rho;
const double rhob_l = prims_l->rho;
const double alpha_face = metric_face->lapse;
const double beta_faceU0 = metric_face->betaU[0];
const double beta_faceU1 = metric_face->betaU[1];
const double beta_faceU2 = metric_face->betaU[2];
const double gamma_faceDD00 = metric_face->gammaDD[0][0];
const double gamma_faceDD01 = metric_face->gammaDD[0][1];
const double gamma_faceDD02 = metric_face->gammaDD[0][2];
const double gamma_faceDD11 = metric_face->gammaDD[1][1];
const double gamma_faceDD12 = metric_face->gammaDD[1][2];
const double gamma_faceDD22 = metric_face->gammaDD[2][2];
  const double _Rational_1_2 = 1.0/2.0;
  const double _Integer_1 = 1.0;
  const double _Integer_2 = 2.0;
  const double _Integer_4 = 4.0;
  const double _Rational_1_4 = 1.0/4.0;
  const double _Rational_1_8 = 1.0/8.0;
  const double _Rational_1_16 = 1.0/16.0;
  const double tmp_1 = (1.0/((alpha_face)*(alpha_face)));
  const double tmp_2 = beta_faceU0*gamma_faceDD00 + beta_faceU1*gamma_faceDD01 + beta_faceU2*gamma_faceDD02;
  const double tmp_3 = beta_faceU0*gamma_faceDD01 + beta_faceU1*gamma_faceDD11 + beta_faceU2*gamma_faceDD12;
  const double tmp_4 = beta_faceU0*gamma_faceDD02 + beta_faceU1*gamma_faceDD12 + beta_faceU2*gamma_faceDD22;
  const double tmp_5 = BlU0*(gamma_faceDD00*u4lU1 + gamma_faceDD01*u4lU2 + gamma_faceDD02*u4lU3 + tmp_2*u4lU0) + BlU1*(gamma_faceDD01*u4lU1 + gamma_faceDD11*u4lU2 + gamma_faceDD12*u4lU3 + tmp_3*u4lU0) + BlU2*(gamma_faceDD02*u4lU1 + gamma_faceDD12*u4lU2 + gamma_faceDD22*u4lU3 + tmp_4*u4lU0);
  const double tmp_6 = BlU0 + tmp_5*u4lU1;
  const double tmp_8 = tmp_1/((SQRT_4_PI)*(SQRT_4_PI));
  const double tmp_9 = tmp_8/((u4lU0)*(u4lU0));
  const double tmp_10 = BlU1 + tmp_5*u4lU2;
  const double tmp_11 = BlU2 + tmp_5*u4lU3;
  const double tmp_12 = tmp_8*(-((alpha_face)*(alpha_face)) + beta_faceU0*tmp_2 + beta_faceU1*tmp_3 + beta_faceU2*tmp_4);
  const double tmp_15 = tmp_5/u4lU0;
  const double tmp_17 = _Integer_2*(gamma_faceDD01*tmp_10*tmp_6*tmp_9 + gamma_faceDD02*tmp_11*tmp_6*tmp_9 + gamma_faceDD12*tmp_10*tmp_11*tmp_9 + tmp_10*tmp_15*tmp_3*tmp_8 + tmp_11*tmp_15*tmp_4*tmp_8 + tmp_15*tmp_2*tmp_6*tmp_8) + gamma_faceDD00*((tmp_6)*(tmp_6))*tmp_9 + gamma_faceDD11*((tmp_10)*(tmp_10))*tmp_9 + gamma_faceDD22*((tmp_11)*(tmp_11))*tmp_9 + tmp_12*((tmp_5)*(tmp_5));
  const double tmp_18 = tmp_17/(h_l*rhob_l + tmp_17);
  const double tmp_19 = cs2_l*(_Integer_1 - tmp_18) + tmp_18;
  const double tmp_21 = _Integer_1 - tmp_19;
  const double tmp_22 = -beta_faceU1*tmp_1*tmp_19 + tmp_21*u4lU0*u4lU2;
  const double tmp_24 = tmp_1*tmp_19 + tmp_21*((u4lU0)*(u4lU0));
  const double tmp_25 = (1.0/(tmp_24));
  const double tmp_26 = -tmp_22*tmp_25;
  const double tmp_27 = BrU0*(gamma_faceDD00*u4rU1 + gamma_faceDD01*u4rU2 + gamma_faceDD02*u4rU3 + tmp_2*u4rU0) + BrU1*(gamma_faceDD01*u4rU1 + gamma_faceDD11*u4rU2 + gamma_faceDD12*u4rU3 + tmp_3*u4rU0) + BrU2*(gamma_faceDD02*u4rU1 + gamma_faceDD12*u4rU2 + gamma_faceDD22*u4rU3 + tmp_4*u4rU0);
  const double tmp_28 = BrU0 + tmp_27*u4rU1;
  const double tmp_30 = tmp_8/((u4rU0)*(u4rU0));
  const double tmp_31 = BrU1 + tmp_27*u4rU2;
  const double tmp_32 = BrU2 + tmp_27*u4rU3;
  const double tmp_34 = tmp_27/u4rU0;
  const double tmp_36 = _Integer_2*(gamma_faceDD01*tmp_28*tmp_30*tmp_31 + gamma_faceDD02*tmp_28*tmp_30*tmp_32 + gamma_faceDD12*tmp_30*tmp_31*tmp_32 + tmp_2*tmp_28*tmp_34*tmp_8 + tmp_3*tmp_31*tmp_34*tmp_8 + tmp_32*tmp_34*tmp_4*tmp_8) + gamma_faceDD00*((tmp_28)*(tmp_28))*tmp_30 + gamma_faceDD11*tmp_30*((tmp_31)*(tmp_31)) + gamma_faceDD22*tmp_30*((tmp_32)*(tmp_32)) + tmp_12*((tmp_27)*(tmp_27));
  const double tmp_37 = tmp_36/(h_r*rhob_r + tmp_36);
  const double tmp_38 = cs2_r*(_Integer_1 - tmp_37) + tmp_37;
  const double tmp_39 = _Integer_1 - tmp_38;
  const double tmp_40 = -beta_faceU1*tmp_1*tmp_38 + tmp_39*u4rU0*u4rU2;
  const double tmp_42 = tmp_1*tmp_38 + tmp_39*((u4rU0)*(u4rU0));
  const double tmp_43 = (1.0/(tmp_42));
  const double tmp_44 = -tmp_40*tmp_43;
  const double tmp_46 = -tmp_40*tmp_43 + tmp_44;
  const double tmp_47 = -tmp_22*tmp_25 + tmp_26 + tmp_46;
  const double tmp_52 = -((beta_faceU1)*(beta_faceU1))*tmp_1 + (gamma_faceDD00*gamma_faceDD22 - ((gamma_faceDD02)*(gamma_faceDD02)))/(_Integer_2*gamma_faceDD01*gamma_faceDD02*gamma_faceDD12 + gamma_faceDD00*gamma_faceDD11*gamma_faceDD22 - gamma_faceDD00*((gamma_faceDD12)*(gamma_faceDD12)) - ((gamma_faceDD01)*(gamma_faceDD01))*gamma_faceDD22 - ((gamma_faceDD02)*(gamma_faceDD02))*gamma_faceDD11);
  const double tmp_53 = ((_Integer_2)*(_Integer_2))*((tmp_22)*(tmp_22)) - _Integer_4*tmp_24*(-tmp_19*tmp_52 + tmp_21*((u4lU2)*(u4lU2)));
  const double tmp_54 = fabs(_Integer_2*_Rational_1_2*(tmp_22*tmp_25 + tmp_26) + tmp_25*sqrt(_Rational_1_2*(tmp_53 + fabs(tmp_53))));
  const double tmp_55 = ((_Integer_2)*(_Integer_2))*((tmp_40)*(tmp_40)) - _Integer_4*tmp_42*(-tmp_38*tmp_52 + tmp_39*((u4rU2)*(u4rU2)));
  const double tmp_56 = fabs(_Integer_2*_Rational_1_2*(tmp_40*tmp_43 + tmp_44) + tmp_43*sqrt(_Rational_1_2*(tmp_55 + fabs(tmp_55))));
  const double tmp_57 = tmp_54 + tmp_56;
  const double tmp_59 = _Integer_2*_Rational_1_4*(tmp_22*tmp_25 - tmp_26 + tmp_46);
  const double tmp_60 = tmp_54 - tmp_56;
  const double tmp_61 = fabs(-_Rational_1_2*tmp_60 + tmp_59);
  const double tmp_62 = -_Integer_2*tmp_47;
  const double tmp_64 = fabs(_Rational_1_2*tmp_60 + tmp_59);
  *cmin_dirn1 = _Integer_2*_Rational_1_16*tmp_47 + _Rational_1_2*fabs(-_Rational_1_2*tmp_61 - _Rational_1_4*tmp_57 + _Rational_1_8*tmp_62) + _Rational_1_4*tmp_61 + _Rational_1_8*tmp_57;
  *cmax_dirn1 = _Rational_1_16*tmp_62 + _Rational_1_2*fabs(_Rational_1_2*tmp_64 + _Rational_1_4*tmp_57 + _Rational_1_8*tmp_62) + _Rational_1_4*tmp_64 + _Rational_1_8*tmp_57;
}
}
