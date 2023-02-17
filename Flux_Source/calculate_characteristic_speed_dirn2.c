#include "flux_source.h"
/*
 * Compute the characteristic speeds in direction 2
 */
void calculate_characteristic_speed_dirn2(const primitive_quantities *restrict prims_r, const primitive_quantities *restrict prims_l, struct eos_parameters const *restrict eos, const metric_quantities *restrict metric_face, double *cmin_dirn2, double *cmax_dirn2) {

{

double h_r, h_l, cs2_r, cs2_l;

eos->compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
eos->compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);
const double u4rU0 = prims_r->ut;
const double u4lU0 = prims_l->ut;
const double u4rU1 = prims_r->vx*u4rU0;
const double u4lU1 = prims_l->vx*u4lU0;
const double u4rU2 = prims_r->vy*u4rU0;
const double u4lU2 = prims_l->vy*u4lU0;
const double u4rU3 = prims_r->vz*u4rU0;
const double u4lU3 = prims_l->vz*u4lU0;
const double BrU0 = prims_r->Bx;
const double BlU0 = prims_l->Bx;
const double BrU1 = prims_r->By;
const double BlU1 = prims_l->By;
const double BrU2 = prims_r->Bz;
const double BlU2 = prims_l->Bz;
const double rhob_r = prims_r->rho;
const double rhob_l = prims_l->rho;
const double alpha_face = metric_face->lapse;
const double beta_faceU0 = metric_face->betax;
const double beta_faceU1 = metric_face->betay;
const double beta_faceU2 = metric_face->betaz;
const double gamma_faceDD00 = metric_face->adm_gxx;
const double gamma_faceDD01 = metric_face->adm_gxy;
const double gamma_faceDD02 = metric_face->adm_gxz;
const double gamma_faceDD11 = metric_face->adm_gyy;
const double gamma_faceDD12 = metric_face->adm_gyz;
const double gamma_faceDD22 = metric_face->adm_gzz;
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
  const double tmp_16 = 2*tmp_8;
  const double tmp_18 = tmp_5/u4lU0;
  const double tmp_20 = gamma_faceDD00*((tmp_6)*(tmp_6))*tmp_9 + 2*gamma_faceDD01*tmp_10*tmp_6*tmp_9 + 2*gamma_faceDD02*tmp_11*tmp_6*tmp_9 + gamma_faceDD11*((tmp_10)*(tmp_10))*tmp_9 + 2*gamma_faceDD12*tmp_10*tmp_11*tmp_9 + gamma_faceDD22*((tmp_11)*(tmp_11))*tmp_9 + tmp_10*tmp_16*tmp_18*tmp_3 + tmp_11*tmp_16*tmp_18*tmp_4 + tmp_12*((tmp_5)*(tmp_5)) + tmp_16*tmp_18*tmp_2*tmp_6;
  const double tmp_21 = tmp_20/(h_l*rhob_l + tmp_20);
  const double tmp_23 = -cs2_l*(tmp_21 - 1);
  const double tmp_25 = tmp_1*(tmp_21 + tmp_23);
  const double tmp_26 = -tmp_21 - tmp_23 + 1;
  const double tmp_27 = tmp_26*((u4lU0)*(u4lU0));
  const double tmp_28 = (1.0/(tmp_25 + tmp_27));
  const double tmp_29 = -beta_faceU2*tmp_1*(2*tmp_21 + 2*tmp_23) + 2*tmp_26*u4lU0*u4lU3;
  const double tmp_30 = BrU0*(gamma_faceDD00*u4rU1 + gamma_faceDD01*u4rU2 + gamma_faceDD02*u4rU3 + tmp_2*u4rU0) + BrU1*(gamma_faceDD01*u4rU1 + gamma_faceDD11*u4rU2 + gamma_faceDD12*u4rU3 + tmp_3*u4rU0) + BrU2*(gamma_faceDD02*u4rU1 + gamma_faceDD12*u4rU2 + gamma_faceDD22*u4rU3 + tmp_4*u4rU0);
  const double tmp_31 = BrU0 + tmp_30*u4rU1;
  const double tmp_33 = tmp_8/((u4rU0)*(u4rU0));
  const double tmp_34 = BrU1 + tmp_30*u4rU2;
  const double tmp_35 = BrU2 + tmp_30*u4rU3;
  const double tmp_37 = 2*tmp_35;
  const double tmp_38 = tmp_30/u4rU0;
  const double tmp_39 = gamma_faceDD00*((tmp_31)*(tmp_31))*tmp_33 + 2*gamma_faceDD01*tmp_31*tmp_33*tmp_34 + gamma_faceDD02*tmp_31*tmp_33*tmp_37 + gamma_faceDD11*tmp_33*((tmp_34)*(tmp_34)) + gamma_faceDD12*tmp_33*tmp_34*tmp_37 + gamma_faceDD22*tmp_33*((tmp_35)*(tmp_35)) + tmp_12*((tmp_30)*(tmp_30)) + tmp_16*tmp_2*tmp_31*tmp_38 + tmp_16*tmp_3*tmp_34*tmp_38 + tmp_37*tmp_38*tmp_4*tmp_8;
  const double tmp_40 = tmp_39/(h_r*rhob_r + tmp_39);
  const double tmp_42 = -cs2_r*(tmp_40 - 1);
  const double tmp_44 = tmp_1*(tmp_40 + tmp_42);
  const double tmp_45 = -tmp_40 - tmp_42 + 1;
  const double tmp_46 = tmp_45*((u4rU0)*(u4rU0));
  const double tmp_47 = (1.0/(tmp_44 + tmp_46));
  const double tmp_48 = -beta_faceU2*tmp_1*(2*tmp_40 + 2*tmp_42) + 2*tmp_45*u4rU0*u4rU3;
  const double tmp_52 = -1.0/4.0*tmp_28*tmp_29;
  const double tmp_54 = (1.0/4.0)*tmp_47*tmp_48;
  const double tmp_55 = (1.0/2.0)*tmp_28;
  const double tmp_58 = -((beta_faceU2)*(beta_faceU2))*tmp_1 + (gamma_faceDD00*gamma_faceDD11 - ((gamma_faceDD01)*(gamma_faceDD01)))/(gamma_faceDD00*gamma_faceDD11*gamma_faceDD22 - gamma_faceDD00*((gamma_faceDD12)*(gamma_faceDD12)) - ((gamma_faceDD01)*(gamma_faceDD01))*gamma_faceDD22 + 2*gamma_faceDD01*gamma_faceDD02*gamma_faceDD12 - ((gamma_faceDD02)*(gamma_faceDD02))*gamma_faceDD11);
  const double tmp_59 = (4*tmp_25 + 4*tmp_27)*(tmp_26*((u4lU3)*(u4lU3)) - tmp_58*(tmp_21 + tmp_23));
  const double tmp_60 = fabs(tmp_28*sqrt((1.0/2.0)*((tmp_29)*(tmp_29)) - 1.0/2.0*tmp_59 + (1.0/2.0)*fabs(((tmp_29)*(tmp_29)) - tmp_59)));
  const double tmp_61 = (1.0/2.0)*tmp_47;
  const double tmp_63 = (4*tmp_44 + 4*tmp_46)*(tmp_45*((u4rU3)*(u4rU3)) - tmp_58*(tmp_40 + tmp_42));
  const double tmp_64 = fabs(tmp_47*sqrt((1.0/2.0)*((tmp_48)*(tmp_48)) - 1.0/2.0*tmp_63 + (1.0/2.0)*fabs(((tmp_48)*(tmp_48)) - tmp_63)));
  const double tmp_65 = (1.0/2.0)*tmp_60 - 1.0/2.0*tmp_64;
  const double tmp_66 = fabs((1.0/4.0)*tmp_28*tmp_29 - 1.0/4.0*tmp_47*tmp_48 - tmp_52 - tmp_54 - tmp_65);
  const double tmp_68 = -1.0/8.0*tmp_28*tmp_29;
  const double tmp_70 = -1.0/8.0*tmp_47*tmp_48;
  const double tmp_71 = (1.0/4.0)*tmp_60 + (1.0/4.0)*tmp_64;
  const double tmp_73 = (1.0/16.0)*tmp_28*tmp_29;
  const double tmp_75 = (1.0/16.0)*tmp_47*tmp_48;
  const double tmp_76 = -1.0/16.0*tmp_28*tmp_29;
  const double tmp_77 = -1.0/16.0*tmp_47*tmp_48;
  const double tmp_78 = (1.0/8.0)*tmp_60 + (1.0/8.0)*tmp_64;
  const double tmp_79 = fabs((1.0/4.0)*tmp_28*tmp_29 - 1.0/4.0*tmp_47*tmp_48 - tmp_52 - tmp_54 + tmp_65);
  *cmin_dirn2 = (1.0/4.0)*tmp_66 - tmp_73 - tmp_75 + tmp_76 + tmp_77 + tmp_78 + (1.0/2.0)*fabs((1.0/8.0)*tmp_28*tmp_29 + (1.0/8.0)*tmp_47*tmp_48 - 1.0/2.0*tmp_66 - tmp_68 - tmp_70 - tmp_71);
  *cmax_dirn2 = tmp_73 + tmp_75 - tmp_76 - tmp_77 + tmp_78 + (1.0/4.0)*tmp_79 + (1.0/2.0)*fabs((1.0/8.0)*tmp_28*tmp_29 + (1.0/8.0)*tmp_47*tmp_48 - tmp_68 - tmp_70 + tmp_71 + (1.0/2.0)*tmp_79);
}
}
