#include "flux_source.h"
/*
 * Compute the HLLE-derived fluxes on the left face in the 0direction for all components.
 */
void calculate_HLLE_fluxes_dirn0(const primitive_quantities *restrict prims_r, const primitive_quantities *restrict prims_l, struct eos_parameters const *restrict eos, const metric_quantities *restrict metric_face, conservative_quantities *restrict cons) {

{
double cmin_dirn0, cmax_dirn0;
calculate_characteristic_speed_dirn0(prims_r, prims_l, eos, metric_face, &cmin_dirn0, &cmax_dirn0);


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
const double P_r = prims_r->press;
const double P_l = prims_l->press;
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
  const double tmp_0 = (1.0/(cmax_dirn0 + cmin_dirn0));
  const double tmp_1 = beta_faceU0*gamma_faceDD00 + beta_faceU1*gamma_faceDD01 + beta_faceU2*gamma_faceDD02;
  const double tmp_2 = beta_faceU0*gamma_faceDD01 + beta_faceU1*gamma_faceDD11 + beta_faceU2*gamma_faceDD12;
  const double tmp_3 = beta_faceU0*gamma_faceDD02 + beta_faceU1*gamma_faceDD12 + beta_faceU2*gamma_faceDD22;
  const double tmp_4 = BlU0*(gamma_faceDD00*u4lU1 + gamma_faceDD01*u4lU2 + gamma_faceDD02*u4lU3 + tmp_1*u4lU0) + BlU1*(gamma_faceDD01*u4lU1 + gamma_faceDD11*u4lU2 + gamma_faceDD12*u4lU3 + tmp_2*u4lU0) + BlU2*(gamma_faceDD02*u4lU1 + gamma_faceDD12*u4lU2 + gamma_faceDD22*u4lU3 + tmp_3*u4lU0);
  const double tmp_5 = BlU0 + tmp_4*u4lU1;
  const double tmp_7 = ((alpha_face)*(alpha_face));
  const double tmp_8 = (1.0/(tmp_7));
  const double tmp_9 = tmp_8/((SQRT_4_PI)*(SQRT_4_PI));
  const double tmp_10 = tmp_9/((u4lU0)*(u4lU0));
  const double tmp_11 = tmp_10*((tmp_5)*(tmp_5));
  const double tmp_13 = BlU1 + tmp_4*u4lU2;
  const double tmp_14 = gamma_faceDD11*tmp_10*((tmp_13)*(tmp_13));
  const double tmp_15 = BlU2 + tmp_4*u4lU3;
  const double tmp_16 = gamma_faceDD22*tmp_10*((tmp_15)*(tmp_15));
  const double tmp_17 = beta_faceU0*tmp_1 + beta_faceU1*tmp_2 + beta_faceU2*tmp_3 - tmp_7;
  const double tmp_18 = ((tmp_4)*(tmp_4))*tmp_9;
  const double tmp_22 = tmp_10*tmp_13*tmp_5;
  const double tmp_23 = tmp_10*tmp_15*tmp_5;
  const double tmp_25 = gamma_faceDD12*tmp_10*tmp_13*tmp_15;
  const double tmp_26 = tmp_4*tmp_9/u4lU0;
  const double tmp_28 = tmp_1*tmp_26*tmp_5;
  const double tmp_30 = tmp_13*tmp_2*tmp_26;
  const double tmp_32 = tmp_15*tmp_26*tmp_3;
  const double tmp_33 = gamma_faceDD00*tmp_11 + 2*gamma_faceDD01*tmp_22 + 2*gamma_faceDD02*tmp_23 + h_l*rhob_l + tmp_14 + tmp_16 + tmp_17*tmp_18 + 2*tmp_25 + 2*tmp_28 + 2*tmp_30 + 2*tmp_32;
  const double tmp_35 = gamma_faceDD00*gamma_faceDD11*gamma_faceDD22 - gamma_faceDD00*((gamma_faceDD12)*(gamma_faceDD12)) - ((gamma_faceDD01)*(gamma_faceDD01))*gamma_faceDD22 + 2*gamma_faceDD01*gamma_faceDD02*gamma_faceDD12 - ((gamma_faceDD02)*(gamma_faceDD02))*gamma_faceDD11;
  const double tmp_36 = (1.0/(tmp_35));
  const double tmp_37 = -((beta_faceU0)*(beta_faceU0))*tmp_8 + tmp_36*(gamma_faceDD11*gamma_faceDD22 - ((gamma_faceDD12)*(gamma_faceDD12)));
  const double tmp_38 = P_l + (1.0/2.0)*gamma_faceDD00*tmp_11 + gamma_faceDD01*tmp_22 + gamma_faceDD02*tmp_23 + (1.0/2.0)*tmp_14 + (1.0/2.0)*tmp_16 + (1.0/2.0)*tmp_17*tmp_18 + tmp_25 + tmp_28 + tmp_30 + tmp_32;
  const double tmp_39 = -tmp_11 + tmp_33*((u4lU1)*(u4lU1)) + tmp_37*tmp_38;
  const double tmp_40 = tmp_38*tmp_8;
  const double tmp_41 = tmp_33*u4lU0;
  const double tmp_42 = beta_faceU0*tmp_40 - tmp_26*tmp_5 + tmp_41*u4lU1;
  const double tmp_45 = -beta_faceU0*beta_faceU1*tmp_8 + tmp_36*(-gamma_faceDD01*gamma_faceDD22 + gamma_faceDD02*gamma_faceDD12);
  const double tmp_46 = -tmp_22 + tmp_33*u4lU1*u4lU2 + tmp_38*tmp_45;
  const double tmp_47 = -beta_faceU0*beta_faceU2*tmp_8 + tmp_36*(gamma_faceDD01*gamma_faceDD12 - gamma_faceDD02*gamma_faceDD11);
  const double tmp_48 = -tmp_23 + tmp_33*u4lU1*u4lU3 + tmp_38*tmp_47;
  const double tmp_49 = sqrt(tmp_35);
  const double tmp_50 = alpha_face*tmp_49;
  const double tmp_51 = cmax_dirn0*tmp_50;
  const double tmp_52 = BrU0*(gamma_faceDD00*u4rU1 + gamma_faceDD01*u4rU2 + gamma_faceDD02*u4rU3 + tmp_1*u4rU0) + BrU1*(gamma_faceDD01*u4rU1 + gamma_faceDD11*u4rU2 + gamma_faceDD12*u4rU3 + tmp_2*u4rU0) + BrU2*(gamma_faceDD02*u4rU1 + gamma_faceDD12*u4rU2 + gamma_faceDD22*u4rU3 + tmp_3*u4rU0);
  const double tmp_53 = BrU0 + tmp_52*u4rU1;
  const double tmp_55 = tmp_9/((u4rU0)*(u4rU0));
  const double tmp_56 = ((tmp_53)*(tmp_53))*tmp_55;
  const double tmp_58 = BrU1 + tmp_52*u4rU2;
  const double tmp_59 = gamma_faceDD11*tmp_55*((tmp_58)*(tmp_58));
  const double tmp_60 = BrU2 + tmp_52*u4rU3;
  const double tmp_61 = gamma_faceDD22*tmp_55*((tmp_60)*(tmp_60));
  const double tmp_62 = ((tmp_52)*(tmp_52))*tmp_9;
  const double tmp_65 = tmp_53*tmp_55*tmp_58;
  const double tmp_66 = tmp_53*tmp_55*tmp_60;
  const double tmp_68 = gamma_faceDD12*tmp_55*tmp_58*tmp_60;
  const double tmp_69 = tmp_52*tmp_9/u4rU0;
  const double tmp_71 = tmp_1*tmp_53*tmp_69;
  const double tmp_73 = tmp_2*tmp_58*tmp_69;
  const double tmp_75 = tmp_3*tmp_60*tmp_69;
  const double tmp_76 = gamma_faceDD00*tmp_56 + 2*gamma_faceDD01*tmp_65 + 2*gamma_faceDD02*tmp_66 + h_r*rhob_r + tmp_17*tmp_62 + tmp_59 + tmp_61 + 2*tmp_68 + 2*tmp_71 + 2*tmp_73 + 2*tmp_75;
  const double tmp_77 = P_r + (1.0/2.0)*gamma_faceDD00*tmp_56 + gamma_faceDD01*tmp_65 + gamma_faceDD02*tmp_66 + (1.0/2.0)*tmp_17*tmp_62 + (1.0/2.0)*tmp_59 + (1.0/2.0)*tmp_61 + tmp_68 + tmp_71 + tmp_73 + tmp_75;
  const double tmp_78 = tmp_37*tmp_77 - tmp_56 + tmp_76*((u4rU1)*(u4rU1));
  const double tmp_79 = tmp_77*tmp_8;
  const double tmp_80 = tmp_76*u4rU0;
  const double tmp_81 = beta_faceU0*tmp_79 - tmp_53*tmp_69 + tmp_80*u4rU1;
  const double tmp_83 = tmp_45*tmp_77 - tmp_65 + tmp_76*u4rU1*u4rU2;
  const double tmp_84 = tmp_47*tmp_77 - tmp_66 + tmp_76*u4rU1*u4rU3;
  const double tmp_85 = cmin_dirn0*tmp_50;
  const double tmp_86 = -tmp_18 + tmp_33*((u4lU0)*(u4lU0)) - tmp_40;
  const double tmp_87 = beta_faceU1*tmp_40 - tmp_13*tmp_26 + tmp_41*u4lU2;
  const double tmp_88 = beta_faceU2*tmp_40 - tmp_15*tmp_26 + tmp_41*u4lU3;
  const double tmp_89 = -tmp_62 + tmp_76*((u4rU0)*(u4rU0)) - tmp_79;
  const double tmp_90 = beta_faceU1*tmp_79 - tmp_58*tmp_69 + tmp_80*u4rU2;
  const double tmp_91 = beta_faceU2*tmp_79 - tmp_60*tmp_69 + tmp_80*u4rU3;
  const double tmp_92 = cmax_dirn0*cmin_dirn0;
  const double tmp_94 = rhob_l*tmp_50*u4lU0;
  const double tmp_96 = rhob_r*tmp_50*u4rU0;
  const double tmp_97 = tmp_49*tmp_7;
  cons->S_x = tmp_0*(tmp_51*(gamma_faceDD00*tmp_39 + gamma_faceDD01*tmp_46 + gamma_faceDD02*tmp_48 + tmp_1*tmp_42) + tmp_85*(gamma_faceDD00*tmp_78 + gamma_faceDD01*tmp_83 + gamma_faceDD02*tmp_84 + tmp_1*tmp_81) - tmp_92*(alpha_face*tmp_49*(gamma_faceDD00*tmp_81 + gamma_faceDD01*tmp_90 + gamma_faceDD02*tmp_91 + tmp_1*tmp_89) - tmp_50*(gamma_faceDD00*tmp_42 + gamma_faceDD01*tmp_87 + gamma_faceDD02*tmp_88 + tmp_1*tmp_86)));
  cons->S_y = tmp_0*(tmp_51*(gamma_faceDD01*tmp_39 + gamma_faceDD11*tmp_46 + gamma_faceDD12*tmp_48 + tmp_2*tmp_42) + tmp_85*(gamma_faceDD01*tmp_78 + gamma_faceDD11*tmp_83 + gamma_faceDD12*tmp_84 + tmp_2*tmp_81) - tmp_92*(alpha_face*tmp_49*(gamma_faceDD01*tmp_81 + gamma_faceDD11*tmp_90 + gamma_faceDD12*tmp_91 + tmp_2*tmp_89) - tmp_50*(gamma_faceDD01*tmp_42 + gamma_faceDD11*tmp_87 + gamma_faceDD12*tmp_88 + tmp_2*tmp_86)));
  cons->S_z = tmp_0*(tmp_51*(gamma_faceDD02*tmp_39 + gamma_faceDD12*tmp_46 + gamma_faceDD22*tmp_48 + tmp_3*tmp_42) + tmp_85*(gamma_faceDD02*tmp_78 + gamma_faceDD12*tmp_83 + gamma_faceDD22*tmp_84 + tmp_3*tmp_81) - tmp_92*(alpha_face*tmp_49*(gamma_faceDD02*tmp_81 + gamma_faceDD12*tmp_90 + gamma_faceDD22*tmp_91 + tmp_3*tmp_89) - tmp_50*(gamma_faceDD02*tmp_42 + gamma_faceDD12*tmp_87 + gamma_faceDD22*tmp_88 + tmp_3*tmp_86)));
  cons->rho = tmp_0*(rhob_l*tmp_51*u4lU1 + rhob_r*tmp_85*u4rU1 - tmp_92*(-tmp_94 + tmp_96));
  cons->tau = tmp_0*(cmax_dirn0*(-rhob_l*tmp_50*u4lU1 + tmp_42*tmp_97) + cmin_dirn0*(-rhob_r*tmp_50*u4rU1 + tmp_81*tmp_97) - tmp_92*(tmp_49*tmp_7*tmp_89 - tmp_86*tmp_97 + tmp_94 - tmp_96));
}
}
