#include "ghl_flux_source.h"
/*
 * Compute the HLLE-derived fluxes on the left face in the 2direction for all components.
 */
void ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy(ghl_primitive_quantities *restrict prims_r, ghl_primitive_quantities *restrict prims_l, const ghl_eos_parameters *restrict eos, const ghl_metric_quantities *restrict metric_face, const double cmin_dirn2, const double cmax_dirn2, ghl_conservative_quantities *restrict cons) {

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
const double P_r = prims_r->press;
const double P_l = prims_l->press;
const double rhob_r = prims_r->rho;
const double rhob_l = prims_l->rho;
const double S_r = prims_r->entropy;
const double S_l = prims_l->entropy;
const double Y_e_r = prims_r->Y_e;
const double Y_e_l = prims_l->Y_e;
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
  const double _Integer_2 = 2.0;
  const double _Rational_1_2 = 1.0/2.0;
  const double tmp_0 = (1.0/(cmax_dirn2 + cmin_dirn2));
  const double tmp_1 = beta_faceU0*gamma_faceDD00 + beta_faceU1*gamma_faceDD01 + beta_faceU2*gamma_faceDD02;
  const double tmp_2 = beta_faceU0*gamma_faceDD01 + beta_faceU1*gamma_faceDD11 + beta_faceU2*gamma_faceDD12;
  const double tmp_3 = beta_faceU0*gamma_faceDD02 + beta_faceU1*gamma_faceDD12 + beta_faceU2*gamma_faceDD22;
  const double tmp_4 = BlU0*(gamma_faceDD00*u4lU1 + gamma_faceDD01*u4lU2 + gamma_faceDD02*u4lU3 + tmp_1*u4lU0) + BlU1*(gamma_faceDD01*u4lU1 + gamma_faceDD11*u4lU2 + gamma_faceDD12*u4lU3 + tmp_2*u4lU0) + BlU2*(gamma_faceDD02*u4lU1 + gamma_faceDD12*u4lU2 + gamma_faceDD22*u4lU3 + tmp_3*u4lU0);
  const double tmp_5 = BlU2 + tmp_4*u4lU3;
  const double tmp_7 = ((alpha_face)*(alpha_face));
  const double tmp_8 = (1.0/(tmp_7));
  const double tmp_9 = tmp_8/((SQRT_4_PI)*(SQRT_4_PI));
  const double tmp_10 = tmp_9/((u4lU0)*(u4lU0));
  const double tmp_11 = tmp_10*((tmp_5)*(tmp_5));
  const double tmp_12 = BlU0 + tmp_4*u4lU1;
  const double tmp_13 = BlU1 + tmp_4*u4lU2;
  const double tmp_15 = tmp_10*tmp_12*tmp_5;
  const double tmp_16 = tmp_10*tmp_13*tmp_5;
  const double tmp_17 = tmp_4*tmp_9/u4lU0;
  const double tmp_21 = gamma_faceDD01*tmp_10*tmp_12*tmp_13 + gamma_faceDD02*tmp_15 + gamma_faceDD12*tmp_16 + tmp_1*tmp_12*tmp_17 + tmp_13*tmp_17*tmp_2 + tmp_17*tmp_3*tmp_5;
  const double tmp_22 = beta_faceU0*tmp_1 + beta_faceU1*tmp_2 + beta_faceU2*tmp_3 - tmp_7;
  const double tmp_23 = ((tmp_4)*(tmp_4))*tmp_9;
  const double tmp_24 = gamma_faceDD00*tmp_10*((tmp_12)*(tmp_12)) + gamma_faceDD11*tmp_10*((tmp_13)*(tmp_13)) + gamma_faceDD22*tmp_11 + tmp_22*tmp_23;
  const double tmp_25 = _Integer_2*tmp_21 + h_l*rhob_l + tmp_24;
  const double tmp_27 = _Integer_2*gamma_faceDD01*gamma_faceDD02*gamma_faceDD12 + gamma_faceDD00*gamma_faceDD11*gamma_faceDD22 - gamma_faceDD00*((gamma_faceDD12)*(gamma_faceDD12)) - ((gamma_faceDD01)*(gamma_faceDD01))*gamma_faceDD22 - ((gamma_faceDD02)*(gamma_faceDD02))*gamma_faceDD11;
  const double tmp_28 = (1.0/(tmp_27));
  const double tmp_29 = -((beta_faceU2)*(beta_faceU2))*tmp_8 + tmp_28*(gamma_faceDD00*gamma_faceDD11 - ((gamma_faceDD01)*(gamma_faceDD01)));
  const double tmp_30 = P_l + _Rational_1_2*tmp_24 + tmp_21;
  const double tmp_31 = -tmp_11 + tmp_25*((u4lU3)*(u4lU3)) + tmp_29*tmp_30;
  const double tmp_32 = tmp_30*tmp_8;
  const double tmp_33 = tmp_25*u4lU0;
  const double tmp_34 = beta_faceU2*tmp_32 - tmp_17*tmp_5 + tmp_33*u4lU3;
  const double tmp_37 = -beta_faceU0*beta_faceU2*tmp_8 + tmp_28*(gamma_faceDD01*gamma_faceDD12 - gamma_faceDD02*gamma_faceDD11);
  const double tmp_38 = -tmp_15 + tmp_25*u4lU1*u4lU3 + tmp_30*tmp_37;
  const double tmp_39 = -beta_faceU1*beta_faceU2*tmp_8 + tmp_28*(-gamma_faceDD00*gamma_faceDD12 + gamma_faceDD01*gamma_faceDD02);
  const double tmp_40 = -tmp_16 + tmp_25*u4lU2*u4lU3 + tmp_30*tmp_39;
  const double tmp_41 = sqrt(tmp_27);
  const double tmp_42 = alpha_face*tmp_41;
  const double tmp_43 = cmax_dirn2*tmp_42;
  const double tmp_44 = BrU0*(gamma_faceDD00*u4rU1 + gamma_faceDD01*u4rU2 + gamma_faceDD02*u4rU3 + tmp_1*u4rU0) + BrU1*(gamma_faceDD01*u4rU1 + gamma_faceDD11*u4rU2 + gamma_faceDD12*u4rU3 + tmp_2*u4rU0) + BrU2*(gamma_faceDD02*u4rU1 + gamma_faceDD12*u4rU2 + gamma_faceDD22*u4rU3 + tmp_3*u4rU0);
  const double tmp_45 = BrU2 + tmp_44*u4rU3;
  const double tmp_47 = tmp_9/((u4rU0)*(u4rU0));
  const double tmp_48 = ((tmp_45)*(tmp_45))*tmp_47;
  const double tmp_49 = BrU0 + tmp_44*u4rU1;
  const double tmp_50 = BrU1 + tmp_44*u4rU2;
  const double tmp_52 = tmp_45*tmp_47*tmp_49;
  const double tmp_53 = tmp_45*tmp_47*tmp_50;
  const double tmp_54 = tmp_44*tmp_9/u4rU0;
  const double tmp_58 = gamma_faceDD01*tmp_47*tmp_49*tmp_50 + gamma_faceDD02*tmp_52 + gamma_faceDD12*tmp_53 + tmp_1*tmp_49*tmp_54 + tmp_2*tmp_50*tmp_54 + tmp_3*tmp_45*tmp_54;
  const double tmp_59 = ((tmp_44)*(tmp_44))*tmp_9;
  const double tmp_60 = gamma_faceDD00*tmp_47*((tmp_49)*(tmp_49)) + gamma_faceDD11*tmp_47*((tmp_50)*(tmp_50)) + gamma_faceDD22*tmp_48 + tmp_22*tmp_59;
  const double tmp_61 = _Integer_2*tmp_58 + h_r*rhob_r + tmp_60;
  const double tmp_62 = P_r + _Rational_1_2*tmp_60 + tmp_58;
  const double tmp_63 = tmp_29*tmp_62 - tmp_48 + tmp_61*((u4rU3)*(u4rU3));
  const double tmp_64 = tmp_62*tmp_8;
  const double tmp_65 = tmp_61*u4rU0;
  const double tmp_66 = beta_faceU2*tmp_64 - tmp_45*tmp_54 + tmp_65*u4rU3;
  const double tmp_68 = tmp_37*tmp_62 - tmp_52 + tmp_61*u4rU1*u4rU3;
  const double tmp_69 = tmp_39*tmp_62 - tmp_53 + tmp_61*u4rU2*u4rU3;
  const double tmp_70 = cmin_dirn2*tmp_42;
  const double tmp_71 = -tmp_23 + tmp_25*((u4lU0)*(u4lU0)) - tmp_32;
  const double tmp_72 = beta_faceU0*tmp_32 - tmp_12*tmp_17 + tmp_33*u4lU1;
  const double tmp_73 = beta_faceU1*tmp_32 - tmp_13*tmp_17 + tmp_33*u4lU2;
  const double tmp_74 = -tmp_59 + tmp_61*((u4rU0)*(u4rU0)) - tmp_64;
  const double tmp_75 = beta_faceU0*tmp_64 - tmp_49*tmp_54 + tmp_65*u4rU1;
  const double tmp_76 = beta_faceU1*tmp_64 - tmp_50*tmp_54 + tmp_65*u4rU2;
  const double tmp_77 = cmax_dirn2*cmin_dirn2;
  const double tmp_79 = rhob_l*tmp_43*u4lU3;
  const double tmp_81 = rhob_r*tmp_70*u4rU3;
  const double tmp_83 = rhob_l*tmp_42*u4lU0;
  const double tmp_85 = rhob_r*tmp_42*u4rU0;
  const double tmp_86 = tmp_41*tmp_7;
  cons->SD[0] = tmp_0*(tmp_43*(gamma_faceDD00*tmp_38 + gamma_faceDD01*tmp_40 + gamma_faceDD02*tmp_31 + tmp_1*tmp_34) + tmp_70*(gamma_faceDD00*tmp_68 + gamma_faceDD01*tmp_69 + gamma_faceDD02*tmp_63 + tmp_1*tmp_66) - tmp_77*(alpha_face*tmp_41*(gamma_faceDD00*tmp_75 + gamma_faceDD01*tmp_76 + gamma_faceDD02*tmp_66 + tmp_1*tmp_74) - tmp_42*(gamma_faceDD00*tmp_72 + gamma_faceDD01*tmp_73 + gamma_faceDD02*tmp_34 + tmp_1*tmp_71)));
  cons->SD[1] = tmp_0*(tmp_43*(gamma_faceDD01*tmp_38 + gamma_faceDD11*tmp_40 + gamma_faceDD12*tmp_31 + tmp_2*tmp_34) + tmp_70*(gamma_faceDD01*tmp_68 + gamma_faceDD11*tmp_69 + gamma_faceDD12*tmp_63 + tmp_2*tmp_66) - tmp_77*(alpha_face*tmp_41*(gamma_faceDD01*tmp_75 + gamma_faceDD11*tmp_76 + gamma_faceDD12*tmp_66 + tmp_2*tmp_74) - tmp_42*(gamma_faceDD01*tmp_72 + gamma_faceDD11*tmp_73 + gamma_faceDD12*tmp_34 + tmp_2*tmp_71)));
  cons->SD[2] = tmp_0*(tmp_43*(gamma_faceDD02*tmp_38 + gamma_faceDD12*tmp_40 + gamma_faceDD22*tmp_31 + tmp_3*tmp_34) + tmp_70*(gamma_faceDD02*tmp_68 + gamma_faceDD12*tmp_69 + gamma_faceDD22*tmp_63 + tmp_3*tmp_66) - tmp_77*(alpha_face*tmp_41*(gamma_faceDD02*tmp_75 + gamma_faceDD12*tmp_76 + gamma_faceDD22*tmp_66 + tmp_3*tmp_74) - tmp_42*(gamma_faceDD02*tmp_72 + gamma_faceDD12*tmp_73 + gamma_faceDD22*tmp_34 + tmp_3*tmp_71)));
  cons->rho = tmp_0*(-tmp_77*(-tmp_83 + tmp_85) + tmp_79 + tmp_81);
  cons->tau = tmp_0*(cmax_dirn2*(-rhob_l*tmp_42*u4lU3 + tmp_34*tmp_86) + cmin_dirn2*(-rhob_r*tmp_42*u4rU3 + tmp_66*tmp_86) - tmp_77*(tmp_41*tmp_7*tmp_74 - tmp_71*tmp_86 + tmp_83 - tmp_85));
  cons->entropy = tmp_0*(S_l*tmp_43*u4lU3 + S_r*tmp_70*u4rU3 - tmp_77*(-S_l*tmp_42*u4lU0 + S_r*alpha_face*tmp_41*u4rU0));
  cons->Y_e = tmp_0*(Y_e_l*tmp_79 + Y_e_r*tmp_81 - tmp_77*(-Y_e_l*tmp_83 + Y_e_r*alpha_face*rhob_r*tmp_41*u4rU0));
}
}
