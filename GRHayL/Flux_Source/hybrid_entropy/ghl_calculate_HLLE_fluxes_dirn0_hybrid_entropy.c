#include "ghl_flux_source.h"
/*
 * Compute the HLLE-derived fluxes on the left face in direction 0 for all components.
 */
ghl_error_codes_t ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy_checked(ghl_primitive_quantities *restrict prims_r, ghl_primitive_quantities *restrict prims_l, const ghl_eos_parameters *restrict eos, const ghl_metric_quantities *restrict metric_face, const double cmin_dirn0, const double cmax_dirn0, ghl_conservative_quantities *restrict cons) {

{

const double wavespeed_scale =
      fmax(1.0, fmax(fabs(cmin_dirn0), fabs(cmax_dirn0)));
if(!isfinite(cmin_dirn0) || !isfinite(cmax_dirn0) ||
   cmin_dirn0 < -DBL_EPSILON*wavespeed_scale ||
   cmax_dirn0 < -DBL_EPSILON*wavespeed_scale)
  return ghl_error_invalid_hlle_wavespeeds;
const double cmin_clamped = fmax(cmin_dirn0, 0.0);
const double cmax_clamped = fmax(cmax_dirn0, 0.0);
if(cmin_clamped > DBL_MAX - cmax_clamped ||
   (cmin_clamped > 1.0 && cmax_clamped > DBL_MAX/cmin_clamped))
  return ghl_error_invalid_hlle_wavespeeds;
const double wavespeed_sum = cmin_clamped + cmax_clamped;
if(wavespeed_sum <= 0.0 || wavespeed_sum < 1.0/DBL_MAX)
  return ghl_error_invalid_hlle_wavespeeds;
const double wavespeed_product = cmin_clamped*cmax_clamped;
if(cmin_clamped > 0.0 && cmax_clamped > 0.0 && wavespeed_product == 0.0)
  return ghl_error_invalid_hlle_wavespeeds;
const double cmin_weight = cmin_clamped/wavespeed_sum;
const double cmax_weight = cmax_clamped/wavespeed_sum;
const double dissipation_speed =
      wavespeed_product/wavespeed_sum;

double h_r, h_l, cs2_r, cs2_l;

ghl_error_codes_t error = ghl_compute_h_and_cs2(eos, prims_r, &h_r, &cs2_r);
if(error != ghl_success)
  return error;
error = ghl_compute_h_and_cs2(eos, prims_l, &h_l, &cs2_l);
if(error != ghl_success)
  return error;
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
  const double tmp_0 = beta_faceU0*gamma_faceDD00 + beta_faceU1*gamma_faceDD01 + beta_faceU2*gamma_faceDD02;
  const double tmp_1 = beta_faceU0*gamma_faceDD01 + beta_faceU1*gamma_faceDD11 + beta_faceU2*gamma_faceDD12;
  const double tmp_2 = beta_faceU0*gamma_faceDD02 + beta_faceU1*gamma_faceDD12 + beta_faceU2*gamma_faceDD22;
  const double tmp_3 = BlU0*(gamma_faceDD00*u4lU1 + gamma_faceDD01*u4lU2 + gamma_faceDD02*u4lU3 + tmp_0*u4lU0) + BlU1*(gamma_faceDD01*u4lU1 + gamma_faceDD11*u4lU2 + gamma_faceDD12*u4lU3 + tmp_1*u4lU0) + BlU2*(gamma_faceDD02*u4lU1 + gamma_faceDD12*u4lU2 + gamma_faceDD22*u4lU3 + tmp_2*u4lU0);
  const double tmp_4 = BlU0 + tmp_3*u4lU1;
  const double tmp_6 = ((alpha_face)*(alpha_face));
  const double tmp_7 = (1.0/(tmp_6));
  const double tmp_8 = tmp_7/((SQRT_4_PI)*(SQRT_4_PI));
  const double tmp_9 = tmp_8/((u4lU0)*(u4lU0));
  const double tmp_10 = ((tmp_4)*(tmp_4))*tmp_9;
  const double tmp_11 = BlU1 + tmp_3*u4lU2;
  const double tmp_13 = tmp_11*tmp_4*tmp_9;
  const double tmp_14 = BlU2 + tmp_3*u4lU3;
  const double tmp_15 = tmp_14*tmp_4*tmp_9;
  const double tmp_16 = tmp_3*tmp_8/u4lU0;
  const double tmp_20 = gamma_faceDD01*tmp_13 + gamma_faceDD02*tmp_15 + gamma_faceDD12*tmp_11*tmp_14*tmp_9 + tmp_0*tmp_16*tmp_4 + tmp_1*tmp_11*tmp_16 + tmp_14*tmp_16*tmp_2;
  const double tmp_21 = beta_faceU0*tmp_0 + beta_faceU1*tmp_1 + beta_faceU2*tmp_2 - tmp_6;
  const double tmp_22 = ((tmp_3)*(tmp_3))*tmp_8;
  const double tmp_23 = gamma_faceDD00*tmp_10 + gamma_faceDD11*((tmp_11)*(tmp_11))*tmp_9 + gamma_faceDD22*((tmp_14)*(tmp_14))*tmp_9 + tmp_21*tmp_22;
  const double tmp_24 = _Integer_2*tmp_20 + h_l*rhob_l + tmp_23;
  const double tmp_26 = _Integer_2*gamma_faceDD01*gamma_faceDD02*gamma_faceDD12 + gamma_faceDD00*gamma_faceDD11*gamma_faceDD22 - gamma_faceDD00*((gamma_faceDD12)*(gamma_faceDD12)) - ((gamma_faceDD01)*(gamma_faceDD01))*gamma_faceDD22 - ((gamma_faceDD02)*(gamma_faceDD02))*gamma_faceDD11;
  const double tmp_27 = (1.0/(tmp_26));
  const double tmp_28 = -((beta_faceU0)*(beta_faceU0))*tmp_7 + tmp_27*(gamma_faceDD11*gamma_faceDD22 - ((gamma_faceDD12)*(gamma_faceDD12)));
  const double tmp_29 = P_l + _Rational_1_2*tmp_23 + tmp_20;
  const double tmp_30 = -tmp_10 + tmp_24*((u4lU1)*(u4lU1)) + tmp_28*tmp_29;
  const double tmp_31 = tmp_29*tmp_7;
  const double tmp_32 = tmp_24*u4lU0;
  const double tmp_33 = beta_faceU0*tmp_31 - tmp_16*tmp_4 + tmp_32*u4lU1;
  const double tmp_36 = -beta_faceU0*beta_faceU1*tmp_7 + tmp_27*(-gamma_faceDD01*gamma_faceDD22 + gamma_faceDD02*gamma_faceDD12);
  const double tmp_37 = -tmp_13 + tmp_24*u4lU1*u4lU2 + tmp_29*tmp_36;
  const double tmp_38 = -beta_faceU0*beta_faceU2*tmp_7 + tmp_27*(gamma_faceDD01*gamma_faceDD12 - gamma_faceDD02*gamma_faceDD11);
  const double tmp_39 = -tmp_15 + tmp_24*u4lU1*u4lU3 + tmp_29*tmp_38;
  const double tmp_40 = sqrt(tmp_26);
  const double tmp_41 = alpha_face*tmp_40;
  const double tmp_42 = cmax_weight*tmp_41;
  const double tmp_43 = BrU0*(gamma_faceDD00*u4rU1 + gamma_faceDD01*u4rU2 + gamma_faceDD02*u4rU3 + tmp_0*u4rU0) + BrU1*(gamma_faceDD01*u4rU1 + gamma_faceDD11*u4rU2 + gamma_faceDD12*u4rU3 + tmp_1*u4rU0) + BrU2*(gamma_faceDD02*u4rU1 + gamma_faceDD12*u4rU2 + gamma_faceDD22*u4rU3 + tmp_2*u4rU0);
  const double tmp_44 = BrU0 + tmp_43*u4rU1;
  const double tmp_46 = tmp_8/((u4rU0)*(u4rU0));
  const double tmp_47 = ((tmp_44)*(tmp_44))*tmp_46;
  const double tmp_48 = BrU1 + tmp_43*u4rU2;
  const double tmp_50 = tmp_44*tmp_46*tmp_48;
  const double tmp_51 = BrU2 + tmp_43*u4rU3;
  const double tmp_52 = tmp_44*tmp_46*tmp_51;
  const double tmp_53 = tmp_43*tmp_8/u4rU0;
  const double tmp_57 = gamma_faceDD01*tmp_50 + gamma_faceDD02*tmp_52 + gamma_faceDD12*tmp_46*tmp_48*tmp_51 + tmp_0*tmp_44*tmp_53 + tmp_1*tmp_48*tmp_53 + tmp_2*tmp_51*tmp_53;
  const double tmp_58 = ((tmp_43)*(tmp_43))*tmp_8;
  const double tmp_59 = gamma_faceDD00*tmp_47 + gamma_faceDD11*tmp_46*((tmp_48)*(tmp_48)) + gamma_faceDD22*tmp_46*((tmp_51)*(tmp_51)) + tmp_21*tmp_58;
  const double tmp_60 = _Integer_2*tmp_57 + h_r*rhob_r + tmp_59;
  const double tmp_61 = P_r + _Rational_1_2*tmp_59 + tmp_57;
  const double tmp_62 = tmp_28*tmp_61 - tmp_47 + tmp_60*((u4rU1)*(u4rU1));
  const double tmp_63 = tmp_61*tmp_7;
  const double tmp_64 = tmp_60*u4rU0;
  const double tmp_65 = beta_faceU0*tmp_63 - tmp_44*tmp_53 + tmp_64*u4rU1;
  const double tmp_67 = tmp_36*tmp_61 - tmp_50 + tmp_60*u4rU1*u4rU2;
  const double tmp_68 = tmp_38*tmp_61 - tmp_52 + tmp_60*u4rU1*u4rU3;
  const double tmp_69 = cmin_weight*tmp_41;
  const double tmp_70 = -tmp_22 + tmp_24*((u4lU0)*(u4lU0)) - tmp_31;
  const double tmp_71 = beta_faceU1*tmp_31 - tmp_11*tmp_16 + tmp_32*u4lU2;
  const double tmp_72 = beta_faceU2*tmp_31 - tmp_14*tmp_16 + tmp_32*u4lU3;
  const double tmp_73 = -tmp_58 + tmp_60*((u4rU0)*(u4rU0)) - tmp_63;
  const double tmp_74 = beta_faceU1*tmp_63 - tmp_48*tmp_53 + tmp_64*u4rU2;
  const double tmp_75 = beta_faceU2*tmp_63 - tmp_51*tmp_53 + tmp_64*u4rU3;
  const double tmp_79 = rhob_l*tmp_41*u4lU0;
  const double tmp_81 = rhob_r*tmp_41*u4rU0;
  const double tmp_82 = tmp_40*tmp_6;
  cons->SD[0] = -dissipation_speed*(alpha_face*tmp_40*(gamma_faceDD00*tmp_65 + gamma_faceDD01*tmp_74 + gamma_faceDD02*tmp_75 + tmp_0*tmp_73) - tmp_41*(gamma_faceDD00*tmp_33 + gamma_faceDD01*tmp_71 + gamma_faceDD02*tmp_72 + tmp_0*tmp_70)) + tmp_42*(gamma_faceDD00*tmp_30 + gamma_faceDD01*tmp_37 + gamma_faceDD02*tmp_39 + tmp_0*tmp_33) + tmp_69*(gamma_faceDD00*tmp_62 + gamma_faceDD01*tmp_67 + gamma_faceDD02*tmp_68 + tmp_0*tmp_65);
  cons->SD[1] = -dissipation_speed*(alpha_face*tmp_40*(gamma_faceDD01*tmp_65 + gamma_faceDD11*tmp_74 + gamma_faceDD12*tmp_75 + tmp_1*tmp_73) - tmp_41*(gamma_faceDD01*tmp_33 + gamma_faceDD11*tmp_71 + gamma_faceDD12*tmp_72 + tmp_1*tmp_70)) + tmp_42*(gamma_faceDD01*tmp_30 + gamma_faceDD11*tmp_37 + gamma_faceDD12*tmp_39 + tmp_1*tmp_33) + tmp_69*(gamma_faceDD01*tmp_62 + gamma_faceDD11*tmp_67 + gamma_faceDD12*tmp_68 + tmp_1*tmp_65);
  cons->SD[2] = -dissipation_speed*(alpha_face*tmp_40*(gamma_faceDD02*tmp_65 + gamma_faceDD12*tmp_74 + gamma_faceDD22*tmp_75 + tmp_2*tmp_73) - tmp_41*(gamma_faceDD02*tmp_33 + gamma_faceDD12*tmp_71 + gamma_faceDD22*tmp_72 + tmp_2*tmp_70)) + tmp_42*(gamma_faceDD02*tmp_30 + gamma_faceDD12*tmp_37 + gamma_faceDD22*tmp_39 + tmp_2*tmp_33) + tmp_69*(gamma_faceDD02*tmp_62 + gamma_faceDD12*tmp_67 + gamma_faceDD22*tmp_68 + tmp_2*tmp_65);
  cons->rho = -dissipation_speed*(-tmp_79 + tmp_81) + rhob_l*tmp_42*u4lU1 + rhob_r*tmp_69*u4rU1;
  cons->tau = cmax_weight*(-rhob_l*tmp_41*u4lU1 + tmp_33*tmp_82) + cmin_weight*(-rhob_r*tmp_41*u4rU1 + tmp_65*tmp_82) - dissipation_speed*(tmp_40*tmp_6*tmp_73 - tmp_70*tmp_82 + tmp_79 - tmp_81);
  cons->entropy = S_l*tmp_42*u4lU1 + S_r*tmp_69*u4rU1 - dissipation_speed*(-S_l*tmp_41*u4lU0 + S_r*alpha_face*tmp_40*u4rU0);
}
return ghl_success;
}

void ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy(ghl_primitive_quantities *restrict prims_r, ghl_primitive_quantities *restrict prims_l, const ghl_eos_parameters *restrict eos, const ghl_metric_quantities *restrict metric_face, const double cmin_dirn0, const double cmax_dirn0, ghl_conservative_quantities *restrict cons) {
  ghl_abort_if_error(ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy_checked(prims_r, prims_l, eos, metric_face, cmin_dirn0, cmax_dirn0, cons));
}
