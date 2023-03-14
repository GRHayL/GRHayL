#include "flux_source.h"
/*
 * Add source term for tau_tilde
 */
void calculate_tau_tilde_source_term(const primitive_quantities *restrict prims, struct eos_parameters const *restrict eos, const metric_quantities *restrict metric, const extrinsic_curvature *restrict curv, const metric_derivatives *restrict metric_derivs, conservative_quantities *restrict cons) {

{

double h, cs2;

eos->compute_h_and_cs2(eos, prims, &h, &cs2);
const double u4U0 = prims->u0;
const double u4U1 = prims->vx*u4U0;
const double u4U2 = prims->vy*u4U0;
const double u4U3 = prims->vz*u4U0;
const double BU0 = prims->Bx;
const double BU1 = prims->By;
const double BU2 = prims->Bz;
const double P = prims->press;
const double rhob = prims->rho;
const double alpha = metric->lapse;
const double betaU0 = metric->betax;
const double betaU1 = metric->betay;
const double betaU2 = metric->betaz;
const double gammaDD00 = metric->adm_gxx;
const double gammaDD01 = metric->adm_gxy;
const double gammaDD02 = metric->adm_gxz;
const double gammaDD11 = metric->adm_gyy;
const double gammaDD12 = metric->adm_gyz;
const double gammaDD22 = metric->adm_gzz;
const double KDD00 = curv->Kxx;
const double KDD01 = curv->Kxy;
const double KDD02 = curv->Kxz;
const double KDD11 = curv->Kyy;
const double KDD12 = curv->Kyz;
const double KDD22 = curv->Kzz;
const double alpha_dD0 = metric_derivs->lapse[0];
const double alpha_dD1 = metric_derivs->lapse[1];
const double alpha_dD2 = metric_derivs->lapse[2];
  const double _Integer_2 = 2.0;
  const double _Rational_1_2 = 1.0/2.0;
  const double tmp_3 = _Integer_2*gammaDD01*gammaDD02*gammaDD12 + gammaDD00*gammaDD11*gammaDD22 - gammaDD00*((gammaDD12)*(gammaDD12)) - ((gammaDD01)*(gammaDD01))*gammaDD22 - ((gammaDD02)*(gammaDD02))*gammaDD11;
  const double tmp_4 = betaU0*gammaDD00 + betaU1*gammaDD01 + betaU2*gammaDD02;
  const double tmp_5 = betaU0*gammaDD01 + betaU1*gammaDD11 + betaU2*gammaDD12;
  const double tmp_6 = betaU0*gammaDD02 + betaU1*gammaDD12 + betaU2*gammaDD22;
  const double tmp_7 = BU0*(gammaDD00*u4U1 + gammaDD01*u4U2 + gammaDD02*u4U3 + tmp_4*u4U0) + BU1*(gammaDD01*u4U1 + gammaDD11*u4U2 + gammaDD12*u4U3 + tmp_5*u4U0) + BU2*(gammaDD02*u4U1 + gammaDD12*u4U2 + gammaDD22*u4U3 + tmp_6*u4U0);
  const double tmp_9 = (1.0/((alpha)*(alpha)));
  const double tmp_10 = tmp_9/((SQRT_4_PI)*(SQRT_4_PI));
  const double tmp_11 = tmp_10*((tmp_7)*(tmp_7));
  const double tmp_12 = BU0 + tmp_7*u4U1;
  const double tmp_14 = tmp_10/((u4U0)*(u4U0));
  const double tmp_15 = ((tmp_12)*(tmp_12))*tmp_14;
  const double tmp_16 = BU1 + tmp_7*u4U2;
  const double tmp_17 = tmp_14*((tmp_16)*(tmp_16));
  const double tmp_18 = BU2 + tmp_7*u4U3;
  const double tmp_19 = tmp_14*((tmp_18)*(tmp_18));
  const double tmp_20 = gammaDD00*tmp_15 + gammaDD11*tmp_17 + gammaDD22*tmp_19 + tmp_11*(-((alpha)*(alpha)) + betaU0*tmp_4 + betaU1*tmp_5 + betaU2*tmp_6);
  const double tmp_22 = tmp_12*tmp_14*tmp_16;
  const double tmp_23 = tmp_12*tmp_14*tmp_18;
  const double tmp_24 = tmp_14*tmp_16*tmp_18;
  const double tmp_25 = tmp_10*tmp_7/u4U0;
  const double tmp_29 = gammaDD01*tmp_22 + gammaDD02*tmp_23 + gammaDD12*tmp_24 + tmp_12*tmp_25*tmp_4 + tmp_16*tmp_25*tmp_5 + tmp_18*tmp_25*tmp_6;
  const double tmp_30 = P + _Rational_1_2*tmp_20 + tmp_29;
  const double tmp_31 = tmp_30*tmp_9;
  const double tmp_32 = _Integer_2*tmp_29 + h*rhob + tmp_20;
  const double tmp_33 = -tmp_11 - tmp_31 + tmp_32*((u4U0)*(u4U0));
  const double tmp_34 = betaU0*tmp_33;
  const double tmp_35 = tmp_32*u4U0;
  const double tmp_36 = betaU0*tmp_31 - tmp_12*tmp_25 + tmp_35*u4U1;
  const double tmp_38 = betaU1*tmp_31 - tmp_16*tmp_25 + tmp_35*u4U2;
  const double tmp_39 = betaU2*tmp_31 - tmp_18*tmp_25 + tmp_35*u4U3;
  const double tmp_41 = (1.0/(tmp_3));
  const double tmp_42 = _Integer_2*betaU0;
  const double tmp_49 = betaU1*tmp_34 - tmp_22 + tmp_30*(-betaU0*betaU1*tmp_9 + tmp_41*(-gammaDD01*gammaDD22 + gammaDD02*gammaDD12)) + tmp_32*u4U1*u4U2;
  const double tmp_51 = betaU2*tmp_34 - tmp_23 + tmp_30*(-betaU0*betaU2*tmp_9 + tmp_41*(gammaDD01*gammaDD12 - gammaDD02*gammaDD11)) + tmp_32*u4U1*u4U3;
  const double tmp_52 = betaU1*betaU2*tmp_33 - tmp_24 + tmp_30*(-betaU1*betaU2*tmp_9 + tmp_41*(-gammaDD00*gammaDD12 + gammaDD01*gammaDD02)) + tmp_32*u4U2*u4U3;
  cons->tau = alpha*sqrt(tmp_3)*(KDD00*(((betaU0)*(betaU0))*tmp_33 - tmp_15 + tmp_30*(-((betaU0)*(betaU0))*tmp_9 + tmp_41*(gammaDD11*gammaDD22 - ((gammaDD12)*(gammaDD12)))) + tmp_32*((u4U1)*(u4U1)) + tmp_36*tmp_42) + KDD01*(tmp_38*tmp_42 + tmp_49) + KDD01*(_Integer_2*betaU1*tmp_36 + tmp_49) + KDD02*(tmp_39*tmp_42 + tmp_51) + KDD02*(_Integer_2*betaU2*tmp_36 + tmp_51) + KDD11*(_Integer_2*betaU1*tmp_38 + ((betaU1)*(betaU1))*tmp_33 - tmp_17 + tmp_30*(-((betaU1)*(betaU1))*tmp_9 + tmp_41*(gammaDD00*gammaDD22 - ((gammaDD02)*(gammaDD02)))) + tmp_32*((u4U2)*(u4U2))) + KDD12*(_Integer_2*betaU1*tmp_39 + tmp_52) + KDD12*(_Integer_2*betaU2*tmp_38 + tmp_52) + KDD22*(_Integer_2*betaU2*tmp_39 + ((betaU2)*(betaU2))*tmp_33 - tmp_19 + tmp_30*(-((betaU2)*(betaU2))*tmp_9 + tmp_41*(gammaDD00*gammaDD11 - ((gammaDD01)*(gammaDD01)))) + tmp_32*((u4U3)*(u4U3))) + alpha_dD0*(-tmp_34 - tmp_36) + alpha_dD1*(-betaU1*tmp_33 - tmp_38) + alpha_dD2*(-betaU2*tmp_33 - tmp_39));
}
}
