#include "flux_source.h"
/*
 * Add source term for 2-component of Stilde and tau_tilde
 */
void calculate_source_terms_dirn2(const primitive_quantities *restrict prims, struct eos_parameters const *restrict eos, const metric_quantities *restrict metric, const metric_quantities *restrict metric_derivs, conservative_quantities *restrict cons) {

{

double h, cs2;

eos->compute_h_and_cs2(eos, prims, &h, &cs2);
const double u4U0 = prims->u0;
const double u4U1 = prims->vU[0]*u4U0;
const double u4U2 = prims->vU[1]*u4U0;
const double u4U3 = prims->vU[2]*u4U0;
const double BU0 = prims->BU[0];
const double BU1 = prims->BU[1];
const double BU2 = prims->BU[2];
const double P = prims->press;
const double rhob = prims->rho;
const double alpha = metric->lapse;
const double betaU0 = metric->betaU[0];
const double betaU1 = metric->betaU[1];
const double betaU2 = metric->betaU[2];
const double gammaDD00 = metric->gammaDD[0][0];
const double gammaDD01 = metric->gammaDD[0][1];
const double gammaDD02 = metric->gammaDD[0][2];
const double gammaDD11 = metric->gammaDD[1][1];
const double gammaDD12 = metric->gammaDD[1][2];
const double gammaDD22 = metric->gammaDD[2][2];
const double alpha_dD2 = metric_derivs[2].lapse;
const double betaU_dD02 = metric_derivs[2].betaU[0];
const double betaU_dD12 = metric_derivs[2].betaU[1];
const double betaU_dD22 = metric_derivs[2].betaU[2];
const double gammaDD_dD002 = metric_derivs[2].gammaDD[0][0];
const double gammaDD_dD012 = metric_derivs[2].gammaDD[0][1];
const double gammaDD_dD022 = metric_derivs[2].gammaDD[0][2];
const double gammaDD_dD112 = metric_derivs[2].gammaDD[1][1];
const double gammaDD_dD122 = metric_derivs[2].gammaDD[1][2];
const double gammaDD_dD222 = metric_derivs[2].gammaDD[2][2];
  const double _Integer_2 = 2.0;
  const double _Rational_1_2 = 1.0/2.0;
  const double tmp_0 = betaU0*gammaDD_dD002 + betaU1*gammaDD_dD012 + betaU2*gammaDD_dD022 + betaU_dD02*gammaDD00 + betaU_dD12*gammaDD01 + betaU_dD22*gammaDD02;
  const double tmp_1 = betaU0*gammaDD00 + betaU1*gammaDD01 + betaU2*gammaDD02;
  const double tmp_2 = betaU0*gammaDD01 + betaU1*gammaDD11 + betaU2*gammaDD12;
  const double tmp_3 = betaU0*gammaDD02 + betaU1*gammaDD12 + betaU2*gammaDD22;
  const double tmp_4 = BU0*(gammaDD00*u4U1 + gammaDD01*u4U2 + gammaDD02*u4U3 + tmp_1*u4U0) + BU1*(gammaDD01*u4U1 + gammaDD11*u4U2 + gammaDD12*u4U3 + tmp_2*u4U0) + BU2*(gammaDD02*u4U1 + gammaDD12*u4U2 + gammaDD22*u4U3 + tmp_3*u4U0);
  const double tmp_5 = BU0 + tmp_4*u4U1;
  const double tmp_7 = (1.0/((alpha)*(alpha)));
  const double tmp_8 = tmp_7/((SQRT_4_PI)*(SQRT_4_PI));
  const double tmp_9 = tmp_4*tmp_8/u4U0;
  const double tmp_12 = tmp_8/((u4U0)*(u4U0));
  const double tmp_13 = tmp_12*((tmp_5)*(tmp_5));
  const double tmp_14 = BU1 + tmp_4*u4U2;
  const double tmp_15 = tmp_12*((tmp_14)*(tmp_14));
  const double tmp_16 = BU2 + tmp_4*u4U3;
  const double tmp_17 = tmp_12*((tmp_16)*(tmp_16));
  const double tmp_18 = ((tmp_4)*(tmp_4))*tmp_8;
  const double tmp_19 = gammaDD00*tmp_13 + gammaDD11*tmp_15 + gammaDD22*tmp_17 + tmp_18*(-((alpha)*(alpha)) + betaU0*tmp_1 + betaU1*tmp_2 + betaU2*tmp_3);
  const double tmp_21 = tmp_12*tmp_14*tmp_5;
  const double tmp_22 = tmp_12*tmp_16*tmp_5;
  const double tmp_23 = tmp_12*tmp_14*tmp_16;
  const double tmp_26 = gammaDD01*tmp_21 + gammaDD02*tmp_22 + gammaDD12*tmp_23 + tmp_1*tmp_5*tmp_9 + tmp_14*tmp_2*tmp_9 + tmp_16*tmp_3*tmp_9;
  const double tmp_27 = P + _Rational_1_2*tmp_19 + tmp_26;
  const double tmp_28 = tmp_27*tmp_7;
  const double tmp_29 = _Integer_2*tmp_26 + h*rhob + tmp_19;
  const double tmp_30 = tmp_29*u4U1;
  const double tmp_34 = _Integer_2*gammaDD01*gammaDD02*gammaDD12 + gammaDD00*gammaDD11*gammaDD22 - gammaDD00*((gammaDD12)*(gammaDD12)) - ((gammaDD01)*(gammaDD01))*gammaDD22 - ((gammaDD02)*(gammaDD02))*gammaDD11;
  const double tmp_35 = sqrt(tmp_34);
  const double tmp_36 = alpha*tmp_35;
  const double tmp_37 = betaU0*gammaDD_dD012 + betaU1*gammaDD_dD112 + betaU2*gammaDD_dD122 + betaU_dD02*gammaDD01 + betaU_dD12*gammaDD11 + betaU_dD22*gammaDD12;
  const double tmp_39 = betaU0*gammaDD_dD022 + betaU1*gammaDD_dD122 + betaU2*gammaDD_dD222 + betaU_dD02*gammaDD02 + betaU_dD12*gammaDD12 + betaU_dD22*gammaDD22;
  const double tmp_40 = betaU2*tmp_28 - tmp_16*tmp_9 + tmp_29*u4U0*u4U3;
  const double tmp_42 = (1.0/(tmp_34));
  const double tmp_43 = -tmp_18 - tmp_28 + tmp_29*((u4U0)*(u4U0));
  cons->SD[2] = _Rational_1_2*(gammaDD_dD002*tmp_36*(-tmp_13 + tmp_27*(-((betaU0)*(betaU0))*tmp_7 + tmp_42*(gammaDD11*gammaDD22 - ((gammaDD12)*(gammaDD12)))) + tmp_29*((u4U1)*(u4U1))) + gammaDD_dD112*tmp_36*(-tmp_15 + tmp_27*(-((betaU1)*(betaU1))*tmp_7 + tmp_42*(gammaDD00*gammaDD22 - ((gammaDD02)*(gammaDD02)))) + tmp_29*((u4U2)*(u4U2))) + gammaDD_dD222*tmp_36*(-tmp_17 + tmp_27*(-((betaU2)*(betaU2))*tmp_7 + tmp_42*(gammaDD00*gammaDD11 - ((gammaDD01)*(gammaDD01)))) + tmp_29*((u4U3)*(u4U3))) + tmp_36*tmp_43*(-_Integer_2*alpha*alpha_dD2 + betaU0*tmp_0 + betaU1*tmp_37 + betaU2*tmp_39 + betaU_dD02*tmp_1 + betaU_dD12*tmp_2 + betaU_dD22*tmp_3)) + gammaDD_dD012*tmp_36*(-tmp_21 + tmp_27*(-betaU0*betaU1*tmp_7 + tmp_42*(-gammaDD01*gammaDD22 + gammaDD02*gammaDD12)) + tmp_30*u4U2) + gammaDD_dD022*tmp_36*(-tmp_22 + tmp_27*(-betaU0*betaU2*tmp_7 + tmp_42*(gammaDD01*gammaDD12 - gammaDD02*gammaDD11)) + tmp_30*u4U3) + gammaDD_dD122*tmp_36*(-tmp_23 + tmp_27*(-betaU1*betaU2*tmp_7 + tmp_42*(-gammaDD00*gammaDD12 + gammaDD01*gammaDD02)) + tmp_29*u4U2*u4U3) + tmp_0*tmp_36*(betaU0*tmp_28 + tmp_30*u4U0 - tmp_5*tmp_9) + tmp_36*tmp_37*(betaU1*tmp_28 - tmp_14*tmp_9 + tmp_29*u4U0*u4U2) + tmp_36*tmp_39*tmp_40;
  cons->tau = alpha*alpha_dD2*tmp_35*(-betaU2*tmp_43 - tmp_40);
}
}
