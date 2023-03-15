#include "flux_source.h"
/*
 * Add extrinsic curvature source term for tau_tilde
 */
void calculate_tau_tilde_source_term_extrinsic_curv(const primitive_quantities *restrict prims, struct eos_parameters const *restrict eos, const metric_quantities *restrict metric, const extrinsic_curvature *restrict curv, const metric_derivatives *restrict metric_derivs, conservative_quantities *restrict cons) {

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
  const double _Integer_2 = 2.0;
  const double _Rational_1_2 = 1.0/2.0;
  const double tmp_3 = _Integer_2*gammaDD01*gammaDD02*gammaDD12 + gammaDD00*gammaDD11*gammaDD22 - gammaDD00*((gammaDD12)*(gammaDD12)) - ((gammaDD01)*(gammaDD01))*gammaDD22 - ((gammaDD02)*(gammaDD02))*gammaDD11;
  const double tmp_4 = betaU0*gammaDD00 + betaU1*gammaDD01 + betaU2*gammaDD02;
  const double tmp_5 = betaU0*gammaDD01 + betaU1*gammaDD11 + betaU2*gammaDD12;
  const double tmp_6 = betaU0*gammaDD02 + betaU1*gammaDD12 + betaU2*gammaDD22;
  const double tmp_7 = BU0*(gammaDD00*u4U1 + gammaDD01*u4U2 + gammaDD02*u4U3 + tmp_4*u4U0) + BU1*(gammaDD01*u4U1 + gammaDD11*u4U2 + gammaDD12*u4U3 + tmp_5*u4U0) + BU2*(gammaDD02*u4U1 + gammaDD12*u4U2 + gammaDD22*u4U3 + tmp_6*u4U0);
  const double tmp_8 = BU0 + tmp_7*u4U1;
  const double tmp_11 = (1.0/((alpha)*(alpha)));
  const double tmp_12 = tmp_11/((SQRT_4_PI)*(SQRT_4_PI));
  const double tmp_13 = tmp_12/((u4U0)*(u4U0));
  const double tmp_14 = tmp_13*((tmp_8)*(tmp_8));
  const double tmp_15 = BU1 + tmp_7*u4U2;
  const double tmp_17 = tmp_13*tmp_15*tmp_8;
  const double tmp_18 = BU2 + tmp_7*u4U3;
  const double tmp_19 = tmp_13*tmp_18*tmp_8;
  const double tmp_20 = tmp_13*tmp_15*tmp_18;
  const double tmp_21 = tmp_12*tmp_7/u4U0;
  const double tmp_25 = gammaDD01*tmp_17 + gammaDD02*tmp_19 + gammaDD12*tmp_20 + tmp_15*tmp_21*tmp_5 + tmp_18*tmp_21*tmp_6 + tmp_21*tmp_4*tmp_8;
  const double tmp_26 = tmp_13*((tmp_15)*(tmp_15));
  const double tmp_27 = tmp_13*((tmp_18)*(tmp_18));
  const double tmp_28 = tmp_12*((tmp_7)*(tmp_7));
  const double tmp_29 = gammaDD00*tmp_14 + gammaDD11*tmp_26 + gammaDD22*tmp_27 + tmp_28*(-((alpha)*(alpha)) + betaU0*tmp_4 + betaU1*tmp_5 + betaU2*tmp_6);
  const double tmp_30 = _Integer_2*tmp_25 + h*rhob + tmp_29;
  const double tmp_32 = (1.0/(tmp_3));
  const double tmp_33 = P + _Rational_1_2*tmp_29 + tmp_25;
  const double tmp_34 = tmp_11*tmp_33;
  const double tmp_35 = -tmp_28 + tmp_30*((u4U0)*(u4U0)) - tmp_34;
  const double tmp_36 = tmp_30*u4U0;
  const double tmp_37 = betaU0*tmp_34 - tmp_21*tmp_8 + tmp_36*u4U1;
  const double tmp_38 = _Integer_2*betaU0;
  const double tmp_40 = betaU1*tmp_34 - tmp_15*tmp_21 + tmp_36*u4U2;
  const double tmp_43 = betaU2*tmp_34 - tmp_18*tmp_21 + tmp_36*u4U3;
  const double tmp_47 = betaU0*betaU1*tmp_35 - tmp_17 + tmp_30*u4U1*u4U2 + tmp_33*(-betaU0*betaU1*tmp_11 + tmp_32*(-gammaDD01*gammaDD22 + gammaDD02*gammaDD12));
  const double tmp_50 = betaU0*betaU2*tmp_35 - tmp_19 + tmp_30*u4U1*u4U3 + tmp_33*(-betaU0*betaU2*tmp_11 + tmp_32*(gammaDD01*gammaDD12 - gammaDD02*gammaDD11));
  const double tmp_52 = betaU1*betaU2*tmp_35 - tmp_20 + tmp_30*u4U2*u4U3 + tmp_33*(-betaU1*betaU2*tmp_11 + tmp_32*(-gammaDD00*gammaDD12 + gammaDD01*gammaDD02));
  cons->tau = alpha*sqrt(tmp_3)*(KDD00*(((betaU0)*(betaU0))*tmp_35 - tmp_14 + tmp_30*((u4U1)*(u4U1)) + tmp_33*(-((betaU0)*(betaU0))*tmp_11 + tmp_32*(gammaDD11*gammaDD22 - ((gammaDD12)*(gammaDD12)))) + tmp_37*tmp_38) + KDD01*(tmp_38*tmp_40 + tmp_47) + KDD01*(_Integer_2*betaU1*tmp_37 + tmp_47) + KDD02*(tmp_38*tmp_43 + tmp_50) + KDD02*(_Integer_2*betaU2*tmp_37 + tmp_50) + KDD11*(_Integer_2*betaU1*tmp_40 + ((betaU1)*(betaU1))*tmp_35 - tmp_26 + tmp_30*((u4U2)*(u4U2)) + tmp_33*(-((betaU1)*(betaU1))*tmp_11 + tmp_32*(gammaDD00*gammaDD22 - ((gammaDD02)*(gammaDD02))))) + KDD12*(_Integer_2*betaU1*tmp_43 + tmp_52) + KDD12*(_Integer_2*betaU2*tmp_40 + tmp_52) + KDD22*(_Integer_2*betaU2*tmp_43 + ((betaU2)*(betaU2))*tmp_35 - tmp_27 + tmp_30*((u4U3)*(u4U3)) + tmp_33*(-((betaU2)*(betaU2))*tmp_11 + tmp_32*(gammaDD00*gammaDD11 - ((gammaDD01)*(gammaDD01))))));
}
}
