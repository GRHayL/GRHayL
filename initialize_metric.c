#include "con2prim_header.h"

/* This function fills the struct metric_quantities with data. The g**
   arguments are the ADM variables. For more information on the
   arguments, see the definition of the struct in new_header.h. */

//We could consider allowing the struct to be either BSSN or ADM in the future if
//we ever need a BSSN struct by having initialize_BSSN and initialize_ADM functions.
void initialize_metric(metric_quantities *restrict metric, const double lapse,
                const double gxx, const double gxy, const double gxz,
                const double gyy, const double gyz, const double gzz,
                const double betax, const double betay, const double betaz) {

  metric->adm_gxx = gxx;
  metric->adm_gxy = gxy;
  metric->adm_gxz = gxz;
  metric->adm_gyy = gyy;
  metric->adm_gyz = gyz;
  metric->adm_gzz = gzz;
  metric->lapse = lapse;
  metric->betax = betax;
  metric->betay = betay;
  metric->betaz = betaz;

  double gijdet = fabs(gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz
                     - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz);

  // This is the analytic algorithm for finding the inverse of a 3x3 matrix. See e.g.
  // https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices
  metric->adm_gupxx =  gijdet*( gyy * gzz - gyz * gyz );
  metric->adm_gupxy = -gijdet*( gxy * gzz - gyz * gxz );
  metric->adm_gupxz =  gijdet*( gxy * gyz - gyy * gxz );
  metric->adm_gupyy =  gijdet*( gxx * gzz - gxz * gxz );
  metric->adm_gupyz = -gijdet*( gxx * gyz - gxy * gxz );
  metric->adm_gupzz =  gijdet*( gxx * gyy - gxy * gxy );

  double phi = (1.0/12.0) * log(gijdet);

  metric->psi2 = exp(2.0*phi);
  metric->psi4 = SQR(metric->psi2);
  metric->psi6 = metric->psi4*metric->psi2;
  metric->psi4inv = 1.0/metric->psi4;
  metric->lapseinv = 1.0/metric->lapse;
  metric->lapseinv2=SQR(metric->lapseinv);

  double shift_xL = metric->adm_gxx*metric->betax + metric->adm_gxy*metric->betay + metric->adm_gxz*metric->betaz;
  double shift_yL = metric->adm_gxy*metric->betax + metric->adm_gyy*metric->betay + metric->adm_gyz*metric->betaz;
  double shift_zL = metric->adm_gxz*metric->betax + metric->adm_gyz*metric->betay + metric->adm_gzz*metric->betaz;
  double beta2L   = shift_xL*metric->betax + shift_yL*metric->betay + shift_zL*metric->betaz;

  metric->g4dn[0][0]                      = -SQR(metric->lapse) + beta2L;
  metric->g4dn[0][1] = metric->g4dn[1][0] = shift_xL;
  metric->g4dn[0][2] = metric->g4dn[2][0] = shift_yL;
  metric->g4dn[0][3] = metric->g4dn[3][0] = shift_zL;
  metric->g4dn[1][1]                      = metric->adm_gxx;
  metric->g4dn[1][2] = metric->g4dn[2][1] = metric->adm_gxy;
  metric->g4dn[1][3] = metric->g4dn[3][1] = metric->adm_gxz;
  metric->g4dn[2][2]                      = metric->adm_gyy;
  metric->g4dn[2][3] = metric->g4dn[3][2] = metric->adm_gyz;
  metric->g4dn[3][3]                      = metric->adm_gzz;

  metric->g4up[0][0]                      = -metric->lapseinv2;
  metric->g4up[0][1] = metric->g4up[1][0] = metric->betax*metric->lapseinv2;
  metric->g4up[0][2] = metric->g4up[2][0] = metric->betay*metric->lapseinv2;
  metric->g4up[0][3] = metric->g4up[3][0] = metric->betaz*metric->lapseinv2;
  metric->g4up[1][1]                      = metric->adm_gupxx - metric->betax*metric->betax*metric->lapseinv2;
  metric->g4up[1][2] = metric->g4up[2][1] = metric->adm_gupxy - metric->betax*metric->betay*metric->lapseinv2;
  metric->g4up[1][3] = metric->g4up[3][1] = metric->adm_gupxz - metric->betax*metric->betaz*metric->lapseinv2;
  metric->g4up[2][2]                      = metric->adm_gupyy - metric->betay*metric->betay*metric->lapseinv2;
  metric->g4up[2][3] = metric->g4up[3][2] = metric->adm_gupyz - metric->betay*metric->betaz*metric->lapseinv2;
  metric->g4up[3][3]                      = metric->adm_gupzz - metric->betaz*metric->betaz*metric->lapseinv2;
}
