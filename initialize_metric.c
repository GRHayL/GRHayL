#include "con2prim_header.h"

/* This function fills the struct metric_quantities with data. The phi, psi,
   and g** arguments are the BSSN variables. For more information on the
   arguments, see the definition of the struct in new_header.h. */
void initialize_metric(metric_quantities *restrict metric, 
                const double phi, const double psi, const double lapse,
                const double gxx, const double gxy, const double gxz,
                const double gyy, const double gyz, const double gzz,
                const double betax, const double betay, const double betaz) {
  metric->bssn_phi = phi;
  metric->bssn_psi = psi;
  metric->bssn_gxx = gxx;
  metric->bssn_gxy = gxy;
  metric->bssn_gxz = gxz;
  metric->bssn_gyy = gyy;
  metric->bssn_gyz = gyz;
  metric->bssn_gzz = gzz;
  metric->lapse = lapse;
  metric->betax = betax;
  metric->betay = betay;
  metric->betaz = betaz;
//TODO: the constraint getgij=1 could be applied inside this initialization
//      automatically, possibly with a parameter to turn it off.
  metric->bssn_gupxx =   ( gyy * gzz - gyz * gyz );
  metric->bssn_gupxy = - ( gxy * gzz - gyz * gxz );
  metric->bssn_gupxz =   ( gxy * gyz - gyy * gxz );
  metric->bssn_gupyy =   ( gxx * gzz - gxz * gxz );
  metric->bssn_gupyz = - ( gxx * gyz - gxy * gxz );
  metric->bssn_gupzz =   ( gxx * gyy - gxy * gxy );
//Do we really need lapm1?
  metric->lapm1 = lapse-1.0;
  metric->psi2 = exp(2.0*metric->bssn_phi);
  metric->psi4 = SQR(metric->psi2);
  metric->psi6 = metric->psi4*metric->psi2;
  metric->psi4inv = 1.0/metric->psi4;
  metric->lapseinv = 1.0/metric->lapse;
  metric->lapseinv2=SQR(metric->lapseinv);
  metric->adm_gxx = metric->bssn_gxx*metric->psi4;
  metric->adm_gxy = metric->bssn_gxy*metric->psi4;
  metric->adm_gxz = metric->bssn_gxz*metric->psi4;
  metric->adm_gyy = metric->bssn_gyy*metric->psi4;
  metric->adm_gyz = metric->bssn_gyz*metric->psi4;
  metric->adm_gzz = metric->bssn_gzz*metric->psi4;
  metric->adm_gupxx = metric->bssn_gupxx*metric->psi4inv;
  metric->adm_gupxy = metric->bssn_gupxy*metric->psi4inv;
  metric->adm_gupxz = metric->bssn_gupxz*metric->psi4inv;
  metric->adm_gupyy = metric->bssn_gupyy*metric->psi4inv;
  metric->adm_gupyz = metric->bssn_gupyz*metric->psi4inv;
  metric->adm_gupzz = metric->bssn_gupzz*metric->psi4inv;

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
