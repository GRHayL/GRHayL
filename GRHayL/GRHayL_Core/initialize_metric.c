#include "ghl.h"

/*
 * Function     : ghl_initialize_metric()
 * Description  : Initialize the ghl_metric_quantities struct from user input
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_initialize_metric
*/

GHL_HOST_DEVICE
void ghl_initialize_metric(
      const double lapse,
      const double betax,
      const double betay,
      const double betaz,
      const double gxx,
      const double gxy,
      const double gxz,
      const double gyy,
      const double gyz,
      const double gzz,
      ghl_metric_quantities *restrict metric) {

  metric->lapse                                 = lapse;
  metric->betaU[0]                              = betax;
  metric->betaU[1]                              = betay;
  metric->betaU[2]                              = betaz;
  metric->gammaDD[0][0]                         = gxx;
  metric->gammaDD[1][1]                         = gyy;
  metric->gammaDD[2][2]                         = gzz;
  metric->gammaDD[0][1] = metric->gammaDD[1][0] = gxy;
  metric->gammaDD[0][2] = metric->gammaDD[2][0] = gxz;
  metric->gammaDD[1][2] = metric->gammaDD[2][1] = gyz;

  metric->lapseinv  = 1.0/metric->lapse;
  metric->lapseinv2 = SQR(metric->lapseinv);

  // This is the analytic algorithm for finding the inverse of a 3x3 matrix. See e.g.
  // https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices
  metric->detgamma = fabs(gxx * (gyy*gzz - gyz*gyz) + gxy * (gyz*gxz - gxy*gzz) + gxz * (gxy*gyz - gyy*gxz));
  metric->sqrt_detgamma = sqrt(metric->detgamma);

  metric->gammaUU[0][0]                         =   ( gyy * gzz - gyz * gyz )/metric->detgamma;
  metric->gammaUU[1][1]                         =   ( gxx * gzz - gxz * gxz )/metric->detgamma;
  metric->gammaUU[2][2]                         =   ( gxx * gyy - gxy * gxy )/metric->detgamma;
  metric->gammaUU[0][1] = metric->gammaUU[1][0] = - ( gxy * gzz - gyz * gxz )/metric->detgamma;
  metric->gammaUU[0][2] = metric->gammaUU[2][0] =   ( gxy * gyz - gyy * gxz )/metric->detgamma;
  metric->gammaUU[1][2] = metric->gammaUU[2][1] = - ( gxx * gyz - gxy * gxz )/metric->detgamma;
}
