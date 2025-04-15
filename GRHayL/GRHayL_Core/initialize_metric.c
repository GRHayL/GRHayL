#include "ghl.h"

/**
 * @ingroup pack_struct
 * @brief Initialize the metric struct from user input
 *
 * @details
 * This function takes pointwise information about the spacetime and uses it
 * to initialize every element of the given ghl_metric_quantities struct.
 * It also computes all the additional struct elements which are derived
 * from these inputs.
 *
 * @param[in] alpha: lapse \f$ \alpha \f$
 *
 * @param[in] betax, betay, betaz: components of the shift \f$ \beta^i \f$
 *
 * @param[in] gxx, gxy, gxz, gyy, gyz, gzz: components of the 3-metric \f$ g_{i j} \f$
 *
 * @param[out] metric: pointer to ghl_metric_quantities struct
 *
 * @returns void
 */
void ghl_initialize_metric(
      const double alpha,
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

  metric->lapse                                 = alpha;
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
