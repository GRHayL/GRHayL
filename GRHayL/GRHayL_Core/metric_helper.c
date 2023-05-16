#include "grhayl.h"

/**********************************************************************
    raise_vector_4D() and raise_vector_3D():

         -- calculates the contravariant form of a covariant vector
            using the inverse of the metric
              v^\alpha = g^{\alpha\beta} v_\beta
            3D and 4D vectors are currently supported.
***********************************************************************/
void raise_vector_4D(
      const metric_quantities *restrict metric,
      const double vcov[4],
      double vcon[4]) {

  for(int i=0; i<4; i++) {
    vcon[i] = 0.0;
    for(int j=0; j<4; j++)
      vcon[i] += metric->g4up[i][j]*vcov[j];
  }
}

void raise_vector_3D(
      const metric_quantities *restrict metric,
      const double vcov[3],
      double vcon[3]) {

  vcon[0] = metric->adm_gupxx * vcov[0]
          + metric->adm_gupxy * vcov[1]
          + metric->adm_gupxz * vcov[2];

  vcon[1] = metric->adm_gupxy * vcov[0]
          + metric->adm_gupyy * vcov[1]
          + metric->adm_gupyz * vcov[2];

  vcon[2] = metric->adm_gupxz * vcov[0]
          + metric->adm_gupyz * vcov[1]
          + metric->adm_gupzz * vcov[2];
}

/**********************************************************************
     lower_vector():

          -- calculates the covariant form of a contravariant vector
             using the metric
              v_\alpha = g_{\alpha\beta} v^\beta
***********************************************************************/
void lower_vector(
      const metric_quantities *restrict metric,
      const double vcon[4],
      double vcov[4]) {

  for(int i=0; i<4; i++) {
    vcov[i] = 0.0;
    for(int j=0; j<4; j++)
      vcov[i] += metric->g4dn[i][j]*vcon[j];
  }
  return ;
}

/**********************************************************************
     square_Vup():

          -- calculates the magnitude squared of a contravariant 3-vector
             using the ADM metric
              v^2 = g_{\alpha\beta} v^\alpha v^\beta
***********************************************************************/
double compute_vec2_from_vcov(
      const metric_quantities *restrict metric,
      const double *restrict vcov) {

  return metric->adm_gupxx * vcov[0] * vcov[0] +
         metric->adm_gupyy * vcov[1] * vcov[1] +
         metric->adm_gupzz * vcov[2] * vcov[2] +
   2.0 * metric->adm_gupxy * vcov[0] * vcov[1] +
   2.0 * metric->adm_gupxz * vcov[0] * vcov[2] +
   2.0 * metric->adm_gupyz * vcov[1] * vcov[2];
}

/**********************************************************************
     square_Vdn():

          -- calculates the magnitude squared of a contravariant 3-vector
             using the ADM metric
              v^2 = g_{\alpha\beta} v^\alpha v^\beta
***********************************************************************/
double compute_vec2_from_vcon(
      const metric_quantities *restrict metric,
      const double *restrict vcon) {

  double vsq = 0.0;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      vsq += metric->g4dn[i+1][j+1] * vcon[i] * vcon[j];

  return vsq;
}
