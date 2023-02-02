#include "GRHayL.h"

/**********************************************************************
    raise_vector():

         -- calculates the contravariant form of a covariant vector
            using the inverse of the metric
              v^\alpha = g^{\alpha\beta} v_\beta
***********************************************************************/
void raise_vector(
      const metric_quantities *restrict metric,
      const double vcov[4],
      double vcon[4]) {

  for(int i=0; i<4; i++) {
    vcon[i] = 0.0;
    for(int j=0; j<4; j++)
      vcon[i] += metric->g4up[i][j]*vcov[j];
  }
  return ;
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
