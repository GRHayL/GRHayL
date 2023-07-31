#include "ghl.h"

/**********************************************************************
    ghl_raise_lower_vector_4D() and ghl_raise_lower_vector_3D():

         -- calculates the contravariant form of a covariant vector
            using the inverse of the metric
              v^\alpha = g^{\alpha\beta} v_\beta
	    or calculates the covariant form of a contravariant vector
	    using the metric
              v_\alpha = g_{\alpha\beta} v^\beta
            3D and 4D vectors are currently supported.
***********************************************************************/
void ghl_raise_lower_vector_4D(
      const double g4[4][4],
      const double vec[4],
      double vec_inv[4]) {

  for(int i=0; i<4; i++) {
    vec_inv[i] = 0.0;
    for(int j=0; j<4; j++)
      vec_inv[i] += g4[i][j]*vec[j];
  }
}

void ghl_raise_lower_vector_3D(
      const double gamma[3][3],
      const double vec[3],
      double vec_inv[3]) {

  vec_inv[0] = gamma[0][0] * vec[0]
             + gamma[0][1] * vec[1]
             + gamma[0][2] * vec[2];

  vec_inv[1] = gamma[0][1] * vec[0]
             + gamma[1][1] * vec[1]
             + gamma[1][2] * vec[2];

  vec_inv[2] = gamma[0][2] * vec[0]
             + gamma[1][2] * vec[1]
             + gamma[2][2] * vec[2];
}

/**********************************************************************
     square_Vup():

          -- calculates the magnitude squared of a contravariant 3-vector
             using the ADM metric
              v^2 = g_{\alpha\beta} v^\alpha v^\beta
***********************************************************************/
double ghl_compute_vec2_from_vec4D(
      const double g4[4][4],
      const double vec[4]) {

  return g4[0][0] * vec[0] * vec[0] +
         g4[1][1] * vec[1] * vec[1] +
         g4[2][2] * vec[2] * vec[2] +
         g4[3][3] * vec[3] * vec[3] +
  2.0 * (g4[0][1] * vec[0] * vec[1] +
         g4[0][2] * vec[0] * vec[2] +
         g4[0][3] * vec[0] * vec[3] +
         g4[1][2] * vec[1] * vec[2] +
         g4[1][3] * vec[1] * vec[3] +
         g4[2][3] * vec[2] * vec[3]);
}

double ghl_compute_vec2_from_vecD(
      const double gammaUU[3][3],
      const double vecD[3]) {

  return gammaUU[0][0] * vecD[0] * vecD[0] +
         gammaUU[1][1] * vecD[1] * vecD[1] +
         gammaUU[2][2] * vecD[2] * vecD[2] +
  2.0 * (gammaUU[0][1] * vecD[0] * vecD[1] +
         gammaUU[0][2] * vecD[0] * vecD[2] +
         gammaUU[1][2] * vecD[1] * vecD[2]);
}

/**********************************************************************
     square_Vdn():

          -- calculates the magnitude squared of a contravariant 3-vector
             using the ADM metric
              v^2 = g_{\alpha\beta} v^\alpha v^\beta
***********************************************************************/
double ghl_compute_vec2_from_vecU(
      const double gammaDD[3][3],
      const double vecU[3]) {

  double vsq = 0.0;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      vsq += gammaDD[i][j] * vecU[i] * vecU[j];

  return vsq;

  //return gammaDD[0][0] * vecU[0] * vecU[0] +
  //       gammaDD[1][1] * vecU[1] * vecU[1] +
  //       gammaDD[2][2] * vecU[2] * vecU[2] +
  //2.0 * (gammaDD[0][1] * vecU[0] * vecU[1] +
  //       gammaDD[0][2] * vecU[0] * vecU[2] +
  //       gammaDD[1][2] * vecU[1] * vecU[2]);
}
