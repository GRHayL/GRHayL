#include "grhayl.h"

/**********************************************************************
    ghl_raise_vector_4D() and ghl_raise_vector_3D():

         -- calculates the contravariant form of a covariant vector
            using the inverse of the metric
              v^\alpha = g^{\alpha\beta} v_\beta
            3D and 4D vectors are currently supported.
***********************************************************************/
void ghl_raise_vector_4D(
      const double g4UU[4][4],
      const double vecD[4],
      double vecU[4]) {

  for(int i=0; i<4; i++) {
    vecU[i] = 0.0;
    for(int j=0; j<4; j++)
      vecU[i] += g4UU[i][j]*vecD[j];
  }
}

void ghl_raise_vector_3D(
      const double gammaUU[3][3],
      const double vecD[3],
      double vecU[3]) {

  vecU[0] = gammaUU[0][0] * vecD[0]
          + gammaUU[0][1] * vecD[1]
          + gammaUU[0][2] * vecD[2];

  vecU[1] = gammaUU[0][1] * vecD[0]
          + gammaUU[1][1] * vecD[1]
          + gammaUU[1][2] * vecD[2];

  vecU[2] = gammaUU[0][2] * vecD[0]
          + gammaUU[1][2] * vecD[1]
          + gammaUU[2][2] * vecD[2];
}

/**********************************************************************

    ghl_lower_vector_4D() and ghl_lower_vector_3D():

         -- calculates the covariant form of a contravariant vector
            using the metric
              v_\alpha = g_{\alpha\beta} v^\beta
            3D and 4D vectors are currently supported.
***********************************************************************/
void ghl_lower_vector_4D(
      const double g4DD[4][4],
      const double vecU[4],
      double vecD[4]) {

  for(int i=0; i<4; i++) {
    vecD[i] = 0.0;
    for(int j=0; j<4; j++)
      vecD[i] += g4DD[i][j]*vecU[j];
  }
}

void ghl_lower_vector_3D(
      const double gammaDD[3][3],
      const double vecU[3],
      double vecD[3]) {

  vecD[0] = gammaDD[0][0] * vecU[0]
          + gammaDD[0][1] * vecU[1]
          + gammaDD[0][2] * vecU[2];

  vecD[1] = gammaDD[0][1] * vecU[0]
          + gammaDD[1][1] * vecU[1]
          + gammaDD[1][2] * vecU[2];

  vecD[2] = gammaDD[0][2] * vecU[0]
          + gammaDD[1][2] * vecU[1]
          + gammaDD[2][2] * vecU[2];
}

/**********************************************************************
     square_Vup():

          -- calculates the magnitude squared of a contravariant 3-vector
             using the ADM metric
              v^2 = g_{\alpha\beta} v^\alpha v^\beta
***********************************************************************/
double ghl_compute_vec2_from_vecD(
      const double gammaUU[3][3],
      const double *restrict vecD) {

  return gammaUU[0][0] * vecD[0] * vecD[0] +
         gammaUU[1][1] * vecD[1] * vecD[1] +
         gammaUU[2][2] * vecD[2] * vecD[2] +
   2.0 * gammaUU[0][1] * vecD[0] * vecD[1] +
   2.0 * gammaUU[0][2] * vecD[0] * vecD[2] +
   2.0 * gammaUU[1][2] * vecD[1] * vecD[2];
}

/**********************************************************************
     square_Vdn():

          -- calculates the magnitude squared of a contravariant 3-vector
             using the ADM metric
              v^2 = g_{\alpha\beta} v^\alpha v^\beta
***********************************************************************/
double ghl_compute_vec2_from_vecU(
      const double gammaDD[3][3],
      const double *restrict vecU) {

  double vsq = 0.0;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      vsq += gammaDD[i][j] * vecU[i] * vecU[j];

  return vsq;

  //return gammaDD[0][0] * vecU[0] * vecU[0] +
  //       gammaDD[1][1] * vecU[1] * vecU[1] +
  //       gammaDD[2][2] * vecU[2] * vecU[2] +
  // 2.0 * gammaDD[0][1] * vecU[0] * vecU[1] +
  // 2.0 * gammaDD[0][2] * vecU[0] * vecU[2] +
  // 2.0 * gammaDD[1][2] * vecU[1] * vecU[2];
}
