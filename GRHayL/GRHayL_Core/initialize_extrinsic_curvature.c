#include "grhayl.h"

/* Function    : ghl_initialize_extrinsic_curvature()
 * Description : Initialize the extrinsic_curvature struct from user input
 *
 * Inputs      : Kij            - value of the (i,j) component of the
 *                                extrinsic curvature tensor K^ij
 *
 * Outputs     : curv           - returns extrinsic_curvature struct containing
 *                                the inputs
 */

void ghl_initialize_extrinsic_curvature(
      const double Kxx, const double Kxy, const double Kxz,
      const double Kyy, const double Kyz, const double Kzz,
      extrinsic_curvature *restrict curv){

  curv->K[0][0]                 = Kxx;
  curv->K[1][1]                 = Kyy;
  curv->K[2][2]                 = Kzz;
  curv->K[0][1] = curv->K[1][0] = Kxy;
  curv->K[0][2] = curv->K[2][0] = Kxz;
  curv->K[1][2] = curv->K[2][1] = Kyz;
}
