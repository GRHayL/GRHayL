#include "ghl.h"

/*
 * Function     : ghl_initialize_extrinsic_curvature()
 * Description  : Initialize the ghl_extrinsic_curvature struct from user input
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_initialize_extrinsic_curvature
*/

void ghl_initialize_extrinsic_curvature(
      const double Kxx,
      const double Kxy,
      const double Kxz,
      const double Kyy,
      const double Kyz,
      const double Kzz,
      ghl_extrinsic_curvature *restrict curv) {

  curv->K[0][0]                 = Kxx;
  curv->K[1][1]                 = Kyy;
  curv->K[2][2]                 = Kzz;
  curv->K[0][1] = curv->K[1][0] = Kxy;
  curv->K[0][2] = curv->K[2][0] = Kxz;
  curv->K[1][2] = curv->K[2][1] = Kyz;
}
