#include "grhayl.h"

/* Function    : initialize_extrinsic_curvature()
 * Description : Initialize the extrinsic_curvature struct from user input
 *
 * Inputs      : Kij            - value of the (i,j) component of the
 *                                extrinsic curvature tensor K^ij
 *
 * Outputs     : curv           - returns extrinsic_curvature struct containing
 *                                the inputs
 */

void initialize_extrinsic_curvature(
      const double Kxx, const double Kxy, const double Kxz,
      const double Kyy, const double Kyz, const double Kzz,
      extrinsic_curvature *restrict curv){

  curv->Kxx = Kxx;
  curv->Kxy = Kxy;
  curv->Kxz = Kxz;
  curv->Kyy = Kyy;
  curv->Kyz = Kyz;
  curv->Kzz = Kzz;
}
