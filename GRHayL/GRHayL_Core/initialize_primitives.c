#include "ghl.h"

/*
 * Function     : ghl_initialize_primitives()
 * Description  : Initialize the ghl_primitive_quantities struct from user input
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_initialize_primitives
*/

void ghl_initialize_primitives(
      const double rho, const double press, const double epsilon,
      const double vx, const double vy, const double vz,
      const double Bx, const double By, const double Bz,
      const double entropy, const double Y_e, const double temperature,
      ghl_primitive_quantities *restrict prims) {

  prims->rho         = rho;
  prims->press       = press;
  prims->vU[0]       = vx;
  prims->vU[1]       = vy;
  prims->vU[2]       = vz;
  prims->BU[0]       = Bx;
  prims->BU[1]       = By;
  prims->BU[2]       = Bz;
  prims->eps         = epsilon;
  prims->entropy     = entropy;
  prims->Y_e         = Y_e;
  prims->temperature = temperature;
}
