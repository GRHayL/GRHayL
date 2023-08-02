#include "ghl.h"

/*
 * Function      : ghl_return_primitives()
 * Description   : Unpacks ghl_primitive_quantities struct data into provided
 *                 provided memory locations
 * Documentation : https://github.com/GRHayL/GRHayL/wiki/ghl_return_primitives
*/

void ghl_return_primitives(const ghl_primitive_quantities *restrict prims,
      double *restrict rho, double *restrict press, double *restrict epsilon,
      double *restrict vx, double *restrict vy, double *restrict vz,
      double *restrict Bx, double *restrict By, double *restrict Bz,
      double *restrict entropy, double *restrict Y_e, double *restrict temperature) {

  *rho         = prims->rho;
  *press       = prims->press;
  *epsilon     = prims->eps;
  *vx          = prims->vU[0];
  *vy          = prims->vU[1];
  *vz          = prims->vU[2];
  *Bx          = prims->BU[0];
  *By          = prims->BU[1];
  *Bz          = prims->BU[2];
  *entropy     = prims->entropy;
  // Tabulated EOS quantities
  *Y_e         = prims->Y_e;
  *temperature = prims->temperature;
  *epsilon     = prims->eps;
}
