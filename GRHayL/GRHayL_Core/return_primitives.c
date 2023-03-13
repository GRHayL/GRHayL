#include "GRHayL.h"

/* Function    : return_primitives()
 * Description : unpacks primitives struct into variables
 *
 * Inputs      : cons            - primitive_quantities struct to be unpacked
 *
 * Outputs     : rho            - pointer to baryonic density
 *             : press          - pointer to pressure
 *             : epsilon        - pointer to epsilon
 *             : vx             - pointer to 3-velocity component u^1/u^0
 *             : vy             - pointer to 3-velocity component u^2/u^0
 *             : vz             - pointer to 3-velocity component u^3/u^0
 *             : Bx             - pointer to magnetic field B^1 TODO: define choice of B
 *             : By             - pointer to magnetic field B^2
 *             : Bz             - pointer to magnetic field B^3
 *             : entropy        - pointer to entropy
 *             : Y_e            - pointer to electron fraction (for tabulated EOS)
 *             : temp           - pointer to temperature (for tabulated EOS)
 */

void return_primitives(const primitive_quantities *restrict prims,
                      double *restrict rho, double *restrict press, double *restrict epsilon,
                      double *restrict vx, double *restrict vy, double *restrict vz,
                      double *restrict Bx, double *restrict By, double *restrict Bz,
                      double *restrict entropy, double *restrict Y_e, double *restrict temperature) {

  *rho         = prims->rho;
  *press       = prims->press;
  *epsilon     = prims->eps;
  *vx          = prims->vx;
  *vy          = prims->vy;
  *vz          = prims->vz;
  *Bx          = prims->Bx;
  *By          = prims->By;
  *Bz          = prims->Bz;
  *entropy     = prims->entropy;
  // Tabulated EOS quantities
  *Y_e         = prims->Y_e;
  *temperature = prims->temperature;
  *epsilon     = prims->eps;

}
