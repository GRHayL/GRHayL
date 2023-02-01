#include "GRHayL.h"

/* Function    : return_primitives()
 * Authors     : Samuel Cupp
 * Description : unpacks primitives struct into variables
 *
 * Inputs      : cons            - primitive_quantities struct with
 *                                 data from the Con2Prim routine
 *
 * Outputs     : rho            - pointer to the baryonic density
 *             : press          - pointer to the pressure
 *             : epsilon        - pointer to epsilon
 *             : vx             - pointer to 3-velocity component u^1/u^0
 *             : vy             - pointer to 3-velocity component u^2/u^0
 *             : vz             - pointer to 3-velocity component u^3/u^0
 *             : Bx             - pointer to the x component of the magnetic field TODO: which magnetic field?
 *             : By             - pointer to the y component of the magnetic field
 *             : Bz             - pointer to the z component of the magnetic field
 *             : entropy        - pointer to the entropy
 *             : Y_e            - pointer to the electron fraction (for tabulated EOS)
 *             : temp           - pointer to the temperature (for tabulated EOS)
 */

void return_primitives(const primitive_quantities *restrict prims,
                      double *restrict rho, double *restrict press, double *restrict epsilon,
                      double *restrict vx, double *restrict vy, double *restrict vz,
                      double *restrict Bx, double *restrict By, double *restrict Bz,
                      double *restrict entropy, double *restrict Y_e, double *restrict temperature) {

  *rho         = prims->rho;
  *press       = prims->press;
  *vx          = prims->vx;
  *vy          = prims->vy;
  *vz          = prims->vz;
  *entropy     = prims->entropy;
  // Tabulated EOS quantities
  *Y_e         = prims->Y_e;
  *temperature = prims->temperature;
  *epsilon     = prims->eps;

}
