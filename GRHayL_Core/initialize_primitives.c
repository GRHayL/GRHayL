#include "GRHayL.h"

/* Function    : initialize_primitives()
 * Description : Initialize the primitives struct from user
 *               input
 *
 * Inputs      : rho            - value of the baryonic density
 *             : press          - value of the pressure
 *             : epsilon        - value of epsilon
 *             : vx             - value of the x component of the 3-velocity u^i/u^0
 *             : vy             - value of the y component of the 3-velocity u^i/u^0
 *             : vz             - value of the z component of the 3-velocity u^i/u^0
 *             : Bx             - value of the x component of the magnetic field TODO: which magnetic field?
 *             : By             - value of the y component of the magnetic field
 *             : Bz             - value of the z component of the magnetic field
 *             : entropy        - value of the entropy
 *             : Y_e            - value of the electron fraction (for tabulated EOS)
 *             : temperature    - value of the temperature (for tabulated EOS)
 *
 * Outputs     : prims          - fully initialized primitive_quantities
 *                                struct containing the input data
 */

void initialize_primitives(
             const double rho, const double press, const double epsilon,
             const double vx, const double vy, const double vz,
             const double Bx, const double By, const double Bz,
             const double entropy, const double Y_e, const double temperature,
             primitive_quantities *restrict prims) {
  prims->rho         = rho;
  prims->press       = press;
  prims->vx          = vx;
  prims->vy          = vy;
  prims->vz          = vz;
  prims->Bx          = Bx;
  prims->By          = By;
  prims->Bz          = Bz;
  prims->eps         = epsilon;
  prims->entropy     = entropy;
  prims->Y_e         = Y_e;
  prims->temperature = temperature;
}
