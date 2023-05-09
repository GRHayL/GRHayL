#include "grhayl.h"

/* Function    : initialize_primitives()
 * Description : Initialize the primitives struct from user input
 *
 * Inputs      : rho            - value of the baryonic density
 *             : press          - value of the pressure
 *             : epsilon        - value of epsilon
 *             : vx             - value of 3-velocity u^1/u^0
 *             : vy             - value of 3-velocity u^2/u^0
 *             : vz             - value of 3-velocity u^3/u^0
 *             : Bx             - value of magnetic field B^1 TODO: define choice of B
 *             : By             - value of magnetic field B^2
 *             : Bz             - value of magnetic field B^3
 *             : entropy        - value of entropy
 *             : Y_e            - value of electron fraction (for tabulated EOS)
 *             : temperature    - value of temperature (for tabulated EOS)
 *
 * Outputs     : prims          - returns initialized primitive_quantities
 *                                struct containing input
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
