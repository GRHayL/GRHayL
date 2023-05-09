#include "grhayl.h"

/* Function    : initialize_stress_energy()
 * Description : Initialize the stress_energy struct from user input
 *
 * Inputs      : Tij            - value of the (i,j) component of the
 *                                stress energy tensor T^\mu\nu
 *
 * Outputs     : Tmunu          - returns stress_energy struct containing
 *                                the inputs
 */

void initialize_stress_energy(
      const double Ttt,
      const double Ttx, const double Tty, const double Ttz,
      const double Txx, const double Txy, const double Txz,
      const double Tyy, const double Tyz, const double Tzz,
      stress_energy *restrict Tmunu) {

  Tmunu->Ttt = Ttt;
  Tmunu->Ttx = Ttx;
  Tmunu->Tty = Tty;
  Tmunu->Ttz = Ttz;
  Tmunu->Txx = Txx;
  Tmunu->Txy = Txy;
  Tmunu->Txz = Txz;
  Tmunu->Tyy = Tyy;
  Tmunu->Tyz = Tyz;
  Tmunu->Tzz = Tzz;
}
