#include "ghl.h"

/* Function    : ghl_initialize_stress_energy()
 * Description : Initialize the ghl_stress_energy struct from user input
 *
 * Inputs      : Tij            - value of the (i,j) component of the
 *                                stress energy tensor T^\mu\nu
 *
 * Outputs     : Tmunu          - returns ghl_stress_energy struct containing
 *                                the inputs
 */

void ghl_initialize_stress_energy(
      const double Ttt,
      const double Ttx, const double Tty, const double Ttz,
      const double Txx, const double Txy, const double Txz,
      const double Tyy, const double Tyz, const double Tzz,
      ghl_stress_energy *restrict Tmunu) {

  Tmunu->T4[0][0]                   = Ttt;
  Tmunu->T4[1][1]                   = Txx;
  Tmunu->T4[2][2]                   = Tyy;
  Tmunu->T4[3][3]                   = Tzz;
  Tmunu->T4[0][1] = Tmunu->T4[1][0] = Ttx;
  Tmunu->T4[0][2] = Tmunu->T4[2][0] = Tty;
  Tmunu->T4[0][3] = Tmunu->T4[3][0] = Ttz;
  Tmunu->T4[1][2] = Tmunu->T4[2][1] = Txy;
  Tmunu->T4[1][3] = Tmunu->T4[3][1] = Txz;
  Tmunu->T4[2][3] = Tmunu->T4[3][2] = Tyz;
}
