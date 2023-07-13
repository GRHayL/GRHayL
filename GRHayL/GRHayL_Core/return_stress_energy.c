#include "ghl.h"

/* Function    : ghl_return_stress_energy()
 * Description : unpacks ghl_stress_energy struct into variables
 *
 * Inputs      : Tmunu          - ghl_stress_energy struct to be unpacked
 *
 * Outputs     : Tij            - pointer to variable for (i,j) component
 *                                of T^\mu\nu
 */

void ghl_return_stress_energy(
      const ghl_stress_energy *restrict Tmunu,
      double *restrict Ttt, double *restrict Ttx, double *restrict Tty,
      double *restrict Ttz, double *restrict Txx, double *restrict Txy,
      double *restrict Txz, double *restrict Tyy, double *restrict Tyz,
      double *restrict Tzz) {

  *Ttt = Tmunu->T4[0][0];
  *Ttx = Tmunu->T4[0][1];
  *Tty = Tmunu->T4[0][2];
  *Ttz = Tmunu->T4[0][3];
  *Txx = Tmunu->T4[1][1];
  *Txy = Tmunu->T4[1][2];
  *Txz = Tmunu->T4[1][3];
  *Tyy = Tmunu->T4[2][2];
  *Tyz = Tmunu->T4[2][3];
  *Tzz = Tmunu->T4[3][3];
}
