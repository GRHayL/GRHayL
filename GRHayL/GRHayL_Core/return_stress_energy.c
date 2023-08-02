#include "ghl.h"

/*
 * Function     : ghl_return_stress_energy()
 * Description  : Unpacks ghl_stress_energy struct data into provided
 *                provided memory locations
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_return_stress_energy
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
