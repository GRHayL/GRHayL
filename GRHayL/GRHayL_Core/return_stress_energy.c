#include "grhayl.h"

/* Function    : return_stress_energy()
 * Description : unpacks stress_energy struct into variables
 *
 * Inputs      : Tmunu          - stress_energy struct to be unpacked
 *
 * Outputs     : Tij            - pointer to variable for (i,j) component
 *                                of T^\mu\nu
 */

void return_stress_energy(
      const stress_energy *restrict Tmunu,
      double *restrict Ttt, double *restrict Ttx, double *restrict Tty,
      double *restrict Ttz, double *restrict Txx, double *restrict Txy,
      double *restrict Txz, double *restrict Tyy, double *restrict Tyz,
      double *restrict Tzz) {

  *Ttt = Tmunu->Ttt;
  *Ttx = Tmunu->Ttx;
  *Tty = Tmunu->Tty;
  *Ttz = Tmunu->Ttz;
  *Txx = Tmunu->Txx;
  *Txy = Tmunu->Txy;
  *Txz = Tmunu->Txz;
  *Tyy = Tmunu->Tyy;
  *Tyz = Tmunu->Tyz;
  *Tzz = Tmunu->Tzz;
}
