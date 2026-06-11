#include "ghl.h"

/**
 * @ingroup unpack_struct
 * @brief Unpacks stress-energy tensor struct data into provided memory locations
 *
 * @details
 * This function takes data from the ghl_stress_energy struct and unpacks
 * (i.e. copies) the data into the memory locations passed to the function.
 *
 * @param[in] Tmunu pointer to ghl_stress_energy struct
 *
 * @param[out] Ttt pointer to the tt-component of the stress-energy tensor
 *                 \f$ T^{\mu\nu} \f$ or \f$ T_{\mu\nu} \f$
 *
 * @param[out] Ttx pointer to the tx-component of the stress-energy tensor
 *
 * @param[out] Tty pointer to the ty-component of the stress-energy tensor
 *
 * @param[out] Ttz pointer to the tz-component of the stress-energy tensor
 *
 * @param[out] Txx pointer to the xx-component of the stress-energy tensor
 *
 * @param[out] Txy pointer to the xy-component of the stress-energy tensor
 *
 * @param[out] Txz pointer to the xz-component of the stress-energy tensor
 *
 * @param[out] Tyy pointer to the yy-component of the stress-energy tensor
 *
 * @param[out] Tyz pointer to the yz-component of the stress-energy tensor
 *
 * @param[out] Tzz pointer to the zz-component of the stress-energy tensor
 *
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
