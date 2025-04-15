#include "ghl.h"

/**
 * @ingroup pack_struct
 * @brief Initialize the stress-energy tensor struct from user input
 *
 * @details
 * This function takes pointwise information about the stress-energy tensor
 * and uses it to initialize every element of the given ghl_stress_energy
 * struct. This can be used for either \f$ T^{\mu\nu} \f$ or \f$ T_{\mu\nu} \f$,
 * as neither this function nor the struct distinguish between raised or
 * lowered indices.
 *
 * @param[in] Ttt, Ttx, Tty, Ttz, Txx, Txy, Txz, Tyy, Tyz, Tzz:
 *            components of the stress-energy tensor \f$ T^{\mu\nu} \f$ or \f$ T_{\mu\nu} \f$
 *
 * @param[out] Tmunu: pointer to ghl_stress_energy struct
 *
 * @returns void
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
