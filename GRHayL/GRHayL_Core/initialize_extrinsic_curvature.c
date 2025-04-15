#include "ghl.h"

/**
 * @ingroup pack_struct
 * @brief Initialize the extrinsic curvature struct from user input
 *
 * @details
 * This function takes pointwise information about the extrinsic curvature
 * and uses it to initialize every element of the given ghl_extrinsic_curvature
 * struct.
 *
 * @param[in] Kxx, Kxy, Kxz, Kyy, Kyz, Kzz: individual components of the extrinsic curvature \f$ K_{ij} \f$
 *
 * @param[out] curv: pointer to ghl_conservative_quantities struct
 *
 * @returns void
 */
void ghl_initialize_extrinsic_curvature(
      const double Kxx,
      const double Kxy,
      const double Kxz,
      const double Kyy,
      const double Kyz,
      const double Kzz,
      ghl_extrinsic_curvature *restrict curv) {

  curv->K[0][0]                 = Kxx;
  curv->K[1][1]                 = Kyy;
  curv->K[2][2]                 = Kzz;
  curv->K[0][1] = curv->K[1][0] = Kxy;
  curv->K[0][2] = curv->K[2][0] = Kxz;
  curv->K[1][2] = curv->K[2][1] = Kyz;
}
