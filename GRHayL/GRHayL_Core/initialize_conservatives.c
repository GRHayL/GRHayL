#include "ghl.h"

/**
 * @ingroup pack_struct
 * @brief Initialize the conservative variable struct from user input
 *
 * @details
 * This function takes pointwise information about the hydrodynamic
 * conservative variables and uses it to initialize every element of the given
 * ghl_conservative_quantities struct. This struct is used for both densitized
 * and undensitized conservative variables, and Con2Prim functions which expect
 * the undensitized values explicitly state this in their documentation and variable
 * names. All input values use the caller's storage convention.
 *
 * @param[in] rho density conservative value
 *
 * @param[in] tau energy conservative value
 *
 * @param[in] S_x,S_y,S_z momentum conservative components
 *
 * @param[in] entropy entropy conservative value
 *
 * @param[in] Y_e electron-fraction conservative value
 *
 * @param[out] cons pointer to ghl_conservative_quantities struct
 *
 */
void ghl_initialize_conservatives(
      const double rho,
      const double tau,
      const double S_x,
      const double S_y,
      const double S_z,
      const double entropy,
      const double Y_e,
      ghl_conservative_quantities *restrict cons) {

  cons->rho = rho;
  cons->tau = tau;
  cons->SD[0] = S_x;
  cons->SD[1] = S_y;
  cons->SD[2] = S_z;
  cons->entropy = entropy;
  cons->Y_e = Y_e;
}
