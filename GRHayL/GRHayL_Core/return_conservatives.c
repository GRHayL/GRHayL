#include "ghl.h"

/**
 * @ingroup unpack_struct
 * @brief Unpack conservative variable struct into provided memory locations
 *
 * @details
 * This function takes data from the ghl_conservative_quantities struct
 * and unpacks (i.e. copies) the data into the memory locations passed
 * to the function. Output values use the same storage convention as the
 * input struct.
 *
 * @param[in] cons pointer to ghl_conservative_quantities struct
 *
 * @param[out] rho pointer to density conservative value
 *
 * @param[out] tau pointer to energy conservative value
 *
 * @param[out] S_x pointer to the x-component of the momentum conservative value
 *
 * @param[out] S_y pointer to the y-component of the momentum conservative value
 *
 * @param[out] S_z pointer to the z-component of the momentum conservative value
 *
 * @param[out] entropy pointer to entropy conservative value
 *
 * @param[out] Y_e pointer to electron-fraction conservative value
 *
 */
void ghl_return_conservatives(
      const ghl_conservative_quantities *restrict cons,
      double *restrict rho,
      double *restrict tau,
      double *restrict S_x,
      double *restrict S_y,
      double *restrict S_z,
      double *restrict entropy,
      double *restrict Y_e) {

  *rho = cons->rho;
  *S_x = cons->SD[0];
  *S_y = cons->SD[1];
  *S_z = cons->SD[2];
  *tau = cons->tau;
  *entropy = cons->entropy;
  *Y_e = cons->Y_e;
}
