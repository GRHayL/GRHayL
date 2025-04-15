#include "ghl.h"

/**
 * @ingroup pack_struct
 * @brief Unpack conservative variable struct into provided memory locations
 *
 * @details
 * This function takes data from the ghl_conservative_quantities struct
 * and unpacks (i.e. copies) the data into the memory locations passed
 * to the function.
 *
 * @param[in] cons:           pointer to ghl_conservative_quantities struct
 *
 * @param[out] rho:           pointer to density variable \f$ \rho_* \f$ (also called \f$ \tilde{D} \f$
 *
 * @param[out] tau:           pointer to energy variable \f$ \tau \f$
 *
 * @param[out] S_x, S_y, S_z: pointer to components of the momentum variable \f$ S_i \f$
 *
 * @param[out] entropy:       pointer to fluid entropy variable \f$ \tilde{S} \f$
 *
 * @param[out] Y_e:           pointer to electron fraction variable \f$ \tilde{Y}_e \f$
 *
 * @returns void
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
