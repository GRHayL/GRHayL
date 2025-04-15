#include "ghl.h"

/**
 * @ingroup pack_struct
 * @brief Unpacks primitive variable struct into provided provided memory locations
 *
 * @details
 * This function takes data from the ghl_primitive_quantities struct and
 * unpacks (i.e. copies) the data into the memory locations passed to the
 * function.
 *
 * @param[in] prims:        pointer to ghl_primitive_quantities struct
 *
 * @param[out] rho:         pointer to baryonic density \f$ \rho \f$
 *
 * @param[out] press:       pointer to fluid pressure \f$ P \f$
 *
 * @param[out] epsilon:     pointer to TODO: \f$ \epsilon \f$
 *
 * @param[out] vx, vy, vz:  pointer to components of the fluid velocity \f$ v^i \f$
 *
 * @param[out] Bx, By, Bz:  pointer to components of the magnetic field \f$ B^i \f$
 *
 * @param[out] entropy:     pointer to fluid entropy
 *
 * @param[out] Y_e:         pointer to fluid electron fraction
 *
 * @param[out] temperature: pointer to fluid temperature
 *
 * @returns void
 */
void ghl_return_primitives(const ghl_primitive_quantities *restrict prims,
      double *restrict rho, double *restrict press, double *restrict epsilon,
      double *restrict vx, double *restrict vy, double *restrict vz,
      double *restrict Bx, double *restrict By, double *restrict Bz,
      double *restrict entropy, double *restrict Y_e, double *restrict temperature) {

  *rho         = prims->rho;
  *press       = prims->press;
  *epsilon     = prims->eps;
  *vx          = prims->vU[0];
  *vy          = prims->vU[1];
  *vz          = prims->vU[2];
  *Bx          = prims->BU[0];
  *By          = prims->BU[1];
  *Bz          = prims->BU[2];
  *entropy     = prims->entropy;
  // Tabulated EOS quantities
  *Y_e         = prims->Y_e;
  *temperature = prims->temperature;
  *epsilon     = prims->eps;
}
