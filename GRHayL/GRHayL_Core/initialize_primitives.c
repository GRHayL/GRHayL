#include "ghl.h"

/**
 * @ingroup pack_struct
 * @brief Initialize the primitive variable struct from user input
 *
 * @details
 * This function takes pointwise information about the hydrodynamic
 * primitive variables and uses it to initialize the elements of the
 * given ghl_primitive_quantities struct. Note that it does not initialize
 * ghl_primitive_quantities::u0.
 *
 * @param[in] rho:         baryonic density \f$ \rho \f$
 *
 * @param[in] press:       fluid pressure \f$ P \f$
 *
 * @param[in] epsilon:     TODO: \f$ \epsilon \f$
 *
 * @param[in] vx, vy, vz:  components of the fluid velocity \f$ v^i \f$
 *
 * @param[in] Bx, By, Bz:  components of the magnetic field \f$ B^i \f$
 *
 * @param[in] entropy:     fluid entropy
 *
 * @param[in] Y_e:         fluid electron fraction
 *
 * @param[in] temperature: fluid temperature
 *
 * @param[out] prims:      pointer to ghl_primitive_quantities struct
 *
 * @returns void
 */
void ghl_initialize_primitives(
      const double rho, const double press, const double epsilon,
      const double vx, const double vy, const double vz,
      const double Bx, const double By, const double Bz,
      const double entropy, const double Y_e, const double temperature,
      ghl_primitive_quantities *restrict prims) {

  prims->rho         = rho;
  prims->press       = press;
  prims->vU[0]       = vx;
  prims->vU[1]       = vy;
  prims->vU[2]       = vz;
  prims->BU[0]       = Bx;
  prims->BU[1]       = By;
  prims->BU[2]       = Bz;
  prims->eps         = epsilon;
  prims->entropy     = entropy;
  prims->Y_e         = Y_e;
  prims->temperature = temperature;
}
