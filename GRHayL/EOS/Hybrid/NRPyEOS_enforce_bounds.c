#include "ghl_nrpyeos_hybrid.h"

/**
 * @ingroup eos_internal
 * @brief Ensure that the input value of rho is within the set bounds
 *
 * @param[in] eos: pointer to ghl_eos_parameters struct
 *
 * @param[in] rho: density value
 *
 * @returns whether density was in the bounds (T) or not (F)
 */
bool NRPyEOS_hybrid_enforce_bounds__rho(
      const ghl_eos_parameters *restrict eos,
      double *restrict rho) {

  bool in_bounds = true;
  if (*rho < eos->rho_min) {
    *rho = eos->rho_min;
    in_bounds = false;
  } else if (*rho > eos->rho_max) {
    in_bounds = false;
    *rho = eos->rho_max;
  }

  return in_bounds;
}
