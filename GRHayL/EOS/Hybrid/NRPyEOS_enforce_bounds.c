#include "ghl_nrpyeos_hybrid.h"

/*
 * Function     : NRPyEOS_enforce_bounds__rho()
 * Description  : Ensure that the input value of rho is within the set bounds
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_hybrid_enforce_bounds__rho
*/
GHL_DEVICE
bool NRPyEOS_hybrid_enforce_bounds__rho(const ghl_eos_parameters *restrict eos, double *restrict rho) {
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
