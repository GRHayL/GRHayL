#include "nrpyeos_tabulated.h"

void NRPyEOS_tabulated_compute_enthalpy_and_cs2(
      ghl_eos_parameters const *restrict eos,
      ghl_primitive_quantities const *restrict prims,
      double *restrict enthalpy_ptr,
      double *restrict cs2_ptr ) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Unpack primitives struct
  double const rho = prims->rho;
  double const Y_e = prims->Y_e;
  double const P   = prims->press;

  // Step 2: Get P, eps, cs2
  double eps, cs2, T=eos->T_atm;
  ghl_tabulated_compute_eps_cs2_T_from_P(eos, rho, Y_e, P, &eps, &cs2, &T);

  // Step 3: Compute the enthalpy
  double const h = 1 + eps + P/rho;

  // Step 4: Set the output
  *enthalpy_ptr = h;
  *cs2_ptr      = cs2;
#endif
}
