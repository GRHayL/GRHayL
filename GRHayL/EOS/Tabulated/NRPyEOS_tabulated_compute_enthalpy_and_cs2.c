#include "nrpyeos_tabulated.h"

void NRPyEOS_tabulated_compute_enthalpy_and_cs2(
    const ghl_eos_parameters *restrict eos,
    ghl_primitive_quantities *restrict prims,
    double *restrict enthalpy_ptr,
    double *restrict cs2_ptr) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Unpack primitives struct
  ghl_tabulated_enforce_bounds_rho_Ye_P(eos, &prims->rho, &prims->Y_e, &prims->press);

  // Step 2: Get P, eps, cs2
  prims->temperature = 5;
  double cs2;
  ghl_tabulated_compute_eps_cs2_T_from_P(
      eos, prims->rho, prims->Y_e, prims->press, &prims->eps, &cs2, &prims->temperature);

  // Step 3: Compute the enthalpy
  const double h = 1 + prims->eps + prims->press / prims->rho;

  // Step 4: Set the output
  *enthalpy_ptr = h;
  *cs2_ptr = cs2 / h;
#endif
}
