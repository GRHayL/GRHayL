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
  const double rho = prims->rho;
  const double Y_e = prims->Y_e;
  const double P   = prims->press;

  // Step 2: Get P, eps, cs2
  double eps, cs2, T=eos->T_min;
  ghl_tabulated_compute_eps_cs2_T_from_P(eos, rho, Y_e, P, &eps, &cs2, &T);

  // Step 3: Compute the enthalpy
  const double h = 1 + eps + P/rho;
//printf(
//"rho %.16e\n"
//"Y_e %.16e\n"
//"P   %.16e\n"
//"eps %.16e\n"
//"cs2 %.16e\n"
//"T   %.16e\n"
//"h   %.16e\n",
//rho, Y_e, P, eps, cs2, T, h);

  // Step 4: Set the output
  *enthalpy_ptr = h;
  *cs2_ptr      = cs2;
#endif
}
