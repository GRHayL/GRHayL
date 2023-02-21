#include "NRPyEOS_Tabulated.h"

void NRPyEOS_tabulated_compute_enthalpy_and_cs2(
      eos_parameters const *restrict eos,
      primitive_quantities const *restrict prims,
      double *restrict enthalpy_ptr,
      double *restrict cs2_ptr ) {

  // Step 1: Unpack primitives struct
  double const rho = prims->rho;
  double const Y_e = prims->Y_e;
  double const T   = prims->temperature;

  // Step 2: Get P, eps, cs2
  double P, eps, cs2;
  eos->tabulated_compute_P_eps_cs2_from_T(eos, rho, Y_e, T, &P, &eps, &cs2);

  // Step 3: Compute the enthalpy
  double const h = 1 + eps + P/rho;

  // Step 4: Set the output
  *enthalpy_ptr = h;
  *cs2_ptr      = cs2;
}
