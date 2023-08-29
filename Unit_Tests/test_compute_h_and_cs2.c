#include "unit_tests.h"

/*
   Simple function for h and cs2 that constructs
   a simple value for cs2 without using the table.
   It does assume the table bounds of LS220, and
   might fail for others. It also assumes prims
   is full and we can directly compute h.
*/
void ghl_test_compute_h_and_cs2(
      const ghl_eos_parameters *restrict eos,
      ghl_primitive_quantities *restrict prims,
      double *restrict h,
      double *restrict cs2) {

  const double rho = prims->rho;
  const double Y_e = prims->Y_e;
  const double T   = prims->temperature;
  const double P   = prims->press;
  const double eps = prims->eps;

  const double rhomax = 1.6e-02;
  const double Yemax  = 5.5e-01;
  const double Tmax   = 2.5e+02;

  *h = 1 + eps + P/rho;
  *cs2 = (rho + T + Y_e)/(rhomax + Yemax + Tmax);
}
