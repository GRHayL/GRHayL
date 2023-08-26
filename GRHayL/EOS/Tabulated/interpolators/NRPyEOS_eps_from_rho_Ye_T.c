#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_eps_from_rho_Ye_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict eps) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Set EOS table keys
  const int keys[1] = {NRPyEOS_eps_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[1];

  // Step 4: Perform the interpolation
  const int error = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos, 1, rho, Y_e, T,
                                                                   keys, outvars, &report);

  // Step 5: Check for errors
  if(error)
    ghl_error(report.message, error);

  // Step 6: Update output variables
  *eps = outvars[0];
#endif
}
