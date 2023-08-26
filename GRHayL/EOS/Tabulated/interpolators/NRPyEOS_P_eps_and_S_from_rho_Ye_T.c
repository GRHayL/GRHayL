#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_P_eps_and_S_from_rho_Ye_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict S) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Set EOS table keys
  const int keys[3] = {NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_entropy_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[3];

  // Step 4: Perform the interpolation
  const int error = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos, 3, rho, Y_e, T,
                                                                   keys, outvars, &report);

  // Step 5: Check for errors
  if(error)
    ghl_error(report.message, error);

  // Step 6: Update output variables
  *P = outvars[0];
  *eps = outvars[1];
  *S = outvars[2];
#endif
}
