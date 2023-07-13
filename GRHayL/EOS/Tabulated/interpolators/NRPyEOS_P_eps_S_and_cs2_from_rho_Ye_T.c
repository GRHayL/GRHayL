#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T(const ghl_eos_parameters *restrict eos,
                                           const double rho,
                                           const double Y_e,
                                           const double T,
                                           double *restrict P,
                                           double *restrict eps,
                                           double *restrict S,
                                           double *restrict cs2) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Set EOS table keys
  const int keys[4] = {NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_entropy_key,NRPyEOS_cs2_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[4];

  // Step 4: Perform the interpolation
  const int error = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos, 4, rho, Y_e, T,
                                                                   keys, outvars, &report);

  // Step 5: Check for errors
  if( error )
    ghl_error(report.message, error);

  // Step 6: Update output variables
  *P = outvars[0];
  *eps = outvars[1];
  *S = outvars[2];
  *cs2 = outvars[3];
#endif
}
