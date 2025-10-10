#include "ghl_nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_P_eps_and_S_from_rho_Ye_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict S) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
  return ghl_error_hdf5_is_disabled;
#else

  // Step 1: Set EOS table keys
  const int keys[3] = {NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_entropy_key};

  // Step 2: Declare output array
  double outvars[3];

  // Step 3: Perform the interpolation
  const ghl_error_codes_t error = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos, 3, rho, Y_e, T,
                                                                   keys, outvars);

  // Step 4: Check for errors
  if(error)
    return error;

  // Step 5: Update output variables
  *P = outvars[0];
  *eps = outvars[1];
  *S = outvars[2];
  return ghl_success;
#endif
}
