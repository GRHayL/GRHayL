#include "ghl_nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_P_and_T_from_rho_Ye_S(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double S,
      double *restrict P,
      double *restrict T) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
  return ghl_error_hdf5_is_disabled;
#else

  // Step 1: Set EOS table keys
  const int keys[1] = {NRPyEOS_press_key};

  // Step 2: Declare output array
  double outvars[1];

  // Step 3: Perform the interpolation
  const double root_finding_precision = eos->root_finding_precision;
  const ghl_error_codes_t error = NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(eos, 1, root_finding_precision,
                                                                                rho, Y_e, S, NRPyEOS_entropy_key,
                                                                                keys, outvars, T);

  // Step 4: Check for errors
  if(error)
    return error;

  // Step 5: Update output variables
  *P = outvars[0];
  return ghl_success;
#endif
}
