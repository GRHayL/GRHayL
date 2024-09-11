#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_P_from_rho_Ye_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P) {

  // Step 1: Set EOS table keys
  const int keys[1] = {NRPyEOS_press_key};

  // Step 2: Declare output array
  double outvars[1];

  // Step 3: Perform the interpolation
  const ghl_error_codes_t error = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos, 1, rho, Y_e, T,
                                                                   keys, outvars);

  // Step 4: Check for errors
  if(error)
    return error;

  // Step 5: Update output variables
  *P = outvars[0];
  return ghl_success;
}
