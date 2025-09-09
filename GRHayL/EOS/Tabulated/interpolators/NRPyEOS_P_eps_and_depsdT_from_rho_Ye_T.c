#include "ghl_nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_P_eps_and_depsdT_from_rho_Ye_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict depsdT) {

  // Step 1: Set EOS table keys
  const int keys[3] = {NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_depsdT_key};

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
  *depsdT = outvars[2];
  return ghl_success;
}
