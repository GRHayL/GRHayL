#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict S,
      double *restrict cs2) {

  // Step 1: Set EOS table keys
  const int keys[4] = {NRPyEOS_press_key,NRPyEOS_eps_key,NRPyEOS_entropy_key,NRPyEOS_cs2_key};

  // Step 2: Declare output array
  double outvars[4];

  // Step 3: Perform the interpolation
  const ghl_error_codes_t error = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos, 4, rho, Y_e, T,
                                                                   keys, outvars);

  // Step 4: Check for errors
  if(error)
    return error;

  // Step 5: Update output variables
  *P = outvars[0];
  *eps = outvars[1];
  *S = outvars[2];
  *cs2 = outvars[3];
  return ghl_success;
}
