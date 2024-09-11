#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_eps_and_T_from_rho_Ye_P(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double P,
      double *restrict eps,
      double *restrict T) {

  // Step 1: Set EOS table keys
  const int keys[1] = {NRPyEOS_eps_key};

  // Step 2: Declare output array
  double outvars[1];

  // Step 3: Perform the interpolation
  const double root_finding_precision = eos->root_finding_precision;
  const ghl_error_codes_t error = NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(eos, 1, root_finding_precision,
                                                                                rho, Y_e, P, NRPyEOS_press_key,
                                                                                keys, outvars, T);

  // Step 4: Check for errors
  if(error)
    return error;

  // Step 5: Update output variables
  *eps = outvars[0];
  return ghl_success;
}
