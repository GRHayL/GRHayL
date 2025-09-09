#include "ghl_nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_P_eps_and_T_from_rho_Ye_S(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double S,
      double *restrict P,
      double *restrict eps,
      double *restrict T) {

  // Step 1: Set EOS table keys
  const int keys[2] = {NRPyEOS_press_key,NRPyEOS_eps_key};

  // Step 2: Declare output array
  double outvars[2];

  // Step 2: Perform the interpolation
  const double root_finding_precision = eos->root_finding_precision;
  const ghl_error_codes_t error = NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(eos, 2, root_finding_precision,
                                                                                rho, Y_e, S, NRPyEOS_entropy_key,
                                                                                keys, outvars, T);

  // Step 3: Check for errors
  if(error)
    return error;

  // Step 4: Update output variables
  *P = outvars[0];
  *eps = outvars[1];
  return ghl_success;
}
