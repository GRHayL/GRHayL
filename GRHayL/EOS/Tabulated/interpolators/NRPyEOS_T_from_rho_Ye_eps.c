#include "ghl_nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_T_from_rho_Ye_eps(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict T) {

  // Step 1: Set EOS table keys
  const int *keys = NULL;

  // Step 2: Declare output array
  double *outvars = NULL;

  // Step 3: Perform the interpolation
  const double root_finding_precision = eos->root_finding_precision;
  const ghl_error_codes_t error = NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(eos, 0, root_finding_precision,
                                                                                rho, Y_e, eps, NRPyEOS_eps_key,
                                                                                keys, outvars, T);

  // Step 4: Check for errors
  if(error)
    return error;

  return ghl_success;
}
