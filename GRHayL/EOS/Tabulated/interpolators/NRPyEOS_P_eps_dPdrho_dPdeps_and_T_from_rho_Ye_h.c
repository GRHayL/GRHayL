#include "ghl_nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_P_eps_dPdrho_dPdeps_and_T_from_rho_Ye_h(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double h,
      double *restrict P,
      double *restrict eps,
      double *restrict dPdrho,
      double *restrict dPdeps,
      double *restrict T) {

  // Step 1: Set EOS table keys
  const int keys[4] = {
    NRPyEOS_press_key,
    NRPyEOS_eps_key,
    NRPyEOS_dPdrho_key,
    NRPyEOS_dPdeps_key,
  };

  // Step 2: Declare output array
  double outvars[4] = {NAN, NAN, NAN, NAN};

  // Step 3: Perform the interpolation
  const double root_finding_precision = eos->root_finding_precision;
  const ghl_error_codes_t error = NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(
                                      eos, 4, root_finding_precision, rho, Y_e, h,
                                      NRPyEOS_enthalpy_key, keys, outvars, T);

  // Step 4: Check for errors
  if(error)
    return error;

  // Step 5: Update output variables
  *P      = outvars[0];
  *eps    = outvars[1];
  *dPdrho = outvars[2];
  *dPdeps = outvars[3];
  return ghl_success;
}
