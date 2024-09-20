#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p) {

  // Step 1: Set EOS table keys
  const int keys[6] = {NRPyEOS_muhat_key,NRPyEOS_mu_e_key,NRPyEOS_mu_p_key,NRPyEOS_mu_n_key,NRPyEOS_X_n_key,NRPyEOS_X_p_key};

  // Step 2: Declare output array
  double outvars[6];

  // Step 3: Perform the interpolation
  const ghl_error_codes_t error = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos, 6, rho, Y_e, T,
                                                                   keys, outvars);

  // Step 4: Check for errors
  if(error)
    return error;

  // Step 5: Update output variables
  *muhat = outvars[0];
  *mu_e = outvars[1];
  *mu_p = outvars[2];
  *mu_n = outvars[3];
  *X_n = outvars[4];
  *X_p = outvars[5];

  return ghl_success;
}
