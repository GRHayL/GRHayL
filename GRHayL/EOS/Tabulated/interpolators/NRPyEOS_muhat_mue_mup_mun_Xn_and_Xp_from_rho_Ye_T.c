#include "NRPyEOS_Tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T(const eos_parameters *restrict eos_params,
                                                       const double rho,
                                                       const double Y_e,
                                                       const double T,
                                                       double *restrict muhat,
                                                       double *restrict mu_e,
                                                       double *restrict mu_p,
                                                       double *restrict mu_n,
                                                       double *restrict X_n,
                                                       double *restrict X_p) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Set EOS table keys
  const int keys[6] = {NRPyEOS_muhat_key,NRPyEOS_mu_e_key,NRPyEOS_mu_p_key,NRPyEOS_mu_n_key,NRPyEOS_X_n_key,NRPyEOS_X_p_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[6];

  // Step 4: Perform the interpolation
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, 6,rho,Y_e,T, keys,outvars, &report );

  // Step 5: Check for errors
  if( report.error )
    grhayl_Error(report.error_key, report.message, report.error_key);

  // Step 6: Update output variables
  *muhat = outvars[0];
  *mu_e = outvars[1];
  *mu_p = outvars[2];
  *mu_n = outvars[3];
  *X_n = outvars[4];
  *X_p = outvars[5];
#endif
}
