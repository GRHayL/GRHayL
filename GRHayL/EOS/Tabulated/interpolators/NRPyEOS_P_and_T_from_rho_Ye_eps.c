#include "NRPyEOS_Tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_P_and_T_from_rho_Ye_eps(const eos_parameters *restrict eos_params,
                                     const double rho,
                                     const double Y_e,
                                     const double eps,
                                     double *restrict P,
                                     double *restrict T) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Set EOS table keys
  const int keys[1] = {NRPyEOS_press_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[1];

  // Step 4: Perform the interpolation
  const double root_finding_precision = eos_params->root_finding_precision;
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, 1,root_finding_precision,
                                                               rho,Y_e,eps,NRPyEOS_eps_key, keys,outvars, T, &report );

  // Step 5: Check for errors
  if( report.error )
    grhayl_Error(report.error_key, report.message, report.error_key);

  // Step 6: Update output variables
  *P = outvars[0];
#endif
}
