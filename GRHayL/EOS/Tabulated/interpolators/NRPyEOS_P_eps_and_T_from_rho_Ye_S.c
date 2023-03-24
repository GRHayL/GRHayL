#include "NRPyEOS_Tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_P_eps_and_T_from_rho_Ye_S(const eos_parameters *restrict eos_params,
                                       const double rho,
                                       const double Y_e,
                                       const double S,
                                       double *restrict P,
                                       double *restrict eps,
                                       double *restrict T) {
#ifndef USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Set EOS table keys
  const int keys[2] = {NRPyEOS_press_key,NRPyEOS_eps_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[2];

  // Step 4: Perform the interpolation
  const double root_finding_precision = 1e-10;
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, 2,root_finding_precision,
                                                               rho,Y_e,S,NRPyEOS_entropy_key, keys,outvars, T, &report );

  // Step 5: Check for errors
  if( report.error )
    grhayl_Error(report.error_key, report.message);

  // Step 6: Update output variables
  *P = outvars[0];
  *eps = outvars[1];
#endif
}
