#include "NRPyEOS_Tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_P_from_rho_Ye_T(const eos_parameters *restrict eos_params,
                             const double rho,
                             const double Y_e,
                             const double T,
                             double *restrict P) {

  // Step 1: Set EOS table keys
  const int keys[1] = {NRPyEOS_press_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[1];

  // Step 4: Perform the interpolation
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, 1,rho,Y_e,T, keys,outvars, &report );

  // Step 5: Check for errors
  if( report.error )
    grhayl_Error(report.error_key, report.message);

  // Step 6: Update output variables
  *P = outvars[0];
}
