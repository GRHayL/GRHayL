#include "NRPyEOS_Tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_P_eps_and_cs2_from_rho_Ye_T(const eos_parameters *restrict eos_params,
                                         const double rho,
                                         const double Y_e,
                                         const double T,
                                         double *restrict P,
                                         double *restrict eps,
                                         double *restrict cs2) {

  // Step 1: Set EOS table keys
  const int keys[3] = {NRPyEOS_press_key, NRPyEOS_eps_key, NRPyEOS_cs2_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[3];

  // Step 4: Perform the interpolation
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities( eos_params, 3, rho, Y_e, T, keys,outvars, &report );

  // Step 5: Check for errors
  if( report.error )
    grhayl_Error(report.error_key, report.message);

  // Step 6: Update output variables
  *P   = outvars[0];
  *eps = outvars[1];
  *cs2 = outvars[2];
}
