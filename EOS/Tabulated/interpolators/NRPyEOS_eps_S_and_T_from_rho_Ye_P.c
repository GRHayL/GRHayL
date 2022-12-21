#include "EOS_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_eps_S_and_T_from_rho_Ye_P(const eos_parameters *restrict eos_params,
                                       const double rho,
                                       const double Y_e,
                                       const double P,
                                       double *restrict eps,
                                       double *restrict S,
                                       double *restrict T) {

  // Step 1: Set EOS table keys
  const int keys[2] = {NRPyEOS_eps_key,NRPyEOS_entropy_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[2];

  // Step 4: Perform the interpolation
  const double root_finding_precision = 1e-10;
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, 2,root_finding_precision,
                                                               rho,Y_e,P,NRPyEOS_press_key, keys,outvars, T, &report );

  // Step 5: Check for errors
  if( report.error ) {
    fprintf(stderr,"(GRHayL - EOS) Inside NRPyEOS_eps_S_and_T_from_rho_Ye_P. Error message: %s (key = %d)",report.message,report.error_key);
  }

  // Step 6: Update output variables
  *eps = outvars[0];
  *S = outvars[1];
}
