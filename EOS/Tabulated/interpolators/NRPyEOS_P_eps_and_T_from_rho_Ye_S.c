#include "../NRPy_basic_defines.h"
#include "../NRPy_function_prototypes.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_P_eps_and_T_from_rho_Ye_S(const NRPyEOS_params *restrict eos_params,
                                       const double rho,
                                       const double Y_e,
                                       const double S,
                                       double *restrict P,
                                       double *restrict eps,
                                       double *restrict T) {

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
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_eps_and_T_from_rho_Ye_S. Error message: %s (key = %d)",report.message,report.error_key);
  }

  // Step 6: Update output variables
  *P = outvars[0];
  *eps = outvars[1];
}
