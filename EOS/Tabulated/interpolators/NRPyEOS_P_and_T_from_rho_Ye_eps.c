#include "../EOS_tabulated.h"
#include "../NRPy_function_prototypes.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_P_and_T_from_rho_Ye_eps(const eos_parameters *restrict eos_params,
                                     const double rho,
                                     const double Y_e,
                                     const double eps,
                                     double *restrict P,
                                     double *restrict T) {

  // Step 1: Set EOS table keys
  const int keys[1] = {NRPyEOS_press_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[1];

  // Step 4: Perform the interpolation
  const double root_finding_precision = 1e-10;
  NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities( eos_params, 1,root_finding_precision,
                                                               rho,Y_e,eps,NRPyEOS_eps_key, keys,outvars, T, &report );

  // Step 5: Check for errors
  if( report.error ) {
    fprintf(stderr,"(NRPyEOS) Inside NRPyEOS_P_and_T_from_rho_Ye_eps. Error message: %s (key = %d)",report.message,report.error_key);
  }

  // Step 6: Update output variables
  *P = outvars[0];
}
