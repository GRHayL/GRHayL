#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_P_S_depsdT_and_T_from_rho_Ye_eps(const eos_parameters *restrict eos,
                                              const double rho,
                                              const double Y_e,
                                              const double eps,
                                              double *restrict P,
                                              double *restrict S,
                                              double *restrict depsdT,
                                              double *restrict T) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Set EOS table keys
  const int keys[3] = {NRPyEOS_press_key,NRPyEOS_entropy_key,NRPyEOS_depsdT_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[3];

  // Step 4: Perform the interpolation
  const double root_finding_precision = eos->root_finding_precision;
  const int error = NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(eos, 3, root_finding_precision,
                                                                                rho, Y_e, eps, NRPyEOS_eps_key,
                                                                                keys, outvars, T, &report);

  // Step 5: Check for errors
  if( error )
    grhayl_error(report.message, error);

  // Step 6: Update output variables
  *P = outvars[0];
  *S = outvars[1];
  *depsdT = outvars[2];
#endif
}
