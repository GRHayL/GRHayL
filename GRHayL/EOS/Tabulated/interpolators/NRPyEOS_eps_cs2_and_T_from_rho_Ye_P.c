#include "nrpyeos_tabulated.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_eps_cs2_and_T_from_rho_Ye_P(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double P,
      double *restrict eps,
      double *restrict cs2,
      double *restrict T) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else
  // Step 1: Set EOS table keys
  const int keys[2] = {NRPyEOS_eps_key, NRPyEOS_cs2_key};

  // Step 2: Declare EOS error report struct
  NRPyEOS_error_report report;

  // Step 3: Declare output array
  double outvars[2];

  // Step 4: Perform the interpolation
  const double root_finding_precision = eos->root_finding_precision;
  const int error = NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(eos, 2, root_finding_precision,
                                                                                rho, Y_e, P, NRPyEOS_press_key,
                                                                                keys, outvars, T, &report);

  // Step 5: Check for errors
  if(error)
    ghl_error(report.message, error);

  // Step 6: Update output variables
  *eps = outvars[0];
  *cs2 = outvars[1];
#endif
}
