#include "nrpyeos_tabulated.h"
#include "NRPyEOS_tabulated_helpers.h"
/*
 * (c) 2022 Leo Werneck
 */
int NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(
      const ghl_eos_parameters *restrict eos,
      const int n,
      const double rho,
      const double Y_e,
      const double T,
      const int *restrict tablevars_keys,
      double *restrict tablevars,
      NRPyEOS_error_report *restrict report) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else

  if( !n ) return 0;

  // Start by assuming no errors
  int error = 0;

  // This function will interpolate n table quantities from
  // (rho,Ye,T). It replaces EOS_Omni calls with keytemp = 1
  if( n > NRPyEOS_ntablekeys ) {
    sprintf(report->message, "In %s call, number of quantities exceed maximum allowed: %d > %d.\n",
            __func__, n, NRPyEOS_ntablekeys);
    return 1;
  }

  // Check table bounds for input variables
  error = NRPyEOS_checkbounds(eos, rho, T, Y_e);
  if( error != 0 ) {
    char message[256];
    switch(error) {
      case 101:
        sprintf(message, "Input Y_e (%.15e) is too large.", Y_e);
        break;
      case 102:
        sprintf(message, "Input Y_e (%.15e) is too small.", Y_e);
        break;
      case 103:
        sprintf(message, "Input temperature (%.15e) is too large.", T);
        break;
      case 104:
        sprintf(message, "Input temperature (%.15e) is too small.", T);
        break;
      case 105:
        sprintf(message, "Input rho (%.15e) is too large.", rho);
        break;
      case 106:
        sprintf(message, "Input rho (%.15e) is too small.", rho);
        break;
    }
    sprintf(report->message,
            "In %s call, problem with checkbounds: %s\n",
            __func__, message);
    return error;
  }

  // Get interpolation spots
  int idx[8];
  double delx,dely,delz;
  const double lr = log(rho);
  const double lt = log(T);
  NRPyEOS_get_interp_spots(eos, lr, lt, Y_e, &delx, &dely, &delz, idx);

  for(int i=0; i<n; i++) {
    // Now perform the interpolations
    int key = tablevars_keys[i];
    double tablevar_out;
    NRPyEOS_linterp_one(eos, idx, delx, dely, delz, &tablevar_out, key);

    // We have the result, but we must convert appropriately.
    // The only edge cases are P and eps, for which we obtain
    // log(P) and log(eps+eps0). We must check for them here
    if( key == NRPyEOS_press_key ) {
      tablevar_out = exp(tablevar_out);
    }
    else if( key == NRPyEOS_eps_key ) {
      tablevar_out = exp(tablevar_out) - eos->energy_shift;
    }

    // Then update tablevars
    tablevars[i] = tablevar_out;
  }
  return 0;
#endif
}
