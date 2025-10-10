#include "ghl_nrpyeos_tabulated.h"
#include "NRPyEOS_tabulated_helpers.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(
      const ghl_eos_parameters *restrict eos,
      const int n,
      const double rho,
      const double Y_e,
      const double T,
      const int *restrict tablevars_keys,
      double *restrict tablevars) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
  return ghl_error_hdf5_is_disabled;
#else

  if(!n) return ghl_success;

  // This function will interpolate n table quantities from
  // (rho,Ye,T). It replaces EOS_Omni calls with keytemp = 1
  if(n > NRPyEOS_ntablekeys) 
    return ghl_error_exceed_table_vars;

  // Check table bounds for input variables
  const ghl_error_codes_t error = NRPyEOS_checkbounds(eos, rho, T, Y_e);
  if(error)
    return error;

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
    if(key == NRPyEOS_press_key) {
      tablevar_out = exp(tablevar_out);
    } else if(key == NRPyEOS_eps_key) {
      tablevar_out = exp(tablevar_out) - eos->energy_shift;
    }

    // Then update tablevars
    tablevars[i] = tablevar_out;
  }
  return ghl_success;
#endif
}
