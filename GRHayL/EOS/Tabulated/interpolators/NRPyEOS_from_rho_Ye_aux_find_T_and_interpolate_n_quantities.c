#include "ghl_nrpyeos_tabulated.h"
#include "NRPyEOS_tabulated_helpers.h"
/*
 * (c) 2022 Leo Werneck
 */
ghl_error_codes_t NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(
      const ghl_eos_parameters *restrict eos,
      const int n,
      const double prec,
      const double rho,
      const double Y_e,
      const double tablevar_in,
      const int tablevar_in_key,
      const int *restrict tablevars_keys,
      double *restrict tablevars,
      double *restrict T) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
  return ghl_error_hdf5_is_disabled;
#else

  // This function will interpolate n table quantities from
  // (rho,Ye,aux). It replaces EOS_Omni calls with keytemp != 1
  if(n > NRPyEOS_ntablekeys) 
    return ghl_error_exceed_table_vars;

  // Check table bounds for input variables
  ghl_error_codes_t error = NRPyEOS_checkbounds_kt0_noTcheck(eos, rho, Y_e);
  if(error)
    return error;

  // First step is to recover the temperature. The variable
  // tablevar_in is the one used in the temperature recovery.
  // For example, if tablevar_in = eps, then we recover T
  // using (rho,Ye,eps).
  double aux = tablevar_in;

  if(tablevar_in_key == NRPyEOS_press_key) {
    // If aux = P, then we need log(P).
    aux = log(aux);
  }
  else if(tablevar_in_key == NRPyEOS_eps_key) {
    // If aux = eps, then we need log(eps+eps0).
    // Compute eps+eps0
    aux += eos->energy_shift;
    // At this point, aux *must* be positive. If not, error out.
    if(aux < 0.0) 
      return ghl_error_table_neg_energy;

    // Compute log(eps+eps0)
    aux = log(aux);
  }

  // Now compute the temperature
  const double lr  = log(rho);
  const double lt0 = log(*T);
  double lt        = 0.0;
  error = NRPyEOS_findtemp_from_any(eos, tablevar_in_key, lr, lt0, Y_e, aux, prec, &lt);
  if(error)
    return error;

  // Now set the temperature
  *T = exp(lt);

  // Then interpolate the quantities we want from (rho,Ye,T)
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos, n, rho, Y_e, *T, tablevars_keys, tablevars);
  return ghl_success;
#endif
}
