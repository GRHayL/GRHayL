#include "nrpyeos_tabulated.h"
#include "NRPyEOS_tabulated_helpers.h"
/*
 * (c) 2022 Leo Werneck
 */
int NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(
      const eos_parameters *restrict eos,
      const int n,
      const double prec,
      const double rho,
      const double Y_e,
      const double tablevar_in,
      const int tablevar_in_key,
      const int *restrict tablevars_keys,
      double *restrict tablevars,
      double *restrict T,
      NRPyEOS_error_report *restrict report) {
#ifndef GRHAYL_USE_HDF5
  HDF5_ERROR_IF_USED;
#else

  if( !n ) return;

  // Start by assuming no errors
  int error = 0;

  // This function will interpolate n table quantities from
  // (rho,Ye,aux). It replaces EOS_Omni calls with keytemp != 1
  if( n > NRPyEOS_ntablekeys ) {
    sprintf(report->message, "In %s call, number of quantities exceed maximum allowed: %d > %d.\n",
            __func__, n, NRPyEOS_ntablekeys);
    return 1;
  }

  // Check table bounds for input variables
  error = NRPyEOS_checkbounds_kt0_noTcheck(eos,rho,Y_e);
  if( error != 0 ) {
    char message[256];
    switch(error) {
      case 101:
        sprintf(message, "Input Y_e (%.15e) is too large.", Y_e);
        break;
      case 102:
        sprintf(message, "Input Y_e (%.15e) is too small.", Y_e);
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

  // First step is to recover the temperature. The variable
  // tablevar_in is the one used in the temperature recovery.
  // For example, if tablevar_in = eps, then we recover T
  // using (rho,Ye,eps).
  double aux = tablevar_in;

  if( tablevar_in_key == NRPyEOS_press_key ) {
    // If aux = P, then we need log(P).
    aux = log(aux);
  }
  else if( tablevar_in_key == NRPyEOS_eps_key ) {
    // If aux = eps, then we need log(eps+eps0).
    // Compute eps+eps0
    aux += eos->energy_shift;
    // At this point, aux *must* be positive. If not, error out.
    if( aux < 0.0 ) {
      sprintf(report->message, "In %s call, found eps+energy_shift < 0.0 (%e).\n",
              __func__, aux);
      return 2;
    }

    // Compute log(eps+eps0)
    aux = log(aux);
  }

  // Now compute the temperature
  const double lr  = log(rho);
  const double lt0 = log(*T);
  double lt        = 0.0;
  NRPyEOS_findtemp_from_any(eos, tablevar_in_key, lr, lt0, Y_e, aux, prec, &lt, &error);

  // Now set the temperature
  *T = exp(lt);

  // Then interpolate the quantities we want from (rho,Ye,T)
  int anyerr=0;
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos, n, rho, Y_e, *T, tablevars_keys, tablevars, report);
  return error;
#endif
}
