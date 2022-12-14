#include "../EOS_tabulated.h"
#include "NRPyEOS_tabulated_helpers.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(const eos_parameters *restrict eos_params,
                                                    const int n,
                                                    const double rho,
                                                    const double Y_e,
                                                    const double T,
                                                    const int *restrict tablevars_keys,
                                                    double *restrict tablevars,
                                                    NRPyEOS_error_report *restrict report) {


  // This function will interpolate n table quantities from
  // (rho,Ye,T). It replaces EOS_Omni calls with keytemp = 1
  if( n > NRPyEOS_ntablekeys ) {
    fprintf(stderr,"(GRHayL - EOS) from_rho_Ye_T_interpolate_n_quantities: number of quantities exceed maximum allowed: %d > %d. ABORTING.",
            n,NRPyEOS_ntablekeys);
  }

  // Start by assuming no errors
  report->error = false;

  // Check table bounds for input variables
  report->error_key = NRPyEOS_checkbounds(eos_params,rho,T,Y_e);
  if( report->error_key != 0 ) {
    // This should never happen, because we enforce
    // limits before calling this function
    sprintf(report->message,"from_rho_Ye_T_interpolate_n_quantities: problem with checkbounds");
    report->error = true;
    return;
  }

  // Get interpolation spots
  int idx[8];
  double delx,dely,delz;
  const double lr = log(rho);
  const double lt = log(T);
  NRPyEOS_get_interp_spots(eos_params,lr,lt,Y_e,&delx,&dely,&delz,idx);

  for(int i=0;i<n;i++) {
    // Now perform the interpolations
    int key = tablevars_keys[i];
    double tablevar_out;
    NRPyEOS_linterp_one(eos_params,idx,delx,dely,delz,&tablevar_out,key);

    // We have the result, but we must convert appropriately.
    // The only edge cases are P and eps, for which we obtain
    // log(P) and log(eps+eps0). We must check for them here
    if( key == NRPyEOS_press_key ) {
      tablevar_out = exp(tablevar_out);
    }
    else if( key == NRPyEOS_eps_key ) {
      tablevar_out = exp(tablevar_out) - eos_params->energy_shift;
    }

    // Then update tablevars
    tablevars[i] = tablevar_out;
  }
}
