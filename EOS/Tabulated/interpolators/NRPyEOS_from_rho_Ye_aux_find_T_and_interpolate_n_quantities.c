#include "../NRPy_basic_defines.h"
#include "../NRPy_function_prototypes.h"
#include "NRPyEOS_tabulated_helpers.h"
/*
 * (c) 2022 Leo Werneck
 */
void NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(const NRPyEOS_params *restrict eos_params,
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


  // This function will interpolate n table quantities from
  // (rho,Ye,aux). It replaces EOS_Omni calls with keytemp != 1
  if( n > NRPyEOS_ntablekeys ) {
    fprintf(stderr,"(NRPyEOS) NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities: number of quantities exceed maximum allowed: %d > %d. ABORTING.",
            n,NRPyEOS_ntablekeys);
  }

  // Check table bounds for input variables
  report->error_key = NRPyEOS_checkbounds_kt0_noTcheck(eos_params,rho,Y_e);
  if( report->error_key != 0 ) {
    // This should never happen, because we enforce
    // limits before calling this function
    sprintf(report->message,"NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities: problem with checkbounds_kt0_noTcheck");
    report->error = true;
    return;
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
    aux += eos_params->energy_shift;
    // At this point, aux *must* be positive. If not, error out.
    if( aux < 0.0 ) {
      fprintf(stderr,"(NRPyEOS) NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities: found eps+energy_shift < 0.0 (%e). ABORTING.",
              aux);
    }
    // Compute log(eps+eps0)
    aux = log(aux);
  }

  // Now compute the temperature
  const double lr  = log(rho);
  const double lt0 = log(*T);
  double lt        = 0.0;
  int keyerr=0;
  NRPyEOS_findtemp_from_any(eos_params,tablevar_in_key,lr,lt0,Y_e,aux,prec,&lt,&keyerr);

  // Now set the temperature
  *T = exp(lt);

  // Then interpolate the quantities we want from (rho,Ye,T)
  int anyerr=0;
  NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(eos_params,n,rho,Y_e,*T,tablevars_keys,tablevars,report);
  report->error_key = keyerr;
  report->error     = anyerr;
}
