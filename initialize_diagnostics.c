#include "con2prim_header.h"

/* This function fills the struct con2prim_diagnostics with data.
   For more information on the arguments, see the definition of the
   struct in new_header.h. */
void initialize_diagnostics(con2prim_diagnostics *restrict diagnostics) {
  diagnostics->failures=0;
  diagnostics->font_fixes=0;
  diagnostics->vel_limited_ptcount=0;
  diagnostics->atm_resets=0;
  diagnostics->rho_star_fix_applied=0;
  diagnostics->pointcount=0;
  diagnostics->failures_inhoriz=0;
  diagnostics->pointcount_inhoriz=0;
  diagnostics->backup[0]=0;
  diagnostics->backup[1]=0;
  diagnostics->backup[2]=0;
  diagnostics->nan_found=0;
  diagnostics->c2p_fail_flag=0;
  diagnostics->error_int_numer=0;
  diagnostics->error_int_denom=0;
  diagnostics->which_routine = None;
} 
