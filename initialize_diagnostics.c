#include "con2prim_header.h"

/* This function fills the struct con2prim_diagnostics with data.
   For more information on the arguments, see the definition of the
   struct in new_header.h. */
void initialize_diagnostics(con2prim_diagnostics *restrict diagnostic) {
  diagnostic->failures=0;
  diagnostic->font_fixes=0;
  diagnostic->vel_limited_ptcount=0;
  diagnostic->atm_resets=0;
  diagnostic->rho_star_fix_applied=0;
  diagnostic->pointcount=0;
  diagnostic->failures_inhoriz=0;
  diagnostic->pointcount_inhoriz=0;
  diagnostic->backup[0]=0;
  diagnostic->backup[1]=0;
  diagnostic->backup[2]=0;
  diagnostic->nan_found=0;
  diagnostic->c2p_fail_flag=0;
  diagnostic->error_int_numer=0;
  diagnostic->error_int_denom=0;
} 
