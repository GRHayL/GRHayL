#include "con2prim.h"

/* Function    : initialize_diagnostics()
 * Description : Initialize the diagnostics struct for Con2Prim
 *
 * Outputs     : diagnostics   - new con2prim_diagnostics struct is
 *                               initialized with zeros in preparation
 *                               for tracking Con2Prim diagnostics
 */

void initialize_diagnostics(con2prim_diagnostics *restrict diagnostics) {
  diagnostics->failures=0;
  diagnostics->font_fixes=0;
  diagnostics->vel_limited_ptcount=0;
  diagnostics->backup[0]=0;
  diagnostics->backup[1]=0;
  diagnostics->backup[2]=0;
  diagnostics->nan_found=0;
  diagnostics->c2p_fail_flag=0;
  diagnostics->error_int_numer=0;
  diagnostics->error_int_denom=0;
  diagnostics->which_routine = None;
}