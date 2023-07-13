#include "con2prim.h"

/* Function    : ghl_initialize_diagnostics()
 * Description : Initialize the diagnostics struct for Con2Prim
 *
 * Outputs     : diagnostics   - new ghl_con2prim_diagnostics struct is
 *                               initialized with zeros in preparation
 *                               for tracking Con2Prim diagnostics
 */

void ghl_initialize_diagnostics(ghl_con2prim_diagnostics *restrict diagnostics) {
  diagnostics->failures=0;
  diagnostics->speed_limited=0;
  diagnostics->backup[0]=0;
  diagnostics->backup[1]=0;
  diagnostics->backup[2]=0;
  diagnostics->nan_found=0;
  diagnostics->error_int_numer=0;
  diagnostics->error_int_denom=0;
  diagnostics->which_routine = None;
}
