#include "ghl_con2prim.h"

/* Function    : ghl_initialize_diagnostics()
 * Description : Initialize the diagnostics struct for Con2Prim
 *
 * Outputs     : diagnostics   - new ghl_con2prim_diagnostics struct is
 *                               initialized with zeros in preparation
 *                               for tracking Con2Prim diagnostics
 */

void ghl_initialize_diagnostics(ghl_con2prim_diagnostics *restrict diagnostics) {
  diagnostics->tau_fix = false;
  diagnostics->Stilde_fix = false;
  diagnostics->speed_limited = 0;
  diagnostics->backup[0] = false;
  diagnostics->backup[1] = false;
  diagnostics->backup[2] = false;
  diagnostics->which_routine = None;
}
