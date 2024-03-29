#include "ghl.h"

/*
 * Function     : ghl_initialize_conservatives()
 * Description  : Initialize the ghl_conservative_quantities struct from user input
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_initialize_conservatives
*/

void ghl_initialize_conservatives(
      const double rho,
      const double tau,
      const double S_x,
      const double S_y,
      const double S_z,
      const double entropy,
      const double Y_e,
      ghl_conservative_quantities *restrict cons) {

  cons->rho = rho;
  cons->tau = tau;
  cons->SD[0] = S_x;
  cons->SD[1] = S_y;
  cons->SD[2] = S_z;
  cons->entropy = entropy;
  cons->Y_e = Y_e;
}
