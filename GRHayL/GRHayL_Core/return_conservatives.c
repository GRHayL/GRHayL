#include "ghl.h"

/*
 * Function     : ghl_return_conservatives()
 * Description  : Unpacks ghl_conservative_quantities struct data into provided
 *                provided memory locations
 * Documentation: https://github.com/GRHayL/GRHayL/wiki/ghl_return_conservatives
*/

void ghl_return_conservatives(
      const ghl_conservative_quantities *restrict cons,
      double *restrict rho, double *restrict tau,
      double *restrict S_x, double *restrict S_y, double *restrict S_z,
      double *restrict Y_e, double *restrict entropy) {

  *rho = cons->rho;
  *S_x = cons->SD[0];
  *S_y = cons->SD[1];
  *S_z = cons->SD[2];
  *tau = cons->tau;
  *Y_e = cons->Y_e;
  *entropy = cons->entropy;
}
