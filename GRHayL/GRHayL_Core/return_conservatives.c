#include "GRHayL.h"

/* Function    : return_conservatives()
 * Description : unpacks conservative struct into variables
 *
 * Intputs      : cons          - conservative_quantities struct to be unpacked
 *
 * Outputs     : rho            - pointer to rho_star
 *             : tau            - pointer to \tilde{tau}
 *             : S_x            - pointer to \tilde{S}_x
 *             : S_y            - pointer to \tilde{S}_y
 *             : S_z            - pointer to \tilde{S}_z
 *             : entropy        - pointer to densitized entropy
 */

void return_conservatives(
             const conservative_quantities *restrict cons,
             double *restrict rho, double *restrict tau,
             double *restrict S_x, double *restrict S_y, double *restrict S_z,
             double *restrict Y_e, double *restrict entropy) {

  *rho = cons->rho;
  *S_x = cons->S_x;
  *S_y = cons->S_y;
  *S_z = cons->S_z;
  *tau = cons->tau;
  *Y_e = cons->Y_e;
  *entropy = cons->entropy;
}
