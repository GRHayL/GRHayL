#include "con2prim_header.h"

/* Function    : return_conservatives()
 * Authors     : Samuel Cupp
 * Description : unpacks conservative struct into variables
 * Dependencies: None
 *
 * Intputs      : cons           - conservative_quantities struct with
 *                                data from the Con2Prim routine
 *
 * Outputs     : rho            - pointer to densitized density variable
 *             : tau            - pointer to densitized energy variable
 *             : S_x            - pointer to the x component of densitized
 *                                momentum variable
 *             : S_y            - pointer to the y component of densitized
                                  momentum variable
 *             : S_z            - pointer to the z component of densitized
 *                                densitized momentum variable
 *             : entropy        - pointer to densitized entropy variable
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
