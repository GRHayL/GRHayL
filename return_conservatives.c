#include "con2prim_header.h"

/* This function fills the struct conservative_quantities with data. 
   For more information on the arguments, see the definition of the
   struct in new_header.h. */
void return_conservatives(
             const GRMHD_parameters *restrict params,
             const eos_parameters *restrict eos,
             const conservative_quantities *restrict cons,
             double *restrict rho, double *restrict tau,
             double *restrict S_x, double *restrict S_y, double *restrict S_z,
             double *restrict Y_e, double *restrict entropy) {

  *rho = cons->rho;
  *S_x = cons->S_x;
  *S_y = cons->S_y;
  *S_z = cons->S_z;
  *tau = cons->tau;
  if( eos->type == 1) {
    *Y_e = cons->Y_e;
  }
  if( params->evolve_entropy ) {
    *entropy = cons->entropy;
  }
}
