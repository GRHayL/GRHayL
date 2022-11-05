#include "con2prim_header.h"

/* This function fills the struct conservative_quantities with data. 
   For more information on the arguments, see the definition of the
   struct in new_header.h. */
void initialize_conservatives(
             const GRMHD_parameters *restrict params,
             const eos_parameters *restrict eos,
             const double rho, const double tau,
             const double S_x, const double S_y, const double S_z,
             const double Y_e, const double entropy,
             conservative_quantities *restrict cons) {
  cons->rho = rho;
  cons->S_x = S_x;
  cons->S_y = S_y;
  cons->S_z = S_z;
  cons->tau = tau;
  if( eos->eos_type == 1) {
    cons->Y_e = Y_e;
  }
  if( params->evolve_entropy ) {
    cons->entropy = entropy;
  }
}
