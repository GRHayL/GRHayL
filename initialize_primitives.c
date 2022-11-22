#include "con2prim_header.h"

/* This function fills the struct primitive_quantities with data. Some codes,
   like IllinoisGRMHD, do not have initial guesses for most primitives. In this
   case, the parameter guess_prims should be set to 1. These values should all
   be initialized regardless, and set to zero if they are unavailable.
   For more information on the arguments, see the definition of the
   struct in new_header.h. */

void initialize_primitives(
             const eos_parameters *restrict eos,
             const metric_quantities *restrict metric,
             const double rho, const double press, const double epsilon,
             const double vx, const double vy, const double vz,
             const double Bx, const double By, const double Bz,
             const double entropy, const double Y_e, const double temp,
             primitive_quantities *restrict prims) {
  prims->rho = rho;
  prims->press = press;
  prims->vx = vx;
  prims->vy = vy;
  prims->vz = vz;
  prims->Bx = Bx;
  prims->By = By;
  prims->Bz = Bz;
  prims->eps = epsilon;
  prims->entropy = entropy/rho;
  // Tabulated EOS quantities
  if( eos->eos_type == 1 ) {
    prims->Y_e = Y_e/rho;
    prims->temp = temp/rho;
  }
  prims->print=false;
}
