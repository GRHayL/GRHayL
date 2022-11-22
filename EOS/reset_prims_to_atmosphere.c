// These were taken from Leo Werneck's IllinoisGRMHD code for con2prim to work,
// but they should actually be provided by the EOS code.

#include "con2prim_header.h"
#include "EOS_hybrid_header.h"

void reset_prims_to_atmosphere( const eos_parameters *restrict eos,
                                primitive_quantities *restrict prims,
                                con2prim_diagnostics *restrict diagnostics ) {

  // Just a simple reset to atmospheric values.
  // Velocities are set to zero. Keeping it
  // inside a single function ensures that
  // resets are consistent throughout the code.
  prims->rho = eos->rho_atm;
  prims->press = eos->press_atm;
  prims->eps = eos->eps_atm;
  prims->entropy = eos->entropy_atm;
  if( eos->eos_type == 1 ) {
    prims->Y_e = eos->Ye_atm;
    prims->temp = eos->temp_atm;
  }
  prims->vx = 0.0;
  prims->vy = 0.0;
  prims->vz = 0.0;  
}
