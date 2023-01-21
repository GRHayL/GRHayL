#include "unit_tests.h"

/*
  In the function validate_primitives, the primitives are validated
  individually, and the function returns upon finding the first failure.
  The failure codes are

    -1: original value differs by more than tolerance. Since this is the
        value from before the tested function, the source of the difference
        occurs due to an unknown problem before the tested function.
     1: vx differs by more than tolerance
     2: vy differs by more than tolerance
     3: vz differs by more than tolerance
     4: rho_b differs by more than tolerance
     5: pressure differs by more than tolerance
     6: entropy differs by more than tolerance
     7: Y_e differs by more than tolerance
     8: temperature differs by more than tolerance
*/

static int rel_tol(const double tolerance, const double x1, const double x2) {
  double rel_diff;
  if(x1!=0) {
    rel_diff = abs((x1-x2)/x1);
  } else if(x1!=0) {
    rel_diff = abs((x1-x2)/x2);
  } else {
    rel_diff = 0.0;
  }
  if(rel_diff > tolerance) return 1;
  return 0;
}

int validate_primitives(
      const double tolerance,
      const int eos_type,
      const bool velocity_only,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims_orig,
      const primitive_quantities *restrict prims,
      FILE *restrict infile) {

  primitive_quantities prims_before, prims_after;
  read_primitive_binary(eos_type, velocity_only, evolve_entropy, &prims_before, &prims_after, infile);

  if(rel_tol(tolerance, prims_before.vx, prims_orig->vx))
    return -1;

  if(rel_tol(tolerance, prims_before.vy, prims_orig->vy))
    return -1;

  if(rel_tol(tolerance, prims_before.vz, prims_orig->vz))
    return -1;

  if(!velocity_only) {
    if(rel_tol(tolerance, prims_before.rho, prims_orig->rho))
      return -1;

    if(rel_tol(tolerance, prims_before.press, prims_orig->press))
      return -1;

    if(evolve_entropy)
    if(rel_tol(tolerance, prims_before.entropy, prims_orig->entropy))
      return -1;

    if(eos_type == 2) { //Tabulated
      if(rel_tol(tolerance, prims_before.Y_e, prims_orig->Y_e))
        return -1;
  
      if(rel_tol(tolerance, prims_before.temperature, prims_orig->temperature))
        return -1;
    }
  }

  if(rel_tol(tolerance, prims_after.vx, prims->vx))
    return 1;

  if(rel_tol(tolerance, prims_after.vy, prims->vy))
    return 2;

  if(rel_tol(tolerance, prims_after.vz, prims->vz))
    return 3;

  if(!velocity_only) {
    if(rel_tol(tolerance, prims_after.rho, prims->rho))
      return 4;

    if(rel_tol(tolerance, prims_after.press, prims->press))
      return 5;

    if(evolve_entropy)
      if(rel_tol(tolerance, prims_after.entropy, prims->entropy))
        return 6;

    if(eos_type == 2) { //Tabulated
      if(rel_tol(tolerance, prims_after.Y_e, prims->Y_e))
        return 7;

      if(rel_tol(tolerance, prims_after.temperature, prims->temperature))
        return 8;
    }
  }
  return 0;
}
