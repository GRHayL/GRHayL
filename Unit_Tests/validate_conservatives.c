#include "unit_tests.h"

/*
  In the function validate_conservatives, the conservatives are validated
  individually, and the function returns upon finding the first failure.
  The failure codes are

    -1: original value differs by more than tolerance. Since this is the
        value from before the tested function, the source of the difference
        occurs due to an unknown problem before the tested function.
     1: rho_star differs by more than tolerance
     2: tau differs by more than tolerance
     3: S_x differs by more than tolerance
     4: S_y differs by more than tolerance
     5: S_z differs by more than tolerance
     6: entropy differs by more than tolerance
*/

int validate_conservatives(
                     const double tolerance,
                     const bool evolve_entropy,
                     const conservative_quantities *restrict cons_orig,
                     const conservative_quantities *restrict cons,
                     FILE *restrict infile) {

  conservative_quantities cons_before, cons_after;
  read_conservative_binary(evolve_entropy, &cons_before, &cons_after, infile);

  if(rel_tol(tolerance, cons_before.rho, cons_orig->rho))
    return -1;

  if(rel_tol(tolerance, cons_before.tau, cons_orig->tau))
    return -1;

  if(rel_tol(tolerance, cons_before.S_x, cons_orig->S_x))
    return -1;

  if(rel_tol(tolerance, cons_before.S_y, cons_orig->S_y))
    return -1;

  if(rel_tol(tolerance, cons_before.S_z, cons_orig->S_z))
    return -1;

  if(evolve_entropy && rel_tol(tolerance, cons_before.entropy, cons_orig->entropy))
    return -1;

  if(rel_tol(tolerance, cons_after.rho, cons->rho))
    return 1;

  if(rel_tol(tolerance, cons_after.tau, cons->tau))
    return 2;

  if(rel_tol(tolerance, cons_after.S_x, cons->S_x))
    return 3;

  if(rel_tol(tolerance, cons_after.S_y, cons->S_y))
    return 4;

  if(rel_tol(tolerance, cons_after.S_z, cons->S_z))
    return 5;

  if(evolve_entropy && rel_tol(tolerance, cons_after.entropy, cons->entropy))
    return 6;

  return 0;
}
