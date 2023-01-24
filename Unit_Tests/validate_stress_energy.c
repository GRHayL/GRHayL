#include "unit_tests.h"

/*
  In the function validate_stress_energy, the components of Tmunu are
  validated individually, and the function returns upon finding the
  first failure. The failure codes are

    -1: original value differs by more than tolerance. Since this is the
        value from before the tested function, the source of the difference
        occurs due to an unknown problem before the tested function.
     1: Ttt differs by more than tolerance
     2: Ttx differs by more than tolerance
     3: Tty differs by more than tolerance
     4: Ttz differs by more than tolerance
     5: Txx differs by more than tolerance
     6: Txy differs by more than tolerance
     7: Txz differs by more than tolerance
     8: Tyy differs by more than tolerance
     9: Tyz differs by more than tolerance
    10: Tzz differs by more than tolerance
*/

int validate_stress_energy(
                     const double tolerance,
                     const stress_energy *restrict Tmunu_orig,
                     const stress_energy *restrict Tmunu,
                     FILE *restrict infile) {

  stress_energy Tmunu_before, Tmunu_after;
  read_stress_energy_binary(&Tmunu_before, infile);
  read_stress_energy_binary(&Tmunu_after, infile);

  if(rel_tol(tolerance, Tmunu_before.Ttt, Tmunu_orig->Ttt))
    return -1;

  if(rel_tol(tolerance, Tmunu_before.Ttx, Tmunu_orig->Ttx))
    return -1;

  if(rel_tol(tolerance, Tmunu_before.Tty, Tmunu_orig->Tty))
    return -1;

  if(rel_tol(tolerance, Tmunu_before.Ttz, Tmunu_orig->Ttz))
    return -1;

  if(rel_tol(tolerance, Tmunu_before.Txx, Tmunu_orig->Txx))
    return -1;

  if(rel_tol(tolerance, Tmunu_before.Txy, Tmunu_orig->Txy))
    return -1;

  if(rel_tol(tolerance, Tmunu_before.Txz, Tmunu_orig->Txz))
    return -1;

  if(rel_tol(tolerance, Tmunu_before.Tyy, Tmunu_orig->Tyy))
    return -1;

  if(rel_tol(tolerance, Tmunu_before.Tyz, Tmunu_orig->Tyz))
    return -1;

  if(rel_tol(tolerance, Tmunu_before.Tzz, Tmunu_orig->Tzz))
    return -1;

  if(rel_tol(tolerance, Tmunu_after.Ttt, Tmunu->Ttt))
    return 1;

  if(rel_tol(tolerance, Tmunu_after.Ttx, Tmunu->Ttx))
    return 2;

  if(rel_tol(tolerance, Tmunu_after.Tty, Tmunu->Tty))
    return 3;

  if(rel_tol(tolerance, Tmunu_after.Ttz, Tmunu->Ttz))
    return 4;

  if(rel_tol(tolerance, Tmunu_after.Txx, Tmunu->Txx))
    return 5;

  if(rel_tol(tolerance, Tmunu_after.Txy, Tmunu->Txy))
    return 6;

  if(rel_tol(tolerance, Tmunu_after.Txz, Tmunu->Txz))
    return 7;

  if(rel_tol(tolerance, Tmunu_after.Tyy, Tmunu->Tyy))
    return 8;

  if(rel_tol(tolerance, Tmunu_after.Tyz, Tmunu->Tyz))
    return 9;

  if(rel_tol(tolerance, Tmunu_after.Tzz, Tmunu->Tzz))
    return 10;

  return 0;
}
