#include "unit_tests.h"

/*
  In the function validate_stress_energy, the components of Tmunu are
  validated individually. In the case of failure, the code errors out
  and lists all the components that failed to be within the error
  bars of the perturbed version of the code.
*/

void validate_stress_energy(
                     const stress_energy *restrict Tmunu,
                     const stress_energy *restrict Tmunu_trusted,
                     const stress_energy *restrict Tmunu_pert) {


  char fail_msg[100] = "Test has failed!\n The stress-energy variable(s) which failed are ";
  int test_fail = 0;
  if( relative_error(Tmunu_trusted->Ttt, Tmunu->Ttt) > relative_error(Tmunu_trusted->Ttt, Tmunu_pert->Ttt) ) {
    sprintf(fail_msg, "%.90s Ttt", fail_msg);
    test_fail = 1;
  }

  if( relative_error(Tmunu_trusted->Ttx, Tmunu->Ttx) > relative_error(Tmunu_trusted->Ttx, Tmunu_pert->Ttx) ) {
    sprintf(fail_msg, "%.90s Ttx", fail_msg);
    test_fail = 1;
  }

  if( relative_error(Tmunu_trusted->Tty, Tmunu->Tty) > relative_error(Tmunu_trusted->Tty, Tmunu_pert->Tty) ) {
    sprintf(fail_msg, "%.90s Tty", fail_msg);
    test_fail = 1;
  }

  if( relative_error(Tmunu_trusted->Ttz, Tmunu->Ttz) > relative_error(Tmunu_trusted->Ttz, Tmunu_pert->Ttz) ) {
    sprintf(fail_msg, "%.90s Ttz", fail_msg);
    test_fail = 1;
  }

  if( relative_error(Tmunu_trusted->Txx, Tmunu->Txx) > relative_error(Tmunu_trusted->Txx, Tmunu_pert->Txx) ) {
    sprintf(fail_msg, "%.90s Txx", fail_msg);
    test_fail = 1;
  }

  if( relative_error(Tmunu_trusted->Txy, Tmunu->Txy) > relative_error(Tmunu_trusted->Txy, Tmunu_pert->Txy) ) {
    sprintf(fail_msg, "%.90s Txy", fail_msg);
    test_fail = 1;
  }

  if( relative_error(Tmunu_trusted->Txz, Tmunu->Txz) > relative_error(Tmunu_trusted->Txz, Tmunu_pert->Txz) ) {
    sprintf(fail_msg, "%.90s Txz", fail_msg);
    test_fail = 1;
  }

  if( relative_error(Tmunu_trusted->Tyy, Tmunu->Tyy) > relative_error(Tmunu_trusted->Tyy, Tmunu_pert->Tyy) ) {
    sprintf(fail_msg, "%.90s Tyy", fail_msg);
    test_fail = 1;
  }

  if( relative_error(Tmunu_trusted->Tyz, Tmunu->Tyz) > relative_error(Tmunu_trusted->Tyz, Tmunu_pert->Tyz) ) {
    sprintf(fail_msg, "%.90s Tyz", fail_msg);
    test_fail = 1;
  }

  if( relative_error(Tmunu_trusted->Tzz, Tmunu->Tzz) > relative_error(Tmunu_trusted->Tzz, Tmunu_pert->Tzz) ) {
    sprintf(fail_msg, "%.90s Tzz", fail_msg);
    test_fail = 1;
  }

  if(test_fail) {
    grhayl_error("%.100s\n", fail_msg);
  }
  return;
}
