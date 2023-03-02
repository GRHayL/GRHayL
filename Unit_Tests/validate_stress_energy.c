#include "unit_tests.h"

/*
  In the function validate_stress_energy, the components of Tmunu are
  validated individually. In the case of failure, the code errors out
  and lists all the components that failed to be within the error
  bars of the perturbed version of the code.
*/

void validate_stress_energy(
                     const stress_energy *restrict Tmunu_trusted,
                     const stress_energy *restrict Tmunu,
                     const stress_energy *restrict Tmunu_pert) {


  char fail_msg[100] = "Test has failed!\n The stress-energy variable(s) which failed are ";
  int test_fail = 0;
  if( validate(Tmunu_trusted->Ttt, Tmunu->Ttt, Tmunu_pert->Ttt) && fabs(Tmunu->Ttt) > 1.0e-12 ) {
    printf("Ttt trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Ttt, Tmunu->Ttt, Tmunu_pert->Ttt);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Ttt, Tmunu->Ttt), relative_error(Tmunu_trusted->Ttt, Tmunu_pert->Ttt));
    sprintf(fail_msg, "%.90s Ttt", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->Ttx, Tmunu->Ttx, Tmunu_pert->Ttx) && fabs(Tmunu->Ttx) > 1.0e-12 ) {
    printf("Ttx trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Ttx, Tmunu->Ttx, Tmunu_pert->Ttx);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Ttx, Tmunu->Ttx), relative_error(Tmunu_trusted->Ttx, Tmunu_pert->Ttx));
    sprintf(fail_msg, "%.90s Ttx", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->Tty, Tmunu->Tty, Tmunu_pert->Tty) && fabs(Tmunu->Tty) > 1.0e-12 ) {
    printf("Tty trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Tty, Tmunu->Tty, Tmunu_pert->Tty);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Tty, Tmunu->Tty), relative_error(Tmunu_trusted->Tty, Tmunu_pert->Tty));
    sprintf(fail_msg, "%.90s Tty", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->Ttz, Tmunu->Ttz, Tmunu_pert->Ttz) && fabs(Tmunu->Ttz) > 1.0e-12 ) {
    printf("Ttz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Ttz, Tmunu->Ttz, Tmunu_pert->Ttz);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Ttz, Tmunu->Ttz), relative_error(Tmunu_trusted->Ttz, Tmunu_pert->Ttz));
    sprintf(fail_msg, "%.90s Ttz", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->Txx, Tmunu->Txx, Tmunu_pert->Txx) && fabs(Tmunu->Txx) > 1.0e-12 ) {
    printf("Txx trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Txx, Tmunu->Txx, Tmunu_pert->Txx);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Txx, Tmunu->Txx), relative_error(Tmunu_trusted->Txx, Tmunu_pert->Txx));
    sprintf(fail_msg, "%.90s Txx", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->Txy, Tmunu->Txy, Tmunu_pert->Txy) && fabs(Tmunu->Txy) > 1.0e-12 ) {
    printf("Txy trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Txy, Tmunu->Txy, Tmunu_pert->Txy);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Txy, Tmunu->Txy), relative_error(Tmunu_trusted->Txy, Tmunu_pert->Txy));
    sprintf(fail_msg, "%.90s Txy", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->Txz, Tmunu->Txz, Tmunu_pert->Txz) && fabs(Tmunu->Txz) > 1.0e-12 ) {
    printf("Txz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Txz, Tmunu->Txz, Tmunu_pert->Txz);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Txz, Tmunu->Txz), relative_error(Tmunu_trusted->Txz, Tmunu_pert->Txz));
    sprintf(fail_msg, "%.90s Txz", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->Tyy, Tmunu->Tyy, Tmunu_pert->Tyy) && fabs(Tmunu->Tyy) > 1.0e-12 ) {
    printf("Tyy trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Tyy, Tmunu->Tyy, Tmunu_pert->Tyy);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Tyy, Tmunu->Tyy), relative_error(Tmunu_trusted->Tyy, Tmunu_pert->Tyy));
    sprintf(fail_msg, "%.90s Tyy", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->Tyz, Tmunu->Tyz, Tmunu_pert->Tyz) && fabs(Tmunu->Tyz) > 1.0e-12 ) {
    printf("Tyz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Tyz, Tmunu->Tyz, Tmunu_pert->Tyz);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Tyz, Tmunu->Tyz), relative_error(Tmunu_trusted->Tyz, Tmunu_pert->Tyz));
    sprintf(fail_msg, "%.90s Tyz", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->Tzz, Tmunu->Tzz, Tmunu_pert->Tzz) && fabs(Tmunu->Tzz) > 1.0e-12 ) {
    printf("Tzz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->Tzz, Tmunu->Tzz, Tmunu_pert->Tzz);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->Tzz, Tmunu->Tzz), relative_error(Tmunu_trusted->Tzz, Tmunu_pert->Tzz));
    sprintf(fail_msg, "%.90s Tzz", fail_msg);
    test_fail = 1;
  }

  if(test_fail) {
    grhayl_error("%.100s\n", fail_msg);
  }
  return;
}
