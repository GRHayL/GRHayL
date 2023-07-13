#include "unit_tests.h"

/*
  In the function validate_ghl_stress_energy, the components of Tmunu are
  validated individually. In the case of failure, the code errors out
  and lists all the components that failed to be within the error
  bars of the perturbed version of the code.
*/

void ghl_validate_ghl_stress_energy(
      const ghl_stress_energy *restrict Tmunu_trusted,
      const ghl_stress_energy *restrict Tmunu,
      const ghl_stress_energy *restrict Tmunu_pert) {

  char fail_msg[100] = "Test has failed!\n The stress-energy variable(s) which failed are ";
  int test_fail = 0;
  if( validate(Tmunu_trusted->T4[0][0], Tmunu->T4[0][0], Tmunu_pert->T4[0][0]) && fabs(Tmunu->T4[0][0]) > 1.0e-12 ) {
    printf("Ttt trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[0][0], Tmunu->T4[0][0], Tmunu_pert->T4[0][0]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[0][0], Tmunu->T4[0][0]), relative_error(Tmunu_trusted->T4[0][0], Tmunu_pert->T4[0][0]));
    sprintf(fail_msg, "%.90s Ttt", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->T4[0][1], Tmunu->T4[0][1], Tmunu_pert->T4[0][1]) && fabs(Tmunu->T4[0][1]) > 1.0e-12 ) {
    printf("Ttx trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[0][1], Tmunu->T4[0][1], Tmunu_pert->T4[0][1]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[0][1], Tmunu->T4[0][1]), relative_error(Tmunu_trusted->T4[0][1], Tmunu_pert->T4[0][1]));
    sprintf(fail_msg, "%.90s Ttx", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->T4[0][2], Tmunu->T4[0][2], Tmunu_pert->T4[0][2]) && fabs(Tmunu->T4[0][2]) > 1.0e-12 ) {
    printf("Tty trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[0][2], Tmunu->T4[0][2], Tmunu_pert->T4[0][2]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[0][2], Tmunu->T4[0][2]), relative_error(Tmunu_trusted->T4[0][2], Tmunu_pert->T4[0][2]));
    sprintf(fail_msg, "%.90s Tty", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->T4[0][3], Tmunu->T4[0][3], Tmunu_pert->T4[0][3]) && fabs(Tmunu->T4[0][3]) > 1.0e-12 ) {
    printf("Ttz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[0][3], Tmunu->T4[0][3], Tmunu_pert->T4[0][3]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[0][3], Tmunu->T4[0][3]), relative_error(Tmunu_trusted->T4[0][3], Tmunu_pert->T4[0][3]));
    sprintf(fail_msg, "%.90s Ttz", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->T4[1][1], Tmunu->T4[1][1], Tmunu_pert->T4[1][1]) && fabs(Tmunu->T4[1][1]) > 1.0e-12 ) {
    printf("Txx trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[1][1], Tmunu->T4[1][1], Tmunu_pert->T4[1][1]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[1][1], Tmunu->T4[1][1]), relative_error(Tmunu_trusted->T4[1][1], Tmunu_pert->T4[1][1]));
    sprintf(fail_msg, "%.90s Txx", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->T4[1][2], Tmunu->T4[1][2], Tmunu_pert->T4[1][2]) && fabs(Tmunu->T4[1][2]) > 1.0e-12 ) {
    printf("Txy trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[1][2], Tmunu->T4[1][2], Tmunu_pert->T4[1][2]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[1][2], Tmunu->T4[1][2]), relative_error(Tmunu_trusted->T4[1][2], Tmunu_pert->T4[1][2]));
    sprintf(fail_msg, "%.90s Txy", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->T4[1][3], Tmunu->T4[1][3], Tmunu_pert->T4[1][3]) && fabs(Tmunu->T4[1][3]) > 1.0e-12 ) {
    printf("Txz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[1][3], Tmunu->T4[1][3], Tmunu_pert->T4[1][3]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[1][3], Tmunu->T4[1][3]), relative_error(Tmunu_trusted->T4[1][3], Tmunu_pert->T4[1][3]));
    sprintf(fail_msg, "%.90s Txz", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->T4[2][2], Tmunu->T4[2][2], Tmunu_pert->T4[2][2]) && fabs(Tmunu->T4[2][2]) > 1.0e-12 ) {
    printf("Tyy trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[2][2], Tmunu->T4[2][2], Tmunu_pert->T4[2][2]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[2][2], Tmunu->T4[2][2]), relative_error(Tmunu_trusted->T4[2][2], Tmunu_pert->T4[2][2]));
    sprintf(fail_msg, "%.90s Tyy", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->T4[2][3], Tmunu->T4[2][3], Tmunu_pert->T4[2][3]) && fabs(Tmunu->T4[2][3]) > 1.0e-12 ) {
    printf("Tyz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[2][3], Tmunu->T4[2][3], Tmunu_pert->T4[2][3]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[2][3], Tmunu->T4[2][3]), relative_error(Tmunu_trusted->T4[2][3], Tmunu_pert->T4[2][3]));
    sprintf(fail_msg, "%.90s Tyz", fail_msg);
    test_fail = 1;
  }

  if( validate(Tmunu_trusted->T4[3][3], Tmunu->T4[3][3], Tmunu_pert->T4[3][3]) && fabs(Tmunu->T4[3][3]) > 1.0e-12 ) {
    printf("Tzz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[3][3], Tmunu->T4[3][3], Tmunu_pert->T4[3][3]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[3][3], Tmunu->T4[3][3]), relative_error(Tmunu_trusted->T4[3][3], Tmunu_pert->T4[3][3]));
    sprintf(fail_msg, "%.90s Tzz", fail_msg);
    test_fail = 1;
  }

  if(test_fail) {
    ghl_error("%.100s\n", fail_msg);
  }
  return;
}
