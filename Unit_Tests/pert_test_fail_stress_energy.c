#include "ghl_unit_tests.h"

/*
  In the function ghl_pert_test_fail_stress_energy, the components of Tmunu are
  ghl_pert_test_faild individually. In the case of failure, the code errors out
  and lists all the components that failed to be within the error
  bars of the perturbed version of the code.
*/

void ghl_pert_test_fail_stress_energy(
      const ghl_stress_energy *restrict Tmunu_trusted,
      const ghl_stress_energy *restrict Tmunu,
      const ghl_stress_energy *restrict Tmunu_pert) {

  bool test_fail = false;
  if( ghl_pert_test_fail(Tmunu_trusted->T4[0][0], Tmunu->T4[0][0], Tmunu_pert->T4[0][0]) && fabs(Tmunu->T4[0][0]) > 1.0e-12 ) {
    printf("Ttt trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[0][0], Tmunu->T4[0][0], Tmunu_pert->T4[0][0]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[0][0], Tmunu->T4[0][0]), relative_error(Tmunu_trusted->T4[0][0], Tmunu_pert->T4[0][0]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(Tmunu_trusted->T4[0][1], Tmunu->T4[0][1], Tmunu_pert->T4[0][1]) && fabs(Tmunu->T4[0][1]) > 1.0e-12 ) {
    printf("Ttx trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[0][1], Tmunu->T4[0][1], Tmunu_pert->T4[0][1]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[0][1], Tmunu->T4[0][1]), relative_error(Tmunu_trusted->T4[0][1], Tmunu_pert->T4[0][1]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(Tmunu_trusted->T4[0][2], Tmunu->T4[0][2], Tmunu_pert->T4[0][2]) && fabs(Tmunu->T4[0][2]) > 1.0e-12 ) {
    printf("Tty trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[0][2], Tmunu->T4[0][2], Tmunu_pert->T4[0][2]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[0][2], Tmunu->T4[0][2]), relative_error(Tmunu_trusted->T4[0][2], Tmunu_pert->T4[0][2]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(Tmunu_trusted->T4[0][3], Tmunu->T4[0][3], Tmunu_pert->T4[0][3]) && fabs(Tmunu->T4[0][3]) > 1.0e-12 ) {
    printf("Ttz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[0][3], Tmunu->T4[0][3], Tmunu_pert->T4[0][3]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[0][3], Tmunu->T4[0][3]), relative_error(Tmunu_trusted->T4[0][3], Tmunu_pert->T4[0][3]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(Tmunu_trusted->T4[1][1], Tmunu->T4[1][1], Tmunu_pert->T4[1][1]) && fabs(Tmunu->T4[1][1]) > 1.0e-12 ) {
    printf("Txx trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[1][1], Tmunu->T4[1][1], Tmunu_pert->T4[1][1]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[1][1], Tmunu->T4[1][1]), relative_error(Tmunu_trusted->T4[1][1], Tmunu_pert->T4[1][1]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(Tmunu_trusted->T4[1][2], Tmunu->T4[1][2], Tmunu_pert->T4[1][2]) && fabs(Tmunu->T4[1][2]) > 1.0e-12 ) {
    printf("Txy trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[1][2], Tmunu->T4[1][2], Tmunu_pert->T4[1][2]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[1][2], Tmunu->T4[1][2]), relative_error(Tmunu_trusted->T4[1][2], Tmunu_pert->T4[1][2]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(Tmunu_trusted->T4[1][3], Tmunu->T4[1][3], Tmunu_pert->T4[1][3]) && fabs(Tmunu->T4[1][3]) > 1.0e-12 ) {
    printf("Txz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[1][3], Tmunu->T4[1][3], Tmunu_pert->T4[1][3]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[1][3], Tmunu->T4[1][3]), relative_error(Tmunu_trusted->T4[1][3], Tmunu_pert->T4[1][3]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(Tmunu_trusted->T4[2][2], Tmunu->T4[2][2], Tmunu_pert->T4[2][2]) && fabs(Tmunu->T4[2][2]) > 1.0e-12 ) {
    printf("Tyy trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[2][2], Tmunu->T4[2][2], Tmunu_pert->T4[2][2]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[2][2], Tmunu->T4[2][2]), relative_error(Tmunu_trusted->T4[2][2], Tmunu_pert->T4[2][2]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(Tmunu_trusted->T4[2][3], Tmunu->T4[2][3], Tmunu_pert->T4[2][3]) && fabs(Tmunu->T4[2][3]) > 1.0e-12 ) {
    printf("Tyz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[2][3], Tmunu->T4[2][3], Tmunu_pert->T4[2][3]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[2][3], Tmunu->T4[2][3]), relative_error(Tmunu_trusted->T4[2][3], Tmunu_pert->T4[2][3]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(Tmunu_trusted->T4[3][3], Tmunu->T4[3][3], Tmunu_pert->T4[3][3]) && fabs(Tmunu->T4[3][3]) > 1.0e-12 ) {
    printf("Tzz trusted %.14e computed %.14e perturbed %.14e\n", Tmunu_trusted->T4[3][3], Tmunu->T4[3][3], Tmunu_pert->T4[3][3]);
    printf("rel.err. %.14e %.14e\n", relative_error(Tmunu_trusted->T4[3][3], Tmunu->T4[3][3]), relative_error(Tmunu_trusted->T4[3][3], Tmunu_pert->T4[3][3]));
    test_fail = true;
  }

  if(test_fail) {
    ghl_error("Test has failed!");
  }
  return;
}
