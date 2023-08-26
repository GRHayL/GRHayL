#include "unit_tests.h"

/*
  In the function ghl_pert_test_fail_conservatives, the conservatives are ghl_pert_test_faild
  individually. In the case of failure, the code errors out and lists all
  the components that failed to be within the error bars of the perturbed
  version of the code.
*/

void ghl_pert_test_fail_conservatives(
      const bool evolve_entropy,
      const ghl_conservative_quantities *restrict cons_trusted,
      const ghl_conservative_quantities *restrict cons,
      const ghl_conservative_quantities *restrict cons_pert) {

  char fail_msg[100] = "Test has failed!\n The conservative variable(s) which failed are ";
  int test_fail = 0;
  if( ghl_pert_test_fail(cons_trusted->rho, cons->rho, cons_pert->rho) ) {
    printf("rho trusted %.14e computed %.14e perturbed %.14e\n", cons_trusted->rho, cons->rho, cons_pert->rho);
    printf("rel.err. %.14e %.14e\n", relative_error(cons_trusted->rho, cons->rho), relative_error(cons_trusted->rho, cons_pert->rho));
    sprintf(fail_msg, "%.90s rho_star", fail_msg);
    test_fail = 1;
  }

  if( ghl_pert_test_fail(cons_trusted->tau, cons->tau, cons_pert->tau) ) {
    printf("tau trusted %.14e computed %.14e perturbed %.14e\n", cons_trusted->tau, cons->tau, cons_pert->tau);
    printf("rel.err. %.14e %.14e\n", relative_error(cons_trusted->tau, cons->tau), relative_error(cons_trusted->tau, cons_pert->tau));
    sprintf(fail_msg, "%.90s tau", fail_msg);
    test_fail = 1;
  }

  if( ghl_pert_test_fail(cons_trusted->SD[0], cons->SD[0], cons_pert->SD[0]) ) {
    printf("S[0] trusted %.14e computed %.14e perturbed %.14e\n", cons_trusted->SD[0], cons->SD[0], cons_pert->SD[0]);
    printf("rel.err. %.14e %.14e\n", relative_error(cons_trusted->SD[0], cons->SD[0]), relative_error(cons_trusted->SD[0], cons_pert->SD[0]));
    sprintf(fail_msg, "%.90s S[0]", fail_msg);
    test_fail = 1;
  }

  if( ghl_pert_test_fail(cons_trusted->SD[1], cons->SD[1], cons_pert->SD[1]) ) {
    printf("S[1] trusted %.14e computed %.14e perturbed %.14e\n", cons_trusted->SD[1], cons->SD[1], cons_pert->SD[1]);
    printf("rel.err. %.14e %.14e\n", relative_error(cons_trusted->SD[1], cons->SD[1]), relative_error(cons_trusted->SD[1], cons_pert->SD[1]));
    sprintf(fail_msg, "%.90s S[1]", fail_msg);
    test_fail = 1;
  }

  if( ghl_pert_test_fail(cons_trusted->SD[2], cons->SD[2], cons_pert->SD[2]) ) {
    printf("S[2] trusted %.14e computed %.14e perturbed %.14e\n", cons_trusted->SD[2], cons->SD[2], cons_pert->SD[2]);
    printf("rel.err. %.14e %.14e\n", relative_error(cons_trusted->SD[2], cons->SD[2]), relative_error(cons_trusted->SD[2], cons_pert->SD[2]));
    sprintf(fail_msg, "%.90s S[2]", fail_msg);
    test_fail = 1;
  }

  if(evolve_entropy)
    if( ghl_pert_test_fail(cons_trusted->entropy, cons->entropy, cons_pert->entropy) ) {
      printf("entropy trusted %.14e computed %.14e perturbed %.14e\n", cons_trusted->entropy, cons->entropy, cons_pert->entropy);
    printf("rel.err. %.14e %.14e\n", relative_error(cons_trusted->entropy, cons->entropy), relative_error(cons_trusted->entropy, cons_pert->entropy));
      sprintf(fail_msg, "%.90s entropy", fail_msg);
      test_fail = 1;
    }

  if(test_fail) {
    ghl_error("%.100s\n", fail_msg);
  }
  return;
}
