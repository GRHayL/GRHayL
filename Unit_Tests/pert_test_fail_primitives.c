#include "unit_tests.h"

/*
  In the function ghl_pert_test_fail_primitives, the primitives are ghl_pert_test_faild
  individually. In the case of failure, the code errors out and lists all
  the components that failed to be within the error bars of the perturbed
  version of the code.
*/

void ghl_pert_test_fail_primitives(
      const bool evolve_entropy,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims_trusted,
      const ghl_primitive_quantities *restrict prims,
      const ghl_primitive_quantities *restrict prims_pert) {

  char fail_msg[100] = "Test has failed!\n The primitive variable(s) which failed are ";
  int test_fail = 0;
  if( ghl_pert_test_fail(prims_trusted->rho, prims->rho, prims_pert->rho) ) {
    printf("rho_b trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->rho, prims->rho, prims_pert->rho);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->rho, prims->rho), relative_error(prims_trusted->rho, prims_pert->rho));
    sprintf(fail_msg, "%.80s rho_b", fail_msg);
    test_fail = 1;
  }

  if( ghl_pert_test_fail(prims_trusted->press, prims->press, prims_pert->press)) {
    printf("pressure trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->press, prims->press, prims_pert->press);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->press, prims->press), relative_error(prims_trusted->press, prims_pert->press));
    sprintf(fail_msg, "%.80s press", fail_msg);
    test_fail = 1;
  }

  if( ghl_pert_test_fail(prims_trusted->eps, prims->eps, prims_pert->eps)) {
    printf("eps trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->eps, prims->eps, prims_pert->eps);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->eps, prims->eps), relative_error(prims_trusted->eps, prims_pert->eps));
    sprintf(fail_msg, "%.80s eps", fail_msg);
    test_fail = 1;
  }

  if( ghl_pert_test_fail(prims_trusted->vU[0], prims->vU[0], prims_pert->vU[0]) ) {
    printf("vel[0] trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->vU[0], prims->vU[0], prims_pert->vU[0]);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->vU[0], prims->vU[0]), relative_error(prims_trusted->vU[0], prims_pert->vU[0]));
    sprintf(fail_msg, "%.80s vel[0]", fail_msg);
    test_fail = 1;
  }

  if(ghl_pert_test_fail(prims_trusted->vU[1], prims->vU[1], prims_pert->vU[1])) {
    printf("vel[1] trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->vU[1], prims->vU[1], prims_pert->vU[1]);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->vU[1], prims->vU[1]), relative_error(prims_trusted->vU[1], prims_pert->vU[1]));
    sprintf(fail_msg, "%.80s vel[1]", fail_msg);
    test_fail = 1;
  }

  if( ghl_pert_test_fail(prims_trusted->vU[2], prims->vU[2], prims_pert->vU[2]) ) {
    printf("vel[2] trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->vU[2], prims->vU[2], prims_pert->vU[2]);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->vU[2], prims->vU[2]), relative_error(prims_trusted->vU[2], prims_pert->vU[2]));
    sprintf(fail_msg, "%.80s vel[2]", fail_msg);
    test_fail = 1;
  }

  if(evolve_entropy)
    if( ghl_pert_test_fail(prims_trusted->entropy, prims->entropy, prims_pert->entropy) ) {
      printf("entropy trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->entropy, prims->entropy, prims_pert->entropy);
      printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->entropy, prims->entropy), relative_error(prims_trusted->entropy, prims_pert->entropy));
      sprintf(fail_msg, "%.80s entropy", fail_msg);
      test_fail = 1;
    }

  if(eos->eos_type == ghl_eos_tabulated) {
    if( ghl_pert_test_fail(prims_trusted->Y_e, prims->Y_e, prims_pert->Y_e) ) {
      printf("Y_e trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->Y_e, prims->Y_e, prims_pert->Y_e);
      printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->Y_e, prims->Y_e), relative_error(prims_trusted->Y_e, prims_pert->Y_e));
      sprintf(fail_msg, "%.80s Y_e", fail_msg);
      test_fail = 1;
    }
    if( ghl_pert_test_fail(prims_trusted->temperature, prims->temperature, prims_pert->temperature) ) {
      printf("temperature trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->temperature, prims->temperature, prims_pert->temperature);
      printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->temperature, prims->temperature), relative_error(prims_trusted->temperature, prims_pert->temperature));
      sprintf(fail_msg, "%.80s temperature", fail_msg);
      test_fail = 1;
    }
  }

  if(test_fail) {
    ghl_error("%.100s\n", fail_msg);
  }
  return;
}
