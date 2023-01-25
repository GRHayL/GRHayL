#include "unit_tests.h"

/*
  In the function validate_primitives, the primitives are validated
  individually. In the case of failure, the code errors out and lists all
  the components that failed to be within the error bars of the perturbed
  version of the code.
*/

void validate_primitives(
      const int eos_type,
      const bool evolve_entropy,
      const primitive_quantities *restrict prims,
      const primitive_quantities *restrict prims_trusted,
      const primitive_quantities *restrict prims_pert) {

  char fail_msg[100] = "Test has failed!\n The primitive variable(s) which failed are ";
  int test_fail = 0;
  if( validate(prims_trusted->rho, prims->rho, prims_pert->rho) ) {
    sprintf(fail_msg, "%.80s rho_b", fail_msg);
    printf("%e %e %e\n", prims_trusted->rho, prims->rho, prims_pert->rho);
    test_fail = 1;
  }

  if( validate(prims_trusted->press, prims->press, prims_pert->press) ) {
    sprintf(fail_msg, "%.80s press", fail_msg);
    test_fail = 1;
  }

  if( validate(prims_trusted->vx, prims->vx, prims_pert->vx) ) {
    sprintf(fail_msg, "%.80s vx", fail_msg);
    printf("%e %e %e\n", prims_trusted->vx, prims->vx, prims_pert->vx);
    test_fail = 1;
  }

  if(validate(prims_trusted->vy, prims->vy, prims_pert->vy)) {
    sprintf(fail_msg, "%.80s vy", fail_msg);
    printf("%.15e %.15e %.15e\n", prims_trusted->vy, prims->vy, prims_pert->vy);
    test_fail = 1;
  }

  if( validate(prims_trusted->vz, prims->vz, prims_pert->vz) ) {
    sprintf(fail_msg, "%.80s vz", fail_msg);
    test_fail = 1;
  }

  if(evolve_entropy)
    if( validate(prims_trusted->entropy, prims->entropy, prims_pert->entropy) ) {
      sprintf(fail_msg, "%.80s entropy", fail_msg);
      test_fail = 1;
    }

  if(eos_type == 2) { //Tabulated
    if( validate(prims_trusted->Y_e, prims->Y_e, prims_pert->Y_e) ) {
      sprintf(fail_msg, "%.80s Y_e", fail_msg);
      test_fail = 1;
    }
    if( validate(prims_trusted->temperature, prims->temperature, prims_pert->temperature) ) {
      sprintf(fail_msg, "%.80s temperature", fail_msg);
      test_fail = 1;
    }
  }

  if(test_fail) {
    grhayl_error("%.100s\n", fail_msg);
  }
  return;
}
