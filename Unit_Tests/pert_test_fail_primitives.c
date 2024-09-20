#include "ghl_unit_tests.h"

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

  bool test_fail = false;
  if( ghl_pert_test_fail(prims_trusted->rho, prims->rho, prims_pert->rho) ) {
    printf("rho_b trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->rho, prims->rho, prims_pert->rho,
           relative_error(prims_trusted->rho, prims->rho), relative_error(prims_trusted->rho, prims_pert->rho));
    test_fail = true;
  }

  if( ghl_pert_test_fail(prims_trusted->press, prims->press, prims_pert->press)) {
    printf("pressure trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->press, prims->press, prims_pert->press,
           relative_error(prims_trusted->press, prims->press), relative_error(prims_trusted->press, prims_pert->press));
    test_fail = true;
  }

  if( ghl_pert_test_fail(prims_trusted->eps, prims->eps, prims_pert->eps)) {
    printf("eps trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->eps, prims->eps, prims_pert->eps,
           relative_error(prims_trusted->eps, prims->eps), relative_error(prims_trusted->eps, prims_pert->eps));
    test_fail = true;
  }

  if( ghl_pert_test_fail(prims_trusted->vU[0], prims->vU[0], prims_pert->vU[0]) ) {
    printf("vel[0] trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->vU[0], prims->vU[0], prims_pert->vU[0],
           relative_error(prims_trusted->vU[0], prims->vU[0]), relative_error(prims_trusted->vU[0], prims_pert->vU[0]));
    test_fail = true;
  }

  if(ghl_pert_test_fail(prims_trusted->vU[1], prims->vU[1], prims_pert->vU[1])) {
    printf("vel[1] trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->vU[1], prims->vU[1], prims_pert->vU[1],
           relative_error(prims_trusted->vU[1], prims->vU[1]), relative_error(prims_trusted->vU[1], prims_pert->vU[1]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(prims_trusted->vU[2], prims->vU[2], prims_pert->vU[2]) ) {
    printf("vel[2] trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->vU[2], prims->vU[2], prims_pert->vU[2],
           relative_error(prims_trusted->vU[2], prims->vU[2]), relative_error(prims_trusted->vU[2], prims_pert->vU[2]));
    test_fail = true;
  }

  if(evolve_entropy)
    if( ghl_pert_test_fail(prims_trusted->entropy, prims->entropy, prims_pert->entropy) ) {
      printf("entropy trusted %.14e computed %.14e perturbed %.14e\n"
             "rel.err. %.14e %.14e\n",
             prims_trusted->entropy, prims->entropy, prims_pert->entropy,
             relative_error(prims_trusted->entropy, prims->entropy), relative_error(prims_trusted->entropy, prims_pert->entropy));
      test_fail = true;
    }

  if(eos->eos_type == ghl_eos_tabulated) {
    if( ghl_pert_test_fail(prims_trusted->Y_e, prims->Y_e, prims_pert->Y_e) ) {
      printf("Y_e trusted %.14e computed %.14e perturbed %.14e\n"
             "rel.err. %.14e %.14e\n",
             prims_trusted->Y_e, prims->Y_e, prims_pert->Y_e,
             relative_error(prims_trusted->Y_e, prims->Y_e), relative_error(prims_trusted->Y_e, prims_pert->Y_e));
      test_fail = true;
    }
    if( ghl_pert_test_fail(prims_trusted->temperature, prims->temperature, prims_pert->temperature) ) {
      printf("temperature trusted %.14e computed %.14e perturbed %.14e\n"
             "rel.err. %.14e %.14e\n",
             prims_trusted->temperature, prims->temperature, prims_pert->temperature,
             relative_error(prims_trusted->temperature, prims->temperature), relative_error(prims_trusted->temperature, prims_pert->temperature));
      test_fail = true;
    }
  }

  if(test_fail)
    ghl_error("Test has failed!\n");
}

void ghl_pert_test_fail_primitives_with_cutoffs(
      const bool evolve_entropy,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims_trusted,
      const ghl_primitive_quantities *restrict prims,
      const ghl_primitive_quantities *restrict prims_pert,
      const double pressure_cutoff,
      const double eps_cutoff) {

  bool test_fail = false;
  if( ghl_pert_test_fail(prims_trusted->rho, prims->rho, prims_pert->rho) ) {
    printf("rho_b trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->rho, prims->rho, prims_pert->rho,
           relative_error(prims_trusted->rho, prims->rho), relative_error(prims_trusted->rho, prims_pert->rho));
    test_fail = true;
  }

  const double min_rel = 8.0e-14; // This is the default relative tolerance cutoff used by ghl_pert_test_fail()
  if( ghl_pert_test_fail_with_tolerance(prims_trusted->press, prims->press, prims_pert->press, min_rel, pressure_cutoff)) {
    printf("pressure trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->press, prims->press, prims_pert->press,
           relative_error(prims_trusted->press, prims->press), relative_error(prims_trusted->press, prims_pert->press));
    test_fail = true;
  }

  if( ghl_pert_test_fail_with_tolerance(prims_trusted->eps, prims->eps, prims_pert->eps, min_rel, eps_cutoff)) {
    printf("eps trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->eps, prims->eps, prims_pert->eps,
           relative_error(prims_trusted->eps, prims->eps), relative_error(prims_trusted->eps, prims_pert->eps));
    test_fail = true;
  }

  if( ghl_pert_test_fail(prims_trusted->vU[0], prims->vU[0], prims_pert->vU[0]) ) {
    printf("vel[0] trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
            prims_trusted->vU[0], prims->vU[0], prims_pert->vU[0],
           relative_error(prims_trusted->vU[0], prims->vU[0]), relative_error(prims_trusted->vU[0], prims_pert->vU[0]));
    test_fail = true;
  }

  if(ghl_pert_test_fail(prims_trusted->vU[1], prims->vU[1], prims_pert->vU[1])) {
    printf("vel[1] trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->vU[1], prims->vU[1], prims_pert->vU[1],
           relative_error(prims_trusted->vU[1], prims->vU[1]), relative_error(prims_trusted->vU[1], prims_pert->vU[1]));
    test_fail = true;
  }

  if( ghl_pert_test_fail(prims_trusted->vU[2], prims->vU[2], prims_pert->vU[2]) ) {
    printf("vel[2] trusted %.14e computed %.14e perturbed %.14e\n"
           "rel.err. %.14e %.14e\n",
           prims_trusted->vU[2], prims->vU[2], prims_pert->vU[2],
           relative_error(prims_trusted->vU[2], prims->vU[2]), relative_error(prims_trusted->vU[2], prims_pert->vU[2]));
    test_fail = true;
  }

  if(evolve_entropy)
    if( ghl_pert_test_fail(prims_trusted->entropy, prims->entropy, prims_pert->entropy) ) {
      printf("entropy trusted %.14e computed %.14e perturbed %.14e\n"
             "rel.err. %.14e %.14e\n",
             prims_trusted->entropy, prims->entropy, prims_pert->entropy,
             relative_error(prims_trusted->entropy, prims->entropy), relative_error(prims_trusted->entropy, prims_pert->entropy));
      test_fail = true;
    }

  if(eos->eos_type == ghl_eos_tabulated) {
    if( ghl_pert_test_fail(prims_trusted->Y_e, prims->Y_e, prims_pert->Y_e) ) {
      printf("Y_e trusted %.14e computed %.14e perturbed %.14e\n"
             "rel.err. %.14e %.14e\n",
             prims_trusted->Y_e, prims->Y_e, prims_pert->Y_e,
             relative_error(prims_trusted->Y_e, prims->Y_e), relative_error(prims_trusted->Y_e, prims_pert->Y_e));
      test_fail = true;
    }
    if( ghl_pert_test_fail(prims_trusted->temperature, prims->temperature, prims_pert->temperature) ) {
      printf("temperature trusted %.14e computed %.14e perturbed %.14e\n"
             "rel.err. %.14e %.14e\n",
             prims_trusted->temperature, prims->temperature, prims_pert->temperature,
             relative_error(prims_trusted->temperature, prims->temperature), relative_error(prims_trusted->temperature, prims_pert->temperature));
      test_fail = true;
    }
  }

  if(test_fail)
    ghl_error("Test has failed!\n");
}
