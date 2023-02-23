#include "unit_tests.h"

/*
  In the function validate_primitives, the primitives are validated
  individually. In the case of failure, the code errors out and lists all
  the components that failed to be within the error bars of the perturbed
  version of the code.
*/

void validate_primitives(
      const bool evolve_entropy,
      const eos_parameters *restrict eos,
      const primitive_quantities *restrict prims,
      const primitive_quantities *restrict prims_trusted,
      const primitive_quantities *restrict prims_pert) {

  char fail_msg[100] = "Test has failed!\n The primitive variable(s) which failed are ";
  int test_fail = 0;
  if( validate(prims_trusted->rho, prims->rho, prims_pert->rho) ) {
    printf("rho_b trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->rho, prims->rho, prims_pert->rho);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->rho, prims->rho), relative_error(prims_trusted->rho, prims_pert->rho));
    sprintf(fail_msg, "%.80s rho_b", fail_msg);
    test_fail = 1;
  }

  const double pressure_cutoff = 1.0e-18;
  // Pressure has an additional absolute difference check because the pressure can become very small depending on the
  // input values. The pressure coming out of HARM doesn't have the accuracy to preserve the stringent accuracy requirements
  // demanded elsewhere, so this relaxes the demands on the pressure for very small values.
  if( validate(prims_trusted->press, prims->press, prims_pert->press) && fabs(prims_trusted->press-prims->press) > pressure_cutoff ) {
    printf("pressure trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->press, prims->press, prims_pert->press);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->press, prims->press), relative_error(prims_trusted->press, prims_pert->press));
    sprintf(fail_msg, "%.80s press", fail_msg);
    test_fail = 1;
  }

  if( validate(prims_trusted->vx, prims->vx, prims_pert->vx) ) {
    printf("vx trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->vx, prims->vx, prims_pert->vx);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->vx, prims->vx), relative_error(prims_trusted->vx, prims_pert->vx));
    sprintf(fail_msg, "%.80s vx", fail_msg);
    test_fail = 1;
  }

  if(validate(prims_trusted->vy, prims->vy, prims_pert->vy)) {
    printf("vy trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->vy, prims->vy, prims_pert->vy);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->vy, prims->vy), relative_error(prims_trusted->vy, prims_pert->vy));
    sprintf(fail_msg, "%.80s vy", fail_msg);
    test_fail = 1;
  }

  if( validate(prims_trusted->vz, prims->vz, prims_pert->vz) ) {
    printf("vz trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->vz, prims->vz, prims_pert->vz);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->vz, prims->vz), relative_error(prims_trusted->vz, prims_pert->vz));
    sprintf(fail_msg, "%.80s vz", fail_msg);
    test_fail = 1;
  }

  // Epsilon has a similar issue with pressure, so we compute a cutoff that is consistent with the above choice.
  //const double eps_cutoff = pressure_cutoff/(pow(pressure_cutoff/eos->K_ppoly[0], 1.0/eos->Gamma_ppoly[0]) * (eos->Gamma_ppoly[0] - 1.0));
  const double eps_cutoff = 1.0e-11; // Above computed 1e-9, which seemed too large to make sense as a cutoff
  if( validate(prims_trusted->eps, prims->eps, prims_pert->eps) && fabs(prims_trusted->eps-prims->eps) > eps_cutoff ) {
    printf("eps trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->eps, prims->eps, prims_pert->eps);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->eps, prims->eps), relative_error(prims_trusted->eps, prims_pert->eps));
    sprintf(fail_msg, "%.80s eps", fail_msg);
    test_fail = 1;
  }

  if(evolve_entropy)
    if( validate(prims_trusted->entropy, prims->entropy, prims_pert->entropy) ) {
      printf("entropy trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->entropy, prims->entropy, prims_pert->entropy);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->entropy, prims->entropy), relative_error(prims_trusted->entropy, prims_pert->entropy));
      sprintf(fail_msg, "%.80s entropy", fail_msg);
      test_fail = 1;
    }

  if(eos->eos_type == grhayl_eos_tabulated) {
    if( validate(prims_trusted->Y_e, prims->Y_e, prims_pert->Y_e) ) {
      printf("Y_e trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->Y_e, prims->Y_e, prims_pert->Y_e);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->Y_e, prims->Y_e), relative_error(prims_trusted->Y_e, prims_pert->Y_e));
      sprintf(fail_msg, "%.80s Y_e", fail_msg);
      test_fail = 1;
    }
    if( validate(prims_trusted->temperature, prims->temperature, prims_pert->temperature) ) {
      printf("temperature trusted %.14e computed %.14e perturbed %.14e\n", prims_trusted->temperature, prims->temperature, prims_pert->temperature);
    printf("rel.err. %.14e %.14e\n", relative_error(prims_trusted->temperature, prims->temperature), relative_error(prims_trusted->temperature, prims_pert->temperature));
      sprintf(fail_msg, "%.80s temperature", fail_msg);
      test_fail = 1;
    }
  }

  if(test_fail) {
    //grhayl_warn("%.100s\n", fail_msg);
printf(" output %.14e %.14e %.14e %.14e %.14e\n", prims->rho, prims->press, prims->vx, prims->vy, prims->vz);
printf("trusted %.14e %.14e %.14e %.14e %.14e\n", prims_trusted->rho, prims_trusted->press, prims_trusted->vx, prims_trusted->vy, prims_trusted->vz);
printf("perturb %.14e %.14e %.14e %.14e %.14e\n", prims_pert->rho, prims_pert->press, prims_pert->vx, prims_pert->vy, prims_pert->vz);
    grhayl_error("%.100s\n", fail_msg);
  }
  return;
}
