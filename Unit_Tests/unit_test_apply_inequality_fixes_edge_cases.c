#include "unit_tests.h"

int main(int argc, char **argv) {

  FILE* infile = fopen("apply_inequality_fixes_edge_cases_input.bin", "rb");
  check_file_was_successfully_open(infile, "unit_test_apply_inequality_fixes_edge_cases_input.bin");

  int ntest; // Number of cases tested
  int key = fread(&ntest, sizeof(int), 1, infile);

  // Only need the magnetic field primitives for this test
  double Bx[ntest];
  double By[ntest];
  double Bz[ntest];
  double rho_star[ntest];
  double tau[ntest], tau_trusted[ntest], tau_pert[ntest];
  double S_x[ntest], S_x_trusted[ntest], S_x_pert[ntest];
  double S_y[ntest], S_y_trusted[ntest], S_y_pert[ntest];
  double S_z[ntest], S_z_trusted[ntest], S_z_pert[ntest];

  key  = fread(Bx,       sizeof(double), ntest, infile);
  key += fread(By,       sizeof(double), ntest, infile);
  key += fread(Bz,       sizeof(double), ntest, infile);
  key += fread(rho_star, sizeof(double), ntest, infile);
  key += fread(tau,      sizeof(double), ntest, infile);
  key += fread(S_x,      sizeof(double), ntest, infile);
  key += fread(S_y,      sizeof(double), ntest, infile);
  key += fread(S_z,      sizeof(double), ntest, infile);
  fclose(infile);
  if(key != ntest*8)
    grhayl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {None,None,None};
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;
  const int update_Tmunu = 1;

  const int neos = 1;
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during a simulation.
  GRHayL_parameters params;
  initialize_GRHayL(None, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess,
                    Psi6threshold, update_Tmunu, 1 /*Cupp Fix*/, 0.0 /*Lorenz damping factor*/, &params);

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

  // Define the various GRHayL structs for the unit tests
  con2prim_diagnostics diagnostics;
  initialize_diagnostics(&diagnostics);
  metric_quantities metric;
  primitive_quantities prims;
  conservative_quantities cons;

  initialize_metric(1.0,
                    1.0, 0.0, 0.0,
                    1.0, 0.0, 1.0,
                    0.0, 0.0, 0.0,
                    &metric);

  // Read initial data and accompanying trusted output
  FILE* outfile = fopen("apply_inequality_fixes_edge_cases_output.bin", "rb");
  check_file_was_successfully_open(outfile, "unit_test_apply_inequality_fixes_edge_cases_output.bin");
  key  = fread(tau_trusted, sizeof(double), ntest, outfile);
  key += fread(S_x_trusted, sizeof(double), ntest, outfile);
  key += fread(S_y_trusted, sizeof(double), ntest, outfile);
  key += fread(S_z_trusted, sizeof(double), ntest, outfile);
  fclose(outfile);
  if(key != ntest*4)
    grhayl_error("An error has occured with reading in output data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  outfile = fopen("apply_inequality_fixes_edge_cases_output_pert.bin", "rb");
  check_file_was_successfully_open(outfile, "unit_test_apply_inequality_fixes_edge_cases_output_pert.bin");
  key  = fread(tau_pert, sizeof(double), ntest, outfile);
  key += fread(S_x_pert, sizeof(double), ntest, outfile);
  key += fread(S_y_pert, sizeof(double), ntest, outfile);
  key += fread(S_z_pert, sizeof(double), ntest, outfile);
  fclose(outfile);
  if(key != ntest*4)
    grhayl_error("An error has occured with reading in perturbed data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Edge case 1) small B
  //   Subcase a) S fix 1, Bx < By, Bz (also catches edge case of Bbar < 1e-300)
  //   Subcase b) Don't fix S, Bx > By, Bz
  for(int i=0; i<2; i++) {
    initialize_primitives(
                      0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      Bx[i], By[i], Bz[i],
                      0.0, 0.0, 0.0,
                      &prims);

    initialize_conservatives(rho_star[i], tau[i],
                             S_x[i], S_y[i], S_z[i],
                             0.0, 0.0, &cons);

    apply_inequality_fixes(&params, &eos, &metric, &prims, &cons, &diagnostics);

    if( validate(tau_trusted[i], cons.tau, tau_pert[i]) )
      grhayl_error("Test unit_test_apply_inequality_fixes_edge_cases has failed for variable tau.\n"
                   "  tau trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", tau_trusted[i], cons.tau, tau_pert[i],
                                               relative_error(tau_trusted[i], cons.tau),
                                               relative_error(tau_trusted[i], tau_pert[i]));

    if( validate(S_x_trusted[i], cons.S_x, S_x_pert[i]) )
      grhayl_error("Test unit_test_apply_inequality_fixes_edge_cases has failed for variable S_x.\n"
                   "  S_x trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", S_x_trusted[i], cons.S_x, S_x_pert[i],
                                               relative_error(S_x_trusted[i], cons.S_x),
                                               relative_error(S_x_trusted[i], S_x_pert[i]));

    if( validate(S_y_trusted[i], cons.S_y, S_y_pert[i]) )
      grhayl_error("Test unit_test_apply_inequality_fixes_edge_cases has failed for variable S_y.\n"
                   "  S_y trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", S_y_trusted[i], cons.S_y, S_y_pert[i],
                                               relative_error(S_y_trusted[i], cons.S_y),
                                               relative_error(S_y_trusted[i], S_y_pert[i]));

    if( validate(S_z_trusted[i], cons.S_z, S_z_pert[i]) )
      grhayl_error("Test unit_test_apply_inequality_fixes_edge_cases has failed for variable S_z.\n"
                   "  S_z trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", S_z_trusted[i], cons.S_z, S_z_pert[i],
                                               relative_error(S_z_trusted[i], cons.S_z),
                                               relative_error(S_z_trusted[i], S_z_pert[i]));
  }

  // Edge case 2) \psi^6 too large
  //   Subcase a) \tau fix 1 & 2, don't fix S
  //   Subcase b) Don't fix \tau, S fix 2

  // Artificially creating case by lowering threshold
  params.psi6threshold = 1e-1;

  for(int i=2; i<4; i++) {
    initialize_primitives(
                      0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      Bx[i], By[i], Bz[i],
                      0.0, 0.0, 0.0,
                      &prims);

    initialize_conservatives(rho_star[i], tau[i],
                             S_x[i], S_y[i], S_z[i],
                             0.0, 0.0, &cons);

    apply_inequality_fixes(&params, &eos, &metric, &prims, &cons, &diagnostics);

    if( validate(tau_trusted[i], cons.tau, tau_pert[i]) )
      grhayl_error("Test unit_test_apply_inequality_fixes_edge_cases has failed for variable tau.\n"
                   "  tau trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", tau_trusted[i], cons.tau, tau_pert[i],
                                               relative_error(tau_trusted[i], cons.tau),
                                               relative_error(tau_trusted[i], tau_pert[i]));

    if( validate(S_x_trusted[i], cons.S_x, S_x_pert[i]) )
      grhayl_error("Test unit_test_apply_inequality_fixes_edge_cases has failed for variable S_x.\n"
                   "  S_x trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", S_x_trusted[i], cons.S_x, S_x_pert[i],
                                               relative_error(S_x_trusted[i], cons.S_x),
                                               relative_error(S_x_trusted[i], S_x_pert[i]));

    if( validate(S_y_trusted[i], cons.S_y, S_y_pert[i]) )
      grhayl_error("Test unit_test_apply_inequality_fixes_edge_cases has failed for variable S_y.\n"
                   "  S_y trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", S_y_trusted[i], cons.S_y, S_y_pert[i],
                                               relative_error(S_y_trusted[i], cons.S_y),
                                               relative_error(S_y_trusted[i], S_y_pert[i]));

    if( validate(S_z_trusted[i], cons.S_z, S_z_pert[i]) )
      grhayl_error("Test unit_test_apply_inequality_fixes_edge_cases has failed for variable S_z.\n"
                   "  S_z trusted %.14e computed %.14e perturbed %.14e\n"
                   "  rel.err. %.14e %.14e\n", S_z_trusted[i], cons.S_z, S_z_pert[i],
                                               relative_error(S_z_trusted[i], cons.S_z),
                                               relative_error(S_z_trusted[i], S_z_pert[i]));
  }
  grhayl_info("apply_inequality_fixes_edge_cases test has passed!\n");
}
