// Thorn      : GRHayL
// File       : unit_test_con2prim.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.
#include "unit_tests.h"

// Tolerance limit for numerical values
const double relative_tolerance = 1.0e-13;

int main(int argc, char **argv) {

  //Number of sampling points in density and pressure/temperature
  int npoints = 80;

  // Count number of routines tested
  int num_routines_tested = 1;
  int con2prim_test_keys[num_routines_tested];
  char con2prim_test_names[num_routines_tested][50];

  con2prim_test_keys[0] = Noble2D;
  sprintf(con2prim_test_names[0],"%s","Noble2D");

  double poison = 1e200;
  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  int backup_routine[3] = {None,None,None};
  bool calc_prims_guess = true;
  double Psi6threshold = 1e100; //Taken from magnetizedTOV.par
  int update_Tmunu = 1; //IGM default

  int eos_type = 0; // Hybrid=0, Tabulated=1;
  int neos = 1;
  double W_max = 10.0; //IGM default
  double rho_b_min = 1e-12;
  double rho_b_max = 1e300; //IGM default
  double gamma_th = 2.0; //Taken from magnetizedTOV.par
  double rho_ppoly[1] = {0.0};
  double gamma_ppoly[1] = {2.0};
  double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  GRHayL_parameters params;
  initialize_GRHayL(None, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess, Psi6threshold, update_Tmunu, 1 /*Cupp Fix*/, &params);

  eos_parameters eos;
  initialize_general_eos(eos_type, W_max,
             rho_b_min, rho_b_min, rho_b_max,
             &eos);

  initialize_hybrid_functions(&eos);

  initialize_hybrid_eos(neos, rho_ppoly,
             gamma_ppoly, k_ppoly0, gamma_th,
             &eos);

  con2prim_diagnostics diagnostics;
  initialize_diagnostics(&diagnostics);

  char filename[100];
  // Now perform one test for each of the selected routines
  // TODO: unlike main loop, this loop can be parallelized because each routine uses a different file
  for(int which_routine=0;which_routine<num_routines_tested;which_routine++) {
    params.main_routine = con2prim_test_keys[which_routine];

    int failures = 0;

    printf("Beginning test for routine %s\n", con2prim_test_names[which_routine]);

    FILE* initial_data;
    sprintf(filename,"C2P_%.30s_norm_initial_data.bin",con2prim_test_names[which_routine]);
    initial_data = fopen(filename,"rb");
    check_file_was_successfully_open(initial_data, filename);

    FILE* infiles[7];

    sprintf(filename,"C2P_%.30s_norm_limit_v_and_output_u0.bin",con2prim_test_names[which_routine]);
    infiles[0] = fopen(filename,"rb");
    check_file_was_successfully_open(infiles[0], filename);

    sprintf(filename,"C2P_%.30s_norm_apply_inequality_fixes.bin",con2prim_test_names[which_routine]);
    infiles[1] = fopen(filename,"rb");
    check_file_was_successfully_open(infiles[1], filename);

    sprintf(filename,"C2P_%.30s_norm_Hybrid_Multi_Method.bin",con2prim_test_names[which_routine]);
    infiles[2] = fopen(filename,"rb");
    check_file_was_successfully_open(infiles[2], filename);

    sprintf(filename,"C2P_%.30s_norm_font_fix.bin",con2prim_test_names[which_routine]);
    infiles[3] = fopen(filename,"rb");
    check_file_was_successfully_open(infiles[3], filename);

    sprintf(filename,"C2P_%.30s_norm_enforce_primitive_limits_and_output_u0.bin",con2prim_test_names[which_routine]);
    infiles[4] = fopen(filename,"rb");
    check_file_was_successfully_open(infiles[4], filename);

    sprintf(filename,"C2P_%.30s_norm_compute_conservs.bin",con2prim_test_names[which_routine]);
    infiles[5] = fopen(filename,"rb");
    check_file_was_successfully_open(infiles[5], filename);

    sprintf(filename,"C2P_%.30s_norm_compute_Tmunu.bin",con2prim_test_names[which_routine]);
    infiles[6] = fopen(filename,"rb");
    check_file_was_successfully_open(infiles[6], filename);

    srand(0);

    // Can't parallelize because it could change the behavior of reading from the files
    for(int i=0;i<npoints;i++) { // Density loop
      for(int j=0;j<npoints;j++) { // Pressure loop
        // Define the various GRHayL structs for the unit tests
        metric_quantities metric;
        primitive_quantities prims, prims_orig, prims_guess;
        conservative_quantities cons, cons_orig, cons_undens;
        stress_energy Tmunu, Tmunu_orig;

        // Leo says: Initialize cons_orig to silence warning; remove later.
        cons_orig.rho = 1e300;
        cons_orig.S_x = 1e300;
        cons_orig.S_y = 1e300;
        cons_orig.S_z = 1e300;

        // Read initial data accompanying trusted output
        read_metric_binary(&metric, initial_data);
        read_primitive_binary(eos.eos_type, 0, params.evolve_entropy, &prims, initial_data);

        double u0 = poison;
        int test_fail;
        prims_orig = prims;
        limit_v_and_output_u0(&eos, &metric, &prims, &u0, &diagnostics);
        test_fail = validate_primitives(relative_tolerance, eos.eos_type, 1, params.evolve_entropy, &prims_orig, &prims, infiles[0]);
        if(test_fail) {
          printf("Test unit_test_con2prim has failed with error code %d after function limit_v_and_output_u0.\n"
                 "Please check file Unit_Tests/validate_primitives.c for information on the possible exit codes.\n", test_fail);
          exit(1);
        }

        // Compute conservatives based on these primitives
        compute_conservs_and_Tmunu(&params, &eos, &metric, &prims, u0, &cons, &Tmunu);

        int check = 0;
        if(cons.rho > 0.0) {

          //This applies the inequality (or "Faber") fixes on the conservatives
          if(eos.eos_type == 0) { //Hybrid-only
            cons_orig = cons;
            apply_inequality_fixes(&params, &eos, &metric, &prims, &cons, &diagnostics);
            test_fail = validate_conservatives(relative_tolerance, params.evolve_entropy, &cons_orig, &cons, infiles[1]);
            if(test_fail) {
              printf("Test unit_test_con2prim has failed with error code %d after function apply_inequality_fixes.\n"
                     "Please check file Unit_Tests/validate_conservatives.c for information on the possible exit codes.\n", test_fail);
              exit(1);
            }
          }

          // The Con2Prim routines require the undensitized variables, but IGM evolves the densitized variables.
          undensitize_conservatives(&metric, &cons, &cons_undens);

          /************* Conservative-to-primitive recovery ************/
          prims_orig = prims;
          check = Hybrid_Multi_Method(&params, &eos, &metric, &cons_undens, &prims, &prims_guess, &diagnostics);
          test_fail = validate_primitives(relative_tolerance, eos.eos_type, 0, params.evolve_entropy, &prims_orig, &prims_guess, infiles[2]);
          if(test_fail) {
            printf("Test unit_test_con2prim has failed with error code %d after function Hybrid_Multi_Method.\n"
                   "Please check file Unit_Tests/validate_primitives.c for information on the possible exit codes.\n", test_fail);
            exit(1);
          }

          if(check!=0) {
            check = font_fix(&eos, &metric, &cons, &prims, &prims_guess, &diagnostics);
            test_fail = validate_primitives(relative_tolerance, eos.eos_type, 0, params.evolve_entropy, &prims_orig, &prims_guess, infiles[3]);
          } else { //The else is so that Font Fix is tested even when the primary routine succeeds.
            primitive_quantities prims_tmp;
            check = font_fix(&eos, &metric, &cons, &prims, &prims_tmp, &diagnostics);
            test_fail = validate_primitives(relative_tolerance, eos.eos_type, 0, params.evolve_entropy, &prims_orig, &prims_tmp, infiles[3]);
          }
          if(test_fail) {
            printf("Test unit_test_con2prim has failed with error code %d after function font_fix.\n"
                   "Please check file Unit_Tests/validate_primitives.c for information on the possible exit codes.\n", test_fail);
            exit(1);
          }

          /*************************************************************/

          if(check==0) {
            prims = prims_guess;
          } else {
            printf("Con2Prim and Font fix failed!");
          }
        } else {
          diagnostics.failure_checker+=1;
          reset_prims_to_atmosphere(&params, &eos, &metric, &prims, &diagnostics);
          //TODO: Validate reset? (rhob press v)
          printf("Negative rho_* triggering atmospheric reset.\n");
        } // if rho_star > 0

        //--------------------------------------------------
        //---------- Primitive recovery completed ----------
        //--------------------------------------------------
        // Enforce limits on primitive variables and recompute conservatives.
        prims_orig = prims;
        enforce_primitive_limits_and_output_u0(&params, &eos, &metric, &prims, &u0, &diagnostics);
        validate_primitives(relative_tolerance, eos.eos_type, 0, params.evolve_entropy, &prims_orig, &prims, infiles[4]);
          if(test_fail) {
            printf("Test unit_test_con2prim has failed with error code %d after function enforce_primitive_limits_and_output_u0.\n"
                   "Please check file Unit_Tests/validate_primitives.c for information on the possible exit codes.\n", test_fail);
            exit(1);
          }

        cons_orig = cons;
        Tmunu_orig = Tmunu;
        compute_conservs_and_Tmunu(&params, &eos, &metric, &prims, u0, &cons, &Tmunu);
        test_fail = validate_conservatives(relative_tolerance, params.evolve_entropy, &cons_orig, &cons, infiles[5]);
        if(test_fail) {
          printf("Test unit_test_con2prim has failed with error code %d after function compute_conservs_and_Tmunu.\n"
                 "Please check file Unit_Tests/validate_conservatives.c for information on the possible exit codes.\n", test_fail);
          exit(1);
        }
        test_fail = validate_stress_energy(relative_tolerance, &Tmunu_orig, &Tmunu, infiles[6]);
        if(test_fail) {
          printf("Test unit_test_con2prim has failed with error code %d after function compute_conservs_and_Tmunu.\n"
                 "Please check file Unit_Tests/validate_stress_energy.c for information on the possible exit codes.\n", test_fail);
          exit(1);
        }

        if( check != 0 )
          failures++;
      } // Pressure loop
    } // Density loop
    for(int k = 0; k < (sizeof(infiles)/sizeof(infiles[0])); k++) fclose(infiles[k]);

    int ntotal = npoints*npoints;

    printf("Completed test for routine %s\n",con2prim_test_names[which_routine]);
    printf("Final report:\n");
    printf("    Number of recovery attempts: %d\n",ntotal);
    printf("    Number of failed recoveries: %d\n",failures);
    printf("    Recovery failure rate      : %.2lf%%\n",((double)failures)/((double)ntotal)*100.0);
  }
  return 0;
}
