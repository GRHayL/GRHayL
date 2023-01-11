// Thorn      : GRHayL
// File       : con2prim_unit_test.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.
#include "unit_tests.h"

#define check_file_was_successfully_open(fp, filename) \
  if( fp == NULL ) { \
    fprintf(stderr, "(GRHayL) ERROR: Could not open file %s. Terminating.\n", filename); \
    exit(1); \
  }

inline void perturb_data(double *restrict rand_val, primitive_quantities *restrict prims, conservative_quantities *restrict cons) {
  prims->rho   *= rand_val[0];
  prims->press *= rand_val[1];
  prims->vx    *= rand_val[2];
  prims->vy    *= rand_val[3];
  prims->vz    *= rand_val[4];
  prims->Bx    *= rand_val[5];
  prims->By    *= rand_val[6];
  prims->Bz    *= rand_val[7];
  cons->rho    *= rand_val[8];
  cons->S_x    *= rand_val[9];
  cons->S_y    *= rand_val[10];
  cons->S_z    *= rand_val[11];
  cons->tau    *= rand_val[12];
}

int main(int argc, char **argv) {

  // These variables set up the tested range of values
  // and number of sampling points.
  int npoints = 256; //Number of sampling points in density and temperature
  double test_rho_min = 1e-12; //Minimum input density
  double test_rho_max = 1e-3; //Maximum input density
  // double test_T_min = 1e-2; //Minimum input temperature
  // double test_T_max = 1e+2; //Maximum input temperature

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
  double rho_b_max = 1e300; //IGM default
  double gamma_th = 2.0; //Taken from magnetizedTOV.par
  double rho_tab[1] = {0.0};
  double gamma_tab[1] = {2.0};
  double k_tab = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  GRHayL_parameters params;
  initialize_GRHayL(None, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess, Psi6threshold, update_Tmunu, 1 /*Cupp Fix*/, &params);

  eos_parameters eos;
  initialize_general_eos(eos_type, W_max,
             test_rho_min, test_rho_min, rho_b_max,
             &eos);

  initialize_hybrid_functions(&eos);

  initialize_hybrid_eos(neos, rho_tab,
             gamma_tab, k_tab, gamma_th,
             &eos);

  con2prim_diagnostics diagnostics;
  initialize_diagnostics(&diagnostics);

  // We will be performing the tabulated EOS test in the following way:
  //
  //      Y_e      = 0.1
  //       W       = 2
  // log10(Pmag/P) = -5
  //
  // rho will vary between rho_min and rho_max (uniformly in log space)
  //  T  will vary between  T_min  and  T_max  (uniformly in log space)

  // Compute the density step size
  const double lrmin        = log(test_rho_min);
  const double lrmax        = log(test_rho_max);
  const double dlr          = (lrmax - lrmin)/(npoints-1);

  // Compute the temperature step size
  //const double ltmin        = log(test_T_min);
  //const double ltmax        = log(test_T_max);
  //const double dlt          = (ltmax - ltmin)/(npoints-1);

  // tau is given by (see :
  //
  // tau := hW^{2} + B^{2} - P - 0.5*( (B.v)^{2} + (B/W)^{2} )
  // Absolutely minimum allowed tau

  char filename[100];
  FILE* summaryf;
  sprintf(filename,"unit_test/C2P_Summary.asc");
  summaryf = fopen(filename,"w");

  // Now perform one test for each of the selected routines
  for(int which_routine=0;which_routine<num_routines_tested;which_routine++) {
    params.main_routine = con2prim_test_keys[which_routine];

    int failures = 0;
    for(int rand=0;rand<2;rand++) {

      double rand_val[13];
      char suffix[10] = "norm";
      if(rand==1) {
        srand(1000000);
        sprintf(suffix, "rand");
        for(int i=0;i<13;i++) rand_val[i] = 1.0 + randf(-1,1)*1.0e-14;
      }

      printf("Beginning %s test for routine %s\n", suffix, con2prim_test_names[which_routine]);

      FILE* outfiles[7];

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_limit_v_and_output_u0.bin",con2prim_test_names[which_routine], suffix);
      outfiles[0] = fopen(filename,"wb");
      check_file_was_successfully_open(outfiles[0], filename);

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_apply_inequality_fixes.bin",con2prim_test_names[which_routine], suffix);
      outfiles[1] = fopen(filename,"w");
      check_file_was_successfully_open(outfiles[1], filename);

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_Hybrid_Multi_Method.bin",con2prim_test_names[which_routine], suffix);
      outfiles[2] = fopen(filename,"wb");
      check_file_was_successfully_open(outfiles[2], filename);

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_font_fix.bin",con2prim_test_names[which_routine], suffix);
      outfiles[3] = fopen(filename,"wb");
      check_file_was_successfully_open(outfiles[3], filename);

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_enforce_primitive_limits_and_output_u0.bin",con2prim_test_names[which_routine], suffix);
      outfiles[4] = fopen(filename,"wb");
      check_file_was_successfully_open(outfiles[4], filename);

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_compute_conservs.bin",con2prim_test_names[which_routine], suffix);
      outfiles[5] = fopen(filename,"w");
      check_file_was_successfully_open(outfiles[5], filename);

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_compute_Tmunu.bin",con2prim_test_names[which_routine], suffix);
      outfiles[6] = fopen(filename,"w");
      check_file_was_successfully_open(outfiles[6], filename);

      srand(0);

      for(int i=0;i<npoints;i++) { // Density loop
        double xrho  = exp(lrmin + dlr*i);
        double P_cold = 0.0;
        eos.hybrid_compute_P_cold(&eos, xrho, &P_cold);

        // Compute the pressure step size
        const double lpmin        = log(1.0e-30);//-P_cold);
        const double lpmax        = log(10.0*P_cold);
        const double dlp          = (lpmax - lpmin)/(npoints-1);
        for(int j=0;j<npoints;j++) { // Pressure loop
          // Start by setting the prims (rho,Ye,T,P,eps)
          //double xtemp = exp(ltmin + dlt*j);
          //double xye   = Ye_test;
          double xpress  = exp(lpmin + dlp*j);
          //WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xpress,&xeps );

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

          // Generate random data to serve as the 'true' primitive values
          bool random_metric = true;
          initial_random_data(xrho, xpress, random_metric, &metric, &prims);

          double u0 = poison;
          prims_orig = prims;
          limit_v_and_output_u0(&eos, &metric, &prims, &u0, &diagnostics);
          write_primitive_binary(eos.eos_type, 1, params.evolve_entropy, &prims_orig, &prims, outfiles[0]);

          // Compute conservatives based on these primitives
          compute_conservs_and_Tmunu(&params, &eos, &metric, &prims, u0, &cons, &Tmunu);

          //This is meant to simulate some round-off error that deviates from the "true" values that we just computed.
          if(rand) perturb_data(rand_val, &prims, &cons);

          int check = 0;
          if(cons.rho > 0.0) {

            //This applies the inequality (or "Faber") fixes on the conservatives
            if(eos.eos_type == 0) { //Hybrid-only
              cons_orig = cons;
              apply_inequality_fixes(&params, &eos, &metric, &prims, &cons, &diagnostics);
              write_conservative_binary(params.evolve_entropy, &cons_orig, &cons, outfiles[1]);
            }

            // The Con2Prim routines require the undensitized variables, but IGM evolves the densitized variables.
            undensitize_conservatives(&metric, &cons, &cons_undens);

            /************* Conservative-to-primitive recovery ************/

            prims_orig = prims;
            check = Hybrid_Multi_Method(&params, &eos, &metric, &cons_undens, &prims, &prims_guess, &diagnostics);
            write_primitive_binary(eos.eos_type, 0, params.evolve_entropy, &prims_orig, &prims_guess, outfiles[2]);

            if(check!=0) {
              check = font_fix(&eos, &metric, &cons, &prims, &prims_guess, &diagnostics);
              diagnostics.font_fixes++;
              write_primitive_binary(eos.eos_type, 0, params.evolve_entropy, &prims_orig, &prims_guess, outfiles[3]);
            } else { //The else is so that Font Fix is tested even when the primary routine succeeds.
              primitive_quantities prims_tmp;
              check = font_fix(&eos, &metric, &cons, &prims, &prims_tmp, &diagnostics);
              write_primitive_binary(eos.eos_type, 0, params.evolve_entropy, &prims_orig, &prims_tmp, outfiles[3]);
            }

            /*************************************************************/

            if(check==0) {
              prims = prims_guess;
            } else {
              printf("Con2Prim and Font fix failed!");
              printf("diagnostics->failure_checker = %d st_i = %e %e %e, rhostar = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                      diagnostics.failure_checker, cons_orig.S_x, cons_orig.S_y, cons_orig.S_z, cons_orig.rho, prims.Bx, prims.By, prims.Bz,
                      metric.adm_gxx, metric.adm_gxy, metric.adm_gxz, metric.adm_gyy, metric.adm_gyz, metric.adm_gzz, metric.psi6);
            }
          } else {
            diagnostics.failure_checker+=1;
            reset_prims_to_atmosphere(&params, &eos, &metric, &prims, &diagnostics);
            //TODO: Validate reset? (rhob press v)
            printf("Negative rho_* triggering atmospheric reset.\n");
            diagnostics.rho_star_fix_applied++;
          } // if rho_star > 0

          //--------------------------------------------------
          //---------- Primitive recovery completed ----------
          //--------------------------------------------------
          // Enforce limits on primitive variables and recompute conservatives.
          prims_orig = prims;
          enforce_primitive_limits_and_output_u0(&params, &eos, &metric, &prims, &u0, &diagnostics);
          write_primitive_binary(eos.eos_type, 0, params.evolve_entropy, &prims_orig, &prims, outfiles[4]);

          cons_orig = cons;
          Tmunu_orig = Tmunu;
          compute_conservs_and_Tmunu(&params, &eos, &metric, &prims, u0, &cons, &Tmunu);
          write_conservative_binary(params.evolve_entropy, &cons_orig, &cons, outfiles[5]);
          write_stress_energy_binary(&Tmunu_orig, &Tmunu, outfiles[6]);

          if( check != 0 ) {
            failures++;
            printf("Recovery FAILED!\n");
          } else {
            printf("Recovery SUCCEEDED FOR POINT (%d,%d)!\n", i,j);
          }

        } // Pressure loop
      } // Density loop
      for(int k = 0; k < (sizeof(outfiles)/sizeof(outfiles[0])); k++) fclose(outfiles[k]);
    } // perturbation loop

    int ntotal = npoints*npoints;

    printf("Completed test for routine %s\n",con2prim_test_names[which_routine]);
    printf("Final report:\n");
    printf("    Number of recovery attempts: %d\n",ntotal);
    printf("    Number of failed recoveries: %d\n",failures);
    printf("    Recovery failure rate      : %.2lf%%\n",((double)failures)/((double)ntotal)*100.0);

    fprintf(summaryf, "Completed test for routine %s\n",con2prim_test_names[which_routine]);
    fprintf(summaryf, "Final report:\n");
    fprintf(summaryf, "    Number of recovery attempts: %d\n",ntotal);
    fprintf(summaryf, "    Number of failed recoveries: %d\n",failures);
    fprintf(summaryf, "    Recovery failure rate      : %.2lf%%\n",((double)failures)/((double)ntotal)*100.0);
  }
  fprintf(summaryf, "All done! Terminating the run.\n");
  fclose(summaryf);

  return 0;
}
