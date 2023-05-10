// Thorn      : GRHayL
// File       : unit_test_data_con2prim.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.
#include "unit_tests.h"

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
  int npoints = 80; //Number of sampling points in density and temperature
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

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  int backup_routine[3] = {None,None,None};
  bool calc_prims_guess = true;
  double Psi6threshold = 1e100; //Taken from magnetizedTOV.par

  int neos = 1;
  double W_max = 10.0; //IGM default
  double rho_b_min = 1e-12;
  double rho_b_max = 1e300; //IGM default
  double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  double rho_ppoly[1] = {0.0};
  double Gamma_ppoly[1] = {2.0};
  double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  GRHayL_parameters params;
  grhayl_initialize(None, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess, Psi6threshold, 1 /*Cupp Fix*/, 0, &params);

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

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

  const double poison = 0.0/0.0;
  // Compute the temperature step size
  //const double ltmin        = log(test_T_min);
  //const double ltmax        = log(test_T_max);
  //const double dlt          = (ltmax - ltmin)/(npoints-1);

  // tau is given by (see :
  //
  // tau := hW^{2} + B^{2} - P - 0.5*( (B.v)^{2} + (B/W)^{2} )
  // Absolutely minimum allowed tau

  char filename[100];
  // Now perform one test for each of the selected routines
  // TODO: unlike main loop, this loop can be parallelized because each routine uses a different file
  for(int which_routine=0;which_routine<num_routines_tested;which_routine++) {
    params.main_routine = con2prim_test_keys[which_routine];
    FILE* initial_data;
    sprintf(filename,"%.30s_initial_data.bin", con2prim_test_names[which_routine]);
    initial_data = fopen(filename,"wb");
    check_file_was_successfully_open(initial_data, filename);
    fwrite(&npoints, sizeof(int), 1, initial_data);

    int failures = 0;
    for(int perturb=0;perturb<2;perturb++) {

      double rand_val[13];
      char suffix[10] = "";
      if (!perturb) {
        printf("Beginning standard data generation for routine %s\n", con2prim_test_names[which_routine]);
      } else {
        printf("Beginning perturbed data generation for routine %s\n", con2prim_test_names[which_routine]);
        srand(1000000);
        sprintf(suffix, "_pert");
        for(int i=0;i<13;i++) rand_val[i] = 1.0 + randf(-1,1)*1.0e-14;
      }

      FILE* outfiles[5];
      sprintf(filename,"apply_inequality_fixes%.5s.bin", suffix);
      outfiles[0] = fopen(filename,"wb");
      check_file_was_successfully_open(outfiles[0], filename);

      sprintf(filename,"%.30s_Hybrid_Multi_Method%.5s.bin", con2prim_test_names[which_routine], suffix);
      outfiles[1] = fopen(filename,"wb");
      check_file_was_successfully_open(outfiles[1], filename);

      sprintf(filename,"font_fix%.5s.bin", suffix);
      outfiles[2] = fopen(filename,"wb");
      check_file_was_successfully_open(outfiles[2], filename);

      sprintf(filename,"enforce_primitive_limits_and_compute_u0%.5s.bin", suffix);
      outfiles[3] = fopen(filename,"wb");
      check_file_was_successfully_open(outfiles[3], filename);

      sprintf(filename,"compute_conservs_and_Tmunu%.5s.bin", suffix);
      outfiles[4] = fopen(filename,"wb");
      check_file_was_successfully_open(outfiles[4], filename);

      srand(0);

      for(int j=0;j<npoints;j++) { // Density loop
        const double xrho  = exp(lrmin + dlr*j);
        double P_cold = 0.0;
        double eps_cold = 0.0;
        eos.hybrid_compute_P_cold_and_eps_cold(&eos, xrho, &P_cold, &eps_cold);

        // Compute the pressure step size
        const double lpmin        = log(1.0e-30);//-P_cold);
        const double lpmax        = log(10.0*P_cold);
        const double dlp          = (lpmax - lpmin)/(npoints-1);
        for(int i=0;i<npoints;i++) { // Pressure loop
          // Start by setting the prims (rho,Ye,T,P,eps)
          //double xtemp = exp(ltmin + dlt*i);
          //double xye   = Ye_test;
          const double xpress  = exp(lpmin + dlp*i);
          //WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xpress,&xeps );

          // Define the various GRHayL structs for the unit tests
          con2prim_diagnostics diagnostics;
          initialize_diagnostics(&diagnostics);
          metric_quantities metric;
          primitive_quantities prims;
          conservative_quantities cons, cons_undens;
          stress_energy Tmunu;

          // Generate random data to serve as the 'true' primitive values
          // and a randomized metric
          double lapse, betax, betay, betaz, gxx, gxy, gxz, gyy, gyz, gzz;
          randomize_metric(&lapse, &gxx, &gxy, &gxz, &gyy, &gyz, &gzz, &betax, &betay, &betaz);

          double eps, vx, vy, vz, Bx, By, Bz;
          randomize_primitives(&eos, xrho, xpress, &eps, &vx, &vy, &vz, &Bx, &By, &Bz);

      initialize_metric(lapse,
                        gxx, gxy, gxz,
                        gyy, gyz, gzz,
                        betax, betay, betaz,
                        &metric);

      initialize_primitives(xrho, xpress, eps,
                            vx, vy, vz,
                            Bx, By, Bz,
                            poison, poison, poison, // entropy, Y_e, temp
                            &prims);

      int speed_limit=0;
          limit_v_and_compute_u0(&eos, &metric, &prims, &speed_limit);

          // We need epsilon to compute the enthalpy in compute_conservs_and_Tmunu;
          // This normally happens in the enforce_primitive_limits_and_compute_u0 function
          prims.eps = eps_cold + (prims.press-P_cold)/(eos.Gamma_th-1.0)/prims.rho;

          // Compute conservatives based on these primitives
          compute_conservs_and_Tmunu(&params, &metric, &prims, &cons, &Tmunu);

          //This is meant to simulate some round-off error that deviates from the "true" values that we just computed.
          if(perturb) perturb_data(rand_val, &prims, &cons);

          if(!perturb)
            write_metric_binary(&metric, initial_data);

          int check = 0;
          if(cons.rho > 0.0) {

            //This applies the inequality (or "Faber") fixes on the conservatives
            if(eos.eos_type == 0) { //Hybrid-only
              if(!perturb && con2prim_test_keys[which_routine] == Noble2D) {
                write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims, outfiles[0]);
                write_conservative_binary(params.evolve_entropy, &cons, outfiles[0]);
              }
              apply_inequality_fixes(&params, &eos, &metric, &prims, &cons, &diagnostics);
              if(con2prim_test_keys[which_routine] == Noble2D)
                write_conservative_binary(params.evolve_entropy, &cons, outfiles[0]);
            }

            if(!perturb) {
              write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims, outfiles[1]);
              write_conservative_binary(params.evolve_entropy, &cons, outfiles[1]);
            }
            // The Con2Prim routines require the undensitized variables, but IGM evolves the densitized variables.
            undensitize_conservatives(&metric, &cons, &cons_undens);

            /************* Conservative-to-primitive recovery ************/
            check = grhayl_con2prim_multi_method(&params, &eos, &metric, &cons_undens, &prims, &diagnostics);
            // If the returned value is 5, then the Newton-Rapson method converged, but the values were so small
            // that u or rho were negative (usually u). Since the method converged, we only need to fix the values
            // using enforce_primitive_limits_and_compute_u0(). There's no need to trigger a Font fix.
            if(check==5) check = 0;

            write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims, outfiles[1]);

            if(con2prim_test_keys[which_routine] != Noble2D) continue;

            // Technically, font fix only needs the B values from prims; could reduce output by
            // only writing B's instead of all prims
            if(check!=0) {
              if(!perturb) {
                write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims, outfiles[2]);
                write_conservative_binary(params.evolve_entropy, &cons, outfiles[2]);
              }
              check = font_fix(&params, &eos, &metric, &cons, &prims, &diagnostics);
              write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims, outfiles[2]);
            } else { //The else is so that Font Fix is tested even when the primary routine succeeds.
              if(!perturb) {
                write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims, outfiles[2]);
                write_conservative_binary(params.evolve_entropy, &cons, outfiles[2]);
              }
              primitive_quantities prims_tmp;
              check = font_fix(&params, &eos, &metric, &cons, &prims_tmp, &diagnostics);
              if(con2prim_test_keys[which_routine] == Noble2D)
                write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims_tmp, outfiles[2]);
            }

            /*************************************************************/

            if(check)
              printf("Con2Prim and Font fix failed!");

          } else {
            diagnostics.failure_checker+=1;
            reset_prims_to_atmosphere(&eos, &prims);
            //TODO: Validate reset? (rhob press v)
            printf("Negative rho_* triggering atmospheric reset.\n");
          } // if rho_star > 0

          //--------------------------------------------------
          //---------- Primitive recovery completed ----------
          //--------------------------------------------------
          // Enforce limits on primitive variables and recompute conservatives.
          if(!perturb) {
            write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims, outfiles[3]);
          }
          enforce_primitive_limits_and_compute_u0(&params, &eos, &metric, &prims, &diagnostics.failure_checker);
          write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims, outfiles[3]);
          fwrite(&prims.u0, sizeof(double), 1, outfiles[3]);

          if(!perturb) {
            write_primitive_binary(eos.eos_type, params.evolve_entropy, &prims, outfiles[4]);
            fwrite(&prims.u0, sizeof(double), 1, outfiles[4]);
          }
          compute_conservs_and_Tmunu(&params, &metric, &prims, &cons, &Tmunu);
          write_conservative_binary(params.evolve_entropy, &cons, outfiles[4]);
          write_stress_energy_binary(&Tmunu, outfiles[4]);

          if( check != 0 )
            failures++;
        } // Pressure loop
      } // Density loop
      for(int k = 0; k < (sizeof(outfiles)/sizeof(outfiles[0])); k++) fclose(outfiles[k]);
    } // perturbation loop
    fclose(initial_data);

    int ntotal = npoints*npoints;

    printf("Completed data generation for routine %s\n",con2prim_test_names[which_routine]);
    printf("Final report:\n");
    printf("    Number of recovery attempts: %d\n",ntotal);
    printf("    Number of failed recoveries: %d\n",failures);
    printf("    Recovery failure rate      : %.2lf%%\n",((double)failures)/((double)ntotal)*100.0);
  }
  return 0;
}
