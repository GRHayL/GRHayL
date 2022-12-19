// Thorn      : GRHayL
// File       : con2prim_test_suite.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "stdlib.h"
#include "../con2prim_gem.h"
#include "../../EOS/Hybrid/EOS_hybrid.h"

inline double randf(double low,double high) {
    return (rand()/(double)(RAND_MAX))*(high-low)+low;
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

inline double relative_error( const double a, const double b ) {
  if     ( a != 0 ) return( fabs(1.0-b/a) );
  else if( b != 0 ) return( fabs(1.0-a/b) );
  else              return( 0.0 );
}

inline void output_primitive_error(
                     const primitive_quantities *restrict prims_orig, 
                     const primitive_quantities *restrict prims, 
                     const int i, const int j, FILE* outfile) {

  primitive_quantities prims_error;

  prims_error.rho   = relative_error(prims->rho,   prims_orig->rho);
  prims_error.press = relative_error(prims->press, prims_orig->press);
  prims_error.vx    = relative_error(prims->vx,    prims_orig->vx);
  prims_error.vy    = relative_error(prims->vy,    prims_orig->vy);
  prims_error.vz    = relative_error(prims->vz,    prims_orig->vz);
  
  fprintf(outfile, "%d %d %.15e %.15e %.15e"
                   " %.15e %.15e %.15e"
                   " %.15e %.15e %.15e"
                   " %.15e %.15e %.15e "
                   " %.15e %.15e %.15e\n",
                   i, j, prims_orig->rho, prims->rho, prims_error.rho,
                   prims_orig->press, prims->press, prims_error.press,
                   prims_orig->vx, prims->vx, prims_error.vx,
                   prims_orig->vy, prims->vy, prims_error.vy,
                   prims_orig->vz, prims->vz, prims_error.vz);
}

inline void output_conservative_error(
                     const conservative_quantities *restrict cons_orig, 
                     const conservative_quantities *restrict cons, 
                     const int i, const int j, FILE* outfile) {

  conservative_quantities cons_error;

  cons_error.rho = relative_error(cons->rho, cons_orig->rho);
  cons_error.tau = relative_error(cons->tau, cons_orig->tau);
  cons_error.S_x = relative_error(cons->S_x, cons_orig->S_x);
  cons_error.S_y = relative_error(cons->S_y, cons_orig->S_y);
  cons_error.S_z = relative_error(cons->S_z, cons_orig->S_z);

  fprintf(outfile, "%d %d %.15e %.15e %.15e"
                   " %.15e %.15e %.15e"
                   " %.15e %.15e %.15e "
                   " %.15e %.15e %.15e\n",
                   i, j, cons_orig->tau, cons->tau, cons_error.tau,
                   cons_orig->S_x, cons->S_x, cons_error.S_x,
                   cons_orig->S_y, cons->S_y, cons_error.S_y,
                   cons_orig->S_z, cons->S_z, cons_error.S_z);
}

void output_stress_energy_error(
                     const stress_energy *restrict Tmunu_orig, 
                     const stress_energy *restrict Tmunu,
                     const int i, const int j, FILE* outfile) {

  stress_energy Tmunu_error;

  Tmunu_error.Ttt = relative_error(Tmunu->Ttt, Tmunu_orig->Ttt);
  Tmunu_error.Ttx = relative_error(Tmunu->Ttx, Tmunu_orig->Ttx);
  Tmunu_error.Tty = relative_error(Tmunu->Tty, Tmunu_orig->Tty);
  Tmunu_error.Ttz = relative_error(Tmunu->Ttz, Tmunu_orig->Ttz);
  Tmunu_error.Txx = relative_error(Tmunu->Txx, Tmunu_orig->Txx);
  Tmunu_error.Txy = relative_error(Tmunu->Txy, Tmunu_orig->Txy);
  Tmunu_error.Txz = relative_error(Tmunu->Txz, Tmunu_orig->Txz);
  Tmunu_error.Tyy = relative_error(Tmunu->Tyy, Tmunu_orig->Tyy);
  Tmunu_error.Tyz = relative_error(Tmunu->Tyz, Tmunu_orig->Tyz);
  Tmunu_error.Tzz = relative_error(Tmunu->Tzz, Tmunu_orig->Tzz);

  fprintf(outfile, "%d %d %.15e %.15e %.15e %.15e %.15e %.15e "
                   "%.15e %.15e %.15e %.15e %.15e %.15e "
                   "%.15e %.15e %.15e %.15e %.15e %.15e "
                   "%.15e %.15e %.15e %.15e %.15e %.15e "
                   "%.15e %.15e %.15e %.15e %.15e %.15e\n",
                   i, j, Tmunu_orig->Ttt, Tmunu->Ttt, Tmunu_error.Ttt,
                   Tmunu_orig->Ttx, Tmunu->Ttx, Tmunu_error.Ttx,
                   Tmunu_orig->Tty, Tmunu->Tty, Tmunu_error.Tty,
                   Tmunu_orig->Ttz, Tmunu->Ttz, Tmunu_error.Ttz,
                   Tmunu_orig->Txx, Tmunu->Txx, Tmunu_error.Txx,
                   Tmunu_orig->Txy, Tmunu->Txy, Tmunu_error.Txy,
                   Tmunu_orig->Txz, Tmunu->Txz, Tmunu_error.Txz,
                   Tmunu_orig->Tyy, Tmunu->Tyy, Tmunu_error.Tyy,
                   Tmunu_orig->Tyz, Tmunu->Tyz, Tmunu_error.Tyz,
                   Tmunu_orig->Tzz, Tmunu->Tzz, Tmunu_error.Tzz);
}

void con2prim_unit_test( CCTK_ARGUMENTS ) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

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
  double tau_atm = 4.876083025795607e-12; //Taken from magnetizedTOV.par
  double W_max = 10.0; //IGM default
  double rho_b_max = 1e300; //IGM default
  double gamma_th = 2.0; //Taken from magnetizedTOV.par

//TODO: I need to fill in rho_tab, gamma_tab, k_tab, eps_tab based on the code, not like this
  double rho_tab[1] = {0.0};
  double gamma_tab[1] = {2.0};
  double k_tab = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  GRHayL_parameters params;
  initialize_GRHayL(None, backup_routine, false, false, calc_prims_guess, Psi6threshold, update_Tmunu, Cupp_Fix, &params);

  eos_parameters eos;
  initialize_general_eos(eos_type, tau_atm, W_max,
             poison, poison, poison, //entropy
             ut_rho_min, ut_rho_min, rho_b_max,
             &eos);

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

  // Number of points in the discretization of rho and T
  const int npoints       = ut_npoints;

  const double test_rho_min = ut_rho_min;
  const double test_rho_max = ut_rho_max;

//  const double test_T_min   = ut_T_min;
//  const double test_T_max   = ut_T_max;

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

      FILE* outfiles[8];
      char filename[100];

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_Summary.asc",con2prim_test_names[which_routine], suffix);
      outfiles[0] = fopen(filename,"w");

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_limit_v_and_output_u0.asc",con2prim_test_names[which_routine], suffix);
      outfiles[1] = fopen(filename,"w");
      fprintf(outfiles[1], "#Each variable has three columns: original value, new value, relative difference\n");
      fprintf(outfiles[1], "#i j vx vy vz\n");

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_apply_inequality_fixes.asc",con2prim_test_names[which_routine], suffix);
      outfiles[2] = fopen(filename,"w");
      fprintf(outfiles[2], "#Each variable has three columns: original value, new value, relative difference\n");
      fprintf(outfiles[2], "#i j tau S_x S_y S_z\n");

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_Hybrid_Multi_Method.asc",con2prim_test_names[which_routine], suffix);
      outfiles[3] = fopen(filename,"w");
      fprintf(outfiles[3], "#Each variable has three columns: original value, new value, relative difference\n");
      fprintf(outfiles[3], "#i j rho_b press vx vy vz\n");

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_font_fix.asc",con2prim_test_names[which_routine], suffix);
      outfiles[4] = fopen(filename,"w");
      fprintf(outfiles[4], "#Each variable has three columns: original value, new value, relative difference\n");
      fprintf(outfiles[4], "#i j rho_b press vx vy vz\n");

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_enforce_primitive_limits_and_output_u0.asc",con2prim_test_names[which_routine], suffix);
      outfiles[5] = fopen(filename,"w");
      fprintf(outfiles[5], "#Each variable has three columns: original value, new value, relative difference\n");
      fprintf(outfiles[5], "#i j rho_b press vx vy vz\n");

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_compute_conservs.asc",con2prim_test_names[which_routine], suffix);
      outfiles[6] = fopen(filename,"w");
      fprintf(outfiles[6], "#Each variable has three columns: original value, new value, relative difference\n");
      fprintf(outfiles[6], "#i j rho_* tau S_x S_y S_z\n");

      sprintf(filename,"unit_test/C2P_%.30s_%.4s_compute_Tmunu.asc",con2prim_test_names[which_routine], suffix);
      outfiles[7] = fopen(filename,"w");
      fprintf(outfiles[7], "#Each variable has three columns: original value, new value, relative difference\n");
      fprintf(outfiles[7], "#i j Ttt Ttx Tty Ttz Txx Txy Txz Tyy Tyz Tzz\n");

      srand(0);

      for(int i=0;i<npoints;i++) { // Density loop
        double xrho  = exp(lrmin + dlr*i);
        double P_cold = 0.0;
        compute_P_cold(&eos, xrho, &P_cold);

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

          // Generate random data to serve as the 'true' primitive values
          bool random_metric = true;
          initial_random_data(xrho, xpress, random_metric, &metric, &prims);

          double u0 = poison;
          prims_orig = prims;
          limit_v_and_output_u0(&eos, &metric, &prims, &u0, &diagnostics);
          fprintf(outfiles[1], "%d %d %.15e %.15e %.15e"
                               " %.15e %.15e %.15e"
                               " %.15e %.15e %.15e\n",
                               i, j, prims_orig.vx, prims.vx, relative_error(prims.vx, prims_orig.vx),
                               prims_orig.vy, prims.vy, relative_error(prims.vy, prims_orig.vy),
                               prims_orig.vz, prims.vz, relative_error(prims.vz, prims_orig.vz));

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
              output_conservative_error(&cons_orig, &cons, i, j, outfiles[2]);
            }

            // The Con2Prim routines require the undensitized variables, but IGM evolves the densitized variables.
            undensitize_conservatives(&metric, &cons, &cons_undens);

            /************* Conservative-to-primitive recovery ************/

            prims_orig = prims;
            check = Hybrid_Multi_Method(&params, &eos, &metric, &cons_undens, &prims, &prims_guess, &diagnostics);
            output_primitive_error(&prims_orig, &prims_guess, i, j, outfiles[3]);
  
            prims_orig = prims;
            if(check!=0) {
              check = font_fix(&eos, &metric, &cons, &prims, &prims_guess, &diagnostics);
              diagnostics.font_fixes++;
              output_primitive_error(&prims_orig, &prims_guess, i, j, outfiles[4]);
            } else { //The else is so that Font Fix is tested even when the primary routine succeeds.
              primitive_quantities prims_tmp;
              check = font_fix(&eos, &metric, &cons, &prims, &prims_tmp, &diagnostics);
              output_primitive_error(&prims_orig, &prims_tmp, i, j, outfiles[4]);
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
          output_primitive_error(&prims_orig, &prims, i, j, outfiles[5]);

          cons_orig = cons;
          Tmunu_orig = Tmunu;
          compute_conservs_and_Tmunu(&params, &eos, &metric, &prims, u0, &cons, &Tmunu);
          output_conservative_error(&cons_orig, &cons, i, j, outfiles[6]);
          output_stress_energy_error(&Tmunu_orig, &Tmunu, i, j, outfiles[7]);

          if( check != 0 ) {
            failures++;
            fprintf(outfiles[0],"Recovery FAILED!\n");
            printf("Recovery FAILED!\n");
          } else {
            fprintf(outfiles[0], "Recovery SUCCEEDED FOR POINT (%d,%d)!\n", i,j);
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

  }
  printf("All done! Terminating the run.\n");
  exit(1);
}
