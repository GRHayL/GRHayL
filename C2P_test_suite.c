// Thorn      : IllinoisGRMHD
// File       : con2prim_test_suite.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: In this file we provide an extensive test suite of
//              the con2prim routines available in the code.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "stdlib.h"
#include "con2prim_header.h"

inline double randf(double low,double high) {
    return (rand()/(double)(RAND_MAX))*(high-low)+low;
}

inline double relative_error( const double a, const double b ) {
  if     ( a != 0 ) return( fabs(1.0-b/a) );
  else if( b != 0 ) return( fabs(1.0-a/b) );
  else              return( 0.0 );
}

void C2P_test_suite( CCTK_ARGUMENTS ) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  double poison = 1e200;
  int main = Noble2D;
  int backup_routine[3] = {None,None,None};
  int update_Tmunu = 1; //IGM default
  int eos_type = 0; // Hybrid=0, Tabulated=1;
  int neos = 1;

  double Psi6threshold = 1e100; //Taken from magnetizedTOV.par
  double tau_atm = 4.876083025795607e-12; //Taken from magnetizedTOV.par
  double W_max = 10.0; //IGM default
  double rho_b_atm = 1.292852735094440e-10; //Taken from magnetizedTOV.par
  double rho_b_min = rho_b_atm;
  double rho_b_max = 1e300; //IGM default
  double gamma_th = 2.0; //Taken from magnetizedTOV.par

//TODO: I need to fill in rho_tab, gamma_tab, k_tab, eps_tab based on the code, not like this
  double rho_tab[1] = {0.0};
  double gamma_tab[1] = {2.0};
  double k_tab[1] = {1.0};
  double eps_tab[1] = {0.0};

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  GRMHD_parameters params;
  initialize_parameters(&params, main, backup_routine, false, false, true, Psi6threshold, update_Tmunu);

  eos_parameters eos;
  initialize_general_eos(&eos, eos_type, tau_atm, W_max,
             poison, poison,poison, //epsilon
             poison, poison, poison, //pressure
             poison, poison, poison, //entropy
             rho_b_min, rho_b_atm, rho_b_max);

  initialize_hybrid_eos(&eos, neos, rho_tab,
             gamma_tab, k_tab, eps_tab, gamma_th);

  con2prim_diagnostics diagnostics;
  initialize_diagnostics(&diagnostics);

  // Count number of routines tested
  int num_routines_tested = 1;
  int con2prim_test_keys[4];
  char con2prim_test_names[4][50];
  con2prim_test_keys[0] = params.main_routine;
  sprintf(con2prim_test_names[0],"%s","Noble2D");
  if( params.backup_routine[0] != None ) {
    num_routines_tested++;
    con2prim_test_keys[1] = params.backup_routine[0];
//    sprintf(con2prim_test_names[1],"%s",igm_con2prim_backup_routine[0]);
    if( params.backup_routine[1] != None ) {
      num_routines_tested++;
      con2prim_test_keys[2] = params.backup_routine[1];
//      sprintf(con2prim_test_names[2],"%s",igm_con2prim_backup_routine[1]);
      if( params.backup_routine[2] != None ) {
        num_routines_tested++;
        con2prim_test_keys[3] = params.backup_routine[2];
//        sprintf(con2prim_test_names[3],"%s",igm_con2prim_backup_routine[2]);
      }
    }
  }

  // We will be performing the tabulated EOS test in the following way:
  //
  //      Y_e      = 0.1
  //       W       = 2
  // log10(Pmag/P) = -5
  //
  // rho will vary between rho_min and rho_max (uniformly in log space)
  //  T  will vary between  T_min  and  T_max  (uniformly in log space)

  // Number of points in the discretization of rho and T
  const int npoints       = C2P_ut_npoints;

  const double test_rho_min = C2P_ut_rho_min;
  const double test_rho_max = C2P_ut_rho_max;

  const double test_T_min   = C2P_ut_T_min;
  const double test_T_max   = C2P_ut_T_max;

  // Compute the density step size
  const double lrmin        = log(test_rho_min);
  const double lrmax        = log(test_rho_max);
  const double dlr          = (lrmax - lrmin)/(npoints-1);

  // Compute the temperature step size
  const double ltmin        = log(test_T_min);
  const double ltmax        = log(test_T_max);
  const double dlt          = (ltmax - ltmin)/(npoints-1);

  // Fix Y_e
  const double Ye_test = 0.1;
//TODO: make v vary via W
  // Fix W
  const double W_test  = 2.0;
  // Fix log10(Pmag/P)
  const double logPmoP = -5.0;

  // tau is given by (see :
  //
  // tau := hW^{2} + B^{2} - P - 0.5*( (B.v)^{2} + (B/W)^{2} )
  // Absolutely minimum allowed tau

  // Now perform one test for each of the selected routines
  for(int which_routine=0;which_routine<num_routines_tested;which_routine++) {

    int failures = 0;
    for(int rand=0;rand<2;rand++) {

      double rand_val[13];
      char suffix[10] = "norm";
      if(rand==1) {
        srand(1000000);
        sprintf(suffix, "rand");
        for(int i=0;i<13;i++) rand_val[i] = 1.0 + randf(-1,1)*1.0e-15;
      }

    printf("Beginning %s test\n", suffix);
//    printf("Beginning %s test for routine %s\n",con2prim_test_names[which_routine]);

      char filename[100];
      sprintf(filename,"C2P_Testsuite_%s_%s.asc",con2prim_test_names[which_routine], suffix);
      FILE* outfile = fopen(filename,"w");

      srand(0);

      for(int i=0;i<npoints;i++) { // Density loop
//      for(int j=0;j<npoints;j++) { // Temperature loop
        // Start by setting the prims (rho,Ye,T,P,eps)
        double xrho  = exp(lrmin + dlr*i);
//        double xtemp = exp(ltmin + dlt*j);
//        double xye   = Ye_test;
        double xpress  = 0.0;
        double xeps  = 0.0;
//        WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xpress,&xeps );
        compute_P_cold_and_eps_cold(&eos, xrho, &xpress, &xeps);

        // Now set the velocities
        // Velocity magnitude
        const double v = sqrt(1.0-1.0/(W_test*W_test));
        const double vx = v*randf(0.0,1.0);
        const double vy = sqrt(v*v - vx*vx)*randf(0.0,1.0);
        const double vz = sqrt(v*v - vx*vx - vy*vy);

        // Now the magnetic fields. We'll set them aligned
        // with the velocities, for simplicity.
        const double Bhatx = vx/v;
        const double Bhaty = vy/v;
        const double Bhatz = vz/v;
        const double B     = sqrt(2.0*pow(10.0,logPmoP)*xpress);
        const double Bx    = -Bhatx * B;
        const double By    = -Bhaty * B;
        const double Bz    = -Bhatz * B;

        const double gxx = 1.0 + randf(0.0,100.0);
        const double gyy = 1.0 + randf(0.0,100.0);
        const double gzz = 1.0 + randf(0.0,100.0);
        const double gxy = randf(-1.0e-1,1.0e-1);
        const double gxz = randf(-1.0e-1,1.0e-1);
        const double gyz = randf(-1.0e-1,1.0e-1);

        const double lapse = randf(-1.0e-10,1.0);

        const double betax = v*randf(0.0,1.0);
        const double betay = sqrt(v*v - betax*betax)*randf(0.0,1.0);
        const double betaz = sqrt(v*v - betax*betax - betay*betay);
printf("randomized betas %.16e %.16e %.16e", betax, betay, betaz);

//        const double gxx = 1.0;// + randf(0.0,100.0);
//        const double gyy = 1.0;// + randf(0.0,100.0);
//        const double gzz = 1.0;// + randf(0.0,100.0);
//        const double gxy = 0.0;//randf(-1.0e-1,1.0e-1);
//        const double gxz = 0.0;//randf(-1.0e-1,1.0e-1);
//        const double gyz = 0.0;//randf(-1.0e-1,1.0e-1);
//
//        const double lapse = 1.0;//randf(-1.0e-10,1.0);
//
//        const double betax = 0.0;//randf(0.0,1.0);
//        const double betay = 0.0;//sqrt(v*v - betax*betax)*randf(0.0,1.0);
//        const double betaz = 0.0;//sqrt(v*v - betax*betax - betay*betay);

        // Store the metric randomized values into the structs
        metric_quantities metric;
        initialize_metric(&metric, lapse,  // phi, psi, lapse
                          gxx, gxy, gxz,  // gxx, gxy, gxz
                          gyy, gyz, gzz,  // gyy, gyz, gzz
                          betax, betay, betaz); // betax, betay, betaz

        conservative_quantities cons, cons_orig, cons_undens; // Not initialized because it will be filled based on primitive data.

        primitive_quantities prims, prims_orig, prims_guess;
        initialize_primitives(&eos, &metric,
                              xrho, xpress, xeps,
                              vx, vy, vz, Bx, By, Bz,
                              poison, poison, poison, // entropy, Y_e=xye, temp=xtemp
                              &prims);
        prims.print = true;

        // Define the stress_energy struct
        stress_energy Tmunu;

        // Then set the conservative variables array
        double TUPMUNU[10],TDNMUNU[10];
        enforce_limits_on_primitives_and_recompute_conservs(&params, &eos, &metric, &prims, &cons,
                                                            TUPMUNU, TDNMUNU, &Tmunu, &diagnostics);

        // Store original prims
        prims_orig = prims;
        cons_orig = cons;

        printf("Initial primitives:\n %.16e\n %.16e\n"
           "  %.16e\n %.16e\n %.16e\n", prims.rho, prims.press, prims.vx, prims.vy, prims.vz);
        printf("Initial conservatives:\n %.16e\n %.16e\n"
           "  %.16e\n %.16e\n %.16e\n", cons.rho, cons.tau, cons.S_x, cons.S_y, cons.S_z);
        //This is meant to simulate some round-off error that deviates from the "true" values that we just computed.
        if(rand==1) {
          prims.rho   *= rand_val[0];
          prims.press *= rand_val[1];
          prims.vx    *= rand_val[2];
          prims.vy    *= rand_val[3];
          prims.vz    *= rand_val[4];
          prims.Bx    *= rand_val[5];
          prims.By    *= rand_val[6];
          prims.Bz    *= rand_val[7];
          cons.rho    *= rand_val[8];
          cons.S_x    *= rand_val[9];
          cons.S_y    *= rand_val[10];
          cons.S_z    *= rand_val[11];
          cons.tau    *= rand_val[12];
        printf("Perturbed conservatives:\n %.16e\n %.16e\n"
           "  %.16e\n %.16e\n %.16e\n", cons.rho, cons.tau, cons.S_x, cons.S_y, cons.S_z);
        }

        int check = 0;
        if(cons.rho > 0.0) {

          //This applies the inequality (or "Faber") fixes on the conservatives
          if(eos.eos_type == 0) //Hybrid-only
            apply_tau_floor(&params, &eos, &metric, &prims, &cons, &diagnostics);

          // The Con2Prim routines require the undensitized variables, but IGM evolves the densitized variables.
          undensitize_conservatives(&eos, con2prim_test_keys[which_routine], &metric, &prims, &cons, &cons_undens);

          // The con2prim routines require primitive guesses in order to perform
          // the recovery. In IllinoisGRMHD, we do not keep track of the primitives
          // in between time steps, and therefore our guesses are *not* the
          // values of the primitives in the previous time level. Instead, we provide
          // guesses based on the conservative variables.

          /************* Conservative-to-primitive recovery ************/

          for(int which_guess=1;which_guess<=2;which_guess++) {
            guess_primitives(&eos, con2prim_test_keys[which_routine], which_guess,
                             &metric, &prims, &cons, &prims_guess);

/*Why was it getting set to 1e300 right after setting the guess?
          for(int i=0;i<numprims;i++) prim[i] = 1e300;
          // prim[TEMP    ] = eos.T_max;
          // prim[RHO     ] = PRIMS[RHOB       ];
          // prim[YE      ] = PRIMS[YEPRIM     ];
          prim[TEMP    ] = PRIMS[TEMPERATURE]*0.95;
          // prim[UTCON1  ] = PRIMS[VX         ];
          // prim[UTCON2  ] = PRIMS[VY         ];
          // prim[UTCON3  ] = PRIMS[VZ         ];
          prim[WLORENTZ] = W_test*0.95;
          // prim[B1_con  ] = PRIMS[BX_CENTER  ];
          // prim[B2_con  ] = PRIMS[BY_CENTER  ];
          // prim[B3_con  ] = PRIMS[BZ_CENTER  ];
*/
            check = C2P_Select_Hybrid_Method(&eos, con2prim_test_keys[which_routine], &metric, &cons_undens, &prims_guess, &diagnostics);
            //If multiple C2P routines are selected as backups, the backup routine logic would go here
  
            if(check!=0) {
              printf("Applying Font Fix\n");
              check = font_fix(&eos, &metric, &cons, &prims, &prims_guess, &diagnostics);
              diagnostics.font_fixes++;
            }
      /*************************************************************/
  
            if(check==0) {
              prims = prims_guess;
              which_guess=3; //TODO: can we make the multiple guess loop cleaner?
            } else {
              printf("Con2Prim and Font fix failed!");
              printf("diagnostics->failure_checker = %d st_i = %e %e %e, rhostar = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                      diagnostics.failure_checker, cons_orig.S_x, cons_orig.S_y, cons_orig.S_z, cons_orig.rho, prims.Bx, prims.By, prims.Bz,
                      metric.adm_gxx, metric.adm_gxy, metric.adm_gxz, metric.adm_gyy, metric.adm_gyz, metric.adm_gzz, metric.psi6);
            }
          } //If we didn't find a root, then try again with a different guess.
        } else {
          diagnostics.failure_checker+=1;
          reset_prims_to_atmosphere(&eos, &prims, &diagnostics);
          diagnostics.rho_star_fix_applied++;
        } // if rho_star > 0

        //--------------------------------------------------
        //---------- Primitive recovery completed ----------
        //--------------------------------------------------
        // Enforce limits on primitive variables and recompute conservatives.
        printf("C2P primitives:\n %.16e\n %.16e\n"
           "  %.16e\n %.16e\n %.16e\n", prims.rho, prims.press, prims.vx, prims.vy, prims.vz);
        enforce_limits_on_primitives_and_recompute_conservs(&params, &eos, &metric, &prims, &cons,
                                                            TUPMUNU, TDNMUNU, &Tmunu, &diagnostics);
        printf("Enforced primitives:\n %.16e\n %.16e\n"
           "  %.16e\n %.16e\n %.16e\n", prims.rho, prims.press, prims.vx, prims.vy, prims.vz);

        primitive_quantities prims_error;
        double accumulated_error = 0.0;
        if( check != 0 ) {
          failures++;
          fprintf(outfile,"Recovery FAILED!\n");
          printf("Recovery FAILED!\n");
          accumulated_error = 1e300;
        } else {
          fprintf(outfile, "Recovery SUCCEEDED!\n");
          printf("Recovery SUCCEEDED!\n");
          prims_error.rho     = relative_error(prims.rho,     prims_orig.rho);
          prims_error.press   = relative_error(prims.press,   prims_orig.press);
          prims_error.eps     = relative_error(prims.eps,     prims_orig.eps);
          prims_error.vx      = relative_error(prims.vx,      prims_orig.vx);
          prims_error.vy      = relative_error(prims.vy,      prims_orig.vy);
          prims_error.vz      = relative_error(prims.vz,      prims_orig.vz);
          prims_error.Bx      = relative_error(prims.Bx,      prims_orig.Bx);
          prims_error.By      = relative_error(prims.By,      prims_orig.By);
          prims_error.Bz      = relative_error(prims.Bz,      prims_orig.Bz);
          prims_error.entropy = relative_error(prims.entropy, prims_orig.entropy);
          prims_error.Y_e     = relative_error(prims.Y_e,     prims_orig.Y_e);
          prims_error.temp    = relative_error(prims.temp,    prims_orig.temp);

//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "rho_b", prims_orig.rho, prims.rho, prims_error.rho);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "press", prims_orig.press, prims.press, prims_error.press);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "eps", prims_orig.eps, prims.eps, prims_error.eps);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "vx", prims_orig.vx, prims.vx, prims_error.vx);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "vy", prims_orig.vy, prims.vy, prims_error.vy);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "vz", prims_orig.vz, prims.vz, prims_error.vz);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "Bx", prims_orig.Bx, prims.Bx, prims_error.Bx);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "By", prims_orig.By, prims.By, prims_error.By);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "Bz", prims_orig.Bz, prims.Bz, prims_error.Bz);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "entropy", prims_orig.entropy, prims.entropy, prims_error.entropy);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "Y_e", prims_orig.Y_e, prims.Y_e, prims_error.Y_e);
//          fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "temp", prims_orig.temp, prims.temp, prims_error.temp);
          fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "rho_b", prims_orig.rho, prims.rho, prims_error.rho);
          fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "press", prims_orig.press, prims.press, prims_error.press);
          fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "vx", prims_orig.vx, prims.vx, prims_error.vx);
          fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "vy", prims_orig.vy, prims.vy, prims_error.vy);
          fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "vz", prims_orig.vz, prims.vz, prims_error.vz);
//          fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "Bx", prims_orig.Bx, prims.Bx, prims_error.Bx);
//          fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "By", prims_orig.By, prims.By, prims_error.By);
//          fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "Bz", prims_orig.Bz, prims.Bz, prims_error.Bz);

          accumulated_error = prims_error.rho + prims_error.press /*+ prims_error.eps*/ + prims_error.vx + prims_error.vy + prims_error.vz
                            + prims_error.Bx + prims_error.By + prims_error.Bz /*+ prims_error.entropy + prims_error.Y_e + prims_error.temp*/;

          fprintf(outfile, "Total accumulated error    : %e\n",accumulated_error);
        }
        fprintf(outfile,"%e %e\n\n",log10(prims_orig.rho), log10(MAX(accumulated_error,1e-16)));
//      } // temperature loop
      }
    fclose(outfile);
    }


    int ntotal = npoints*npoints;

    printf("Completed test for routine %s\n",con2prim_test_names[which_routine]);
    printf("Final report:\n");
    printf("    Number of recovery attempts: %d\n",ntotal);
    printf("    Number of failed recoveries: %d\n",failures);
    printf("    Recovery failure rate      : %.2lf%%\n",((double)failures)/((double)ntotal)*100.0);

  }

//  int exp_num = 2;
//  int failures = 0;
//
//  double explicit_rho_b[] = {1.2928527350944400e-10, 1.2928527350944400e-10};
//  double explicit_P[] = {1.5043213751770567e-20, 1.5043213751770567e-20};
//  double explicit_vx[] = {0.0, 0.0, 0.0};
//  double explicit_vy[] = {0.0, 0.0, 0.0};
//  double explicit_vz[] = {0.0, 0.0, 0.0};
//  double explicit_Bx[] = {0.0, 0.0, 0.0};
//  double explicit_By[] = {0.0, 0.0, 0.0};
//  double explicit_Bz[] = {0.0, 0.0, 0.0};
//  double explicit_rho_star[] = {1.3234724020689511e-10, 1.3234724020689511e-10};
//  double explicit_tau[] = {1.5399466368020248e-20, 1.5399466368068784e-20};
//  double explicit_Sx[] = {8.7368463342954066e-15, 8.7368463342954066e-15};
//  double explicit_Sy[] = {7.1465425792354919e-15, 7.1465425792387736e-15};
//  double explicit_Sz[] = {8.7368463342954066e-15, 8.7368463342954066e-15};
//
//  double out_rhob=1.2928527350944400e-10;
//  double out_P=1.5043213751770567e-20;
//  double out_vx = 9.1927394657560273e-06;
//  double out_vy = 7.5194528435267824e-06;
//  double out_vz = 9.1927394657560273e-06;
//  double out_Bx = 0.0;
//  double out_By = 0.0;
//  double out_Bz = 0.0;
//  double out_rhos = 1.3234724034565631e-10;
//  double out_tau = 3.0798958584463867e-20;
//  double out_Sx = 1.2454672417374646e-15;
//  double out_Sy = 1.0187640177651899e-15;
//  double out_Sz = 1.2454672417374646e-15;
//  
//  for(int which_routine=0;which_routine<num_routines_tested;which_routine++) {
//
//    CCTK_VINFO("Beginning explicit value test for routine %s",con2prim_test_names[which_routine]);
//
//    char filename[100];
//    sprintf(filename,"C2P_Testsuite_%s.asc",con2prim_test_names[which_routine]);
//    FILE* outfile = fopen(filename,"a");
//    fprintf(outfile, "Beginning explicit data tests...");
//
//    int failures = 0;
//
//    srand(0);
//
//    for(int i=0;i<exp_num;i++) {
//
//      // Set the metric to flat space
//      metric_quantities metric;
//      initialize_metric(&metric, 0.0, 0.0, 1.0,  // phi, psi, lapse
//                                 1.0, 0.0, 0.0,  // gxx, gxy, gxz
//                                 1.0, 0.0, 1.0,  // gyy, gyz, gzz
//                                 0.0, 0.0, 0.0); // betax, betay, betaz
//
//      conservative_quantities cons, cons_undens;
//      initialize_conservatives(&params, &eos,
//                               explicit_rho_star[i], explicit_tau[i],
//                               explicit_Sx[i], explicit_Sy[i], explicit_Sz[i],
//                               poison, poison, // entropy, Y_e
//                               &cons);
//
//      primitive_quantities prims, prims_orig, prims_guess;
//      initialize_primitives(&eos, &metric,
//                            explicit_rho_b[i], explicit_P[i], poison,
//                            explicit_vx[i], explicit_vy[i], explicit_vz[i],
//                            explicit_Bx[i], explicit_By[i], explicit_Bz[i],
//                            poison, poison, poison, // entropy, Y_e=xye, temp=xtemp
//                            &prims);
//
//      // Store original prims
//      prims_orig = prims;
//
//      int check = 0;
//      if(cons.rho>0.0) {
//
//        //This applies the inequality (or "Faber") fixes on the conservatives
//        apply_tau_floor(&params, &eos, &metric, &prims, &cons, &diagnostics);
//
//        // The Con2Prim routines require the undensitized variables, but IGM evolves the densitized variables.
//        undensitize_conservatives(&eos, con2prim_test_keys[which_routine], &metric, &prims, &cons, &cons_undens);
//
//        // The con2prim routines require primitive guesses in order to perform
//        // the recovery. In IllinoisGRMHD, we do not keep track of the primitives
//        // in between time steps, and therefore our guesses are *not* the
//        // values of the primitives in the previous time level. Instead, we provide
//        // guesses based on the conservative variables.
//
//        for(int which_guess=1;which_guess<=2;which_guess++) {
//          guess_primitives(&eos, con2prim_test_keys[which_routine], which_guess,
//                           &metric, &prims, &cons, &prims_guess);
//
//          check = C2P_Select_Hybrid_Method(&eos, con2prim_test_keys[which_routine], &metric, &cons_undens, &prims_guess, &diagnostics);
//          if(check!=0) {
//            CCTK_VINFO("Applying Font fix to explicit point %d", i);
//            check = font_fix(&eos, &metric, &cons_undens, &prims, &prims_guess, &diagnostics);
//            diagnostics.font_fixes++;
//          }
//          if( check == 0 ) {
//            prims = prims_guess;
//            which_guess = 3;
//          }
//        }
//      } else {
//        diagnostics.failure_checker+=1;
//        reset_prims_to_atmosphere(&eos, &prims, &diagnostics);
//        diagnostics.rho_star_fix_applied++;
//        CCTK_VINFO("Applying rho_* fix to explicit point %d with rho %e", i, cons.rho);
//      }
//      if(check != 0) {
//        CCTK_VINFO("C2P and FF failed for explicit point %d! Reseting to atm...", i);
//        // Sigh, reset to atmosphere
//        reset_prims_to_atmosphere( &eos, &prims, &diagnostics);
//        diagnostics.atm_resets++;
//        // Then flag this point as a "success"
//        check = 0;
//      }
//
//      //--------------------------------------------------
//      //---------- Primitive recovery completed ----------
//      //--------------------------------------------------
//      // Enforce limits on primitive variables and recompute conservatives.
//      double TUPMUNU[10],TDNMUNU[10];
//      stress_energy Tmunu;
//      enforce_limits_on_primitives_and_recompute_conservs(&params, &eos, &metric, &prims, &cons,
//                                                          TUPMUNU, TDNMUNU, &Tmunu, &diagnostics);
////CCTK_VINFO("Post-enforce %.16e %.16e\n   %.16e %.16e\n   %.16e", prims.rho, prims.press, prims.vx, prims.vy, prims.vz);
//
//      primitive_quantities prims_error;
//      conservative_quantities cons_error;
//      double total_prims_error = 0.0;
//      double total_cons_error = 0.0;
//      if( check != 0 ) {
//        failures++;
//        CCTK_VINFO("Recovery FAILED!\n");
//        total_prims_error = total_cons_error = 1e300;
//      } else {
//        prims = prims_guess;
//        CCTK_VINFO("Recovery SUCCEEDED!");
//        prims_error.rho   = relative_error(prims.rho,   out_rhob);
//        prims_error.press = relative_error(prims.press, out_P);
//        prims_error.vx    = relative_error(prims.vx,    out_vx);
//        prims_error.vy    = relative_error(prims.vy,    out_vy);
//        prims_error.vz    = relative_error(prims.vz,    out_vz);
//
//        cons_error.rho   = relative_error(cons.rho, out_rhos);
//        cons_error.tau = relative_error(cons.tau, out_tau);
//        cons_error.S_x    = relative_error(cons.S_x,  out_Sx);
//        cons_error.S_y    = relative_error(cons.S_y,  out_Sy);
//        cons_error.S_z    = relative_error(cons.S_z,  out_Sz);
//
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "rho_b", out_rhob, prims.rho, prims_error.rho);
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "press", out_P, prims.press, prims_error.press);
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "vx", out_vx, prims.vx, prims_error.vx);
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "vy", out_vy, prims.vy, prims_error.vy);
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "vz", out_vz, prims.vz, prims_error.vz);
//
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "rho_*", out_rhos, cons.rho, cons_error.rho);
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "tau", out_tau, cons.tau, cons_error.tau);
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "Sx", out_Sx, cons.S_x, cons_error.S_x);
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "Sy", out_Sy, cons.S_y, cons_error.S_y);
//        fprintf(outfile, "Relative error for %s: (%e -> %e) %.3e\n", "Sz", out_Sz, cons.S_z, cons_error.S_z);
//
//        total_prims_error = prims_error.rho + prims_error.press
//                          + prims_error.vx + prims_error.vy + prims_error.vz;
//
//        total_cons_error = cons_error.rho + cons_error.tau
//                          + cons_error.S_x + cons_error.S_y + cons_error.S_z;
//
//        fprintf("Total accumulated errors    : %e %e\n",total_prims_error, total_cons_error);
//      }
//      fprintf(outfile,"\n");
//    }
//
//    fclose(outfile);
//
//    int ntotal = npoints*npoints;
//
//    printf("Completed explicit data test for routine %s",con2prim_test_names[which_routine]);
//    printf("Final report:");
//    printf("    Number of recovery attempts: %d",exp_num);
//    printf("    Number of failed recoveries: %d",failures);
//    printf("    Recovery failure rate      : %.2lf%%",((double)failures)/((double)ntotal)*100.0);
//
//  }

  printf("All done! Terminating the run.\n");
  exit(1);
}
/*
  Conservative ID that should trigger Font fix.
Incoming conservatives:
  rho_*=2.3995320450572215e-06
  tau=3.8712996396831345e-09
  S=(-1.2764009790120381e-07,1.2763915658920790e-07,1.2714959053780215e-07)
*/
