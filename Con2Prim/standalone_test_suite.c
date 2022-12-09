// Thorn      : GrHayL
// File       : con2prim_test_suite.cc
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive test suite of
//              the Con2Prim gem.

#include "stdlib.h"
#include "con2prim_gem.h"
#include "../EOS/Hybrid/EOS_hybrid.h"

inline double randf(double low,double high) {
    return (rand()/(double)(RAND_MAX))*(high-low)+low;
}

inline double relative_error( const double a, const double b ) {
  if     ( a != 0 ) return( fabs(1.0-b/a) );
  else if( b != 0 ) return( fabs(1.0-a/b) );
  else              return( 0.0 );
}

inline void primitive_error(const primitive_quantities *restrict prims_orig, 
                        const primitive_quantities *restrict prims, 
                        primitive_quantities *restrict prims_error) {
          prims_error->rho     = relative_error(prims->rho,     prims_orig->rho);
          prims_error->press   = relative_error(prims->press,   prims_orig->press);
          prims_error->eps     = relative_error(prims->eps,     prims_orig->eps);
          prims_error->vx      = relative_error(prims->vx,      prims_orig->vx);
          prims_error->vy      = relative_error(prims->vy,      prims_orig->vy);
          prims_error->vz      = relative_error(prims->vz,      prims_orig->vz);
          prims_error->Bx      = relative_error(prims->Bx,      prims_orig->Bx);
          prims_error->By      = relative_error(prims->By,      prims_orig->By);
          prims_error->Bz      = relative_error(prims->Bz,      prims_orig->Bz);
          prims_error->entropy = relative_error(prims->entropy, prims_orig->entropy);
          prims_error->Y_e     = relative_error(prims->Y_e,     prims_orig->Y_e);
          prims_error->temp    = relative_error(prims->temp,    prims_orig->temp);
  
}

inline void conservative_error(const conservative_quantities *restrict cons_orig, 
                        const conservative_quantities *restrict cons, 
                        conservative_quantities *restrict cons_error) {
          cons_error->rho     = relative_error(cons->rho,     cons_orig->rho);
          cons_error->tau     = relative_error(cons->tau,     cons_orig->tau);
          cons_error->S_x     = relative_error(cons->S_x,     cons_orig->S_x);
          cons_error->S_y     = relative_error(cons->S_y,     cons_orig->S_y);
          cons_error->S_z     = relative_error(cons->S_z,     cons_orig->S_z);
          cons_error->entropy = relative_error(cons->entropy, cons_orig->entropy);
          cons_error->Y_e     = relative_error(cons->Y_e,     cons_orig->Y_e);
}

inline void stress_energy_error(const stress_energy *restrict Tmunu_orig, 
                        const stress_energy *restrict Tmunu, 
                        stress_energy *restrict Tmunu_error) {
          Tmunu_error->Ttt = relative_error(Tmunu->Ttt, Tmunu_orig->Ttt);
          Tmunu_error->Ttx = relative_error(Tmunu->Ttx, Tmunu_orig->Ttx);
          Tmunu_error->Tty = relative_error(Tmunu->Tty, Tmunu_orig->Tty);
          Tmunu_error->Ttz = relative_error(Tmunu->Ttz, Tmunu_orig->Ttz);
          Tmunu_error->Txx = relative_error(Tmunu->Txx, Tmunu_orig->Txx);
          Tmunu_error->Txy = relative_error(Tmunu->Txy, Tmunu_orig->Txy);
          Tmunu_error->Txz = relative_error(Tmunu->Txz, Tmunu_orig->Txz);
          Tmunu_error->Tyy = relative_error(Tmunu->Tyy, Tmunu_orig->Tyy);
          Tmunu_error->Tyz = relative_error(Tmunu->Tyz, Tmunu_orig->Tyz);
          Tmunu_error->Tzz = relative_error(Tmunu->Tzz, Tmunu_orig->Tzz);
}

void main() {

  double poison = 1e200;
  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  int main = Noble2D;
  int backup_routine[3] = {None,None,None};
  int update_Tmunu = 1; //IGM default
  int eos_type = 0; // Hybrid=0, Tabulated=1;
  int neos = 1;

  double Psi6threshold = 1e100; //Taken from magnetizedTOV.par
  double tau_atm = 4.876083025795607e-12; //Taken from magnetizedTOV.par
  double W_max = 10.0; //IGM default
  double rho_b_max = 1e300; //IGM default
  double gamma_th = 2.0; //Taken from magnetizedTOV.par

//TODO: I need to fill in rho_tab, gamma_tab, k_tab, eps_tab based on the code, not like this
  double rho_tab[1] = {0.0};
  double gamma_tab[1] = {2.0};
  double k_tab = 1.0;

  //Setting things that are set by the parfile for ET
  bool Cupp_Fix = "yes";
  int C2P_ut_npoints = 256;
  double C2P_ut_rho_min = 1e-12;
  double C2P_ut_rho_max = 1e-3;
  double C2P_ut_T_min = 1e-2;
  double C2P_ut_T_max = 1e+2;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  GRHayL_parameters params;
  initialize_GRHayL(main, backup_routine, false, false, true, Psi6threshold, update_Tmunu, Cupp_Fix, &params);

  eos_parameters eos;
  initialize_general_eos(eos_type, tau_atm, W_max,
             poison, poison, poison, //entropy
             C2P_ut_rho_min, C2P_ut_rho_min, rho_b_max,
             &eos);

  initialize_hybrid_eos(neos, rho_tab,
             gamma_tab, k_tab, gamma_th,
             &eos);

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
  //const double ltmin        = log(test_T_min);
  //const double ltmax        = log(test_T_max);
  //const double dlt          = (ltmax - ltmin)/(npoints-1);

  // Fix Y_e
  //const double Ye_test = 0.1;
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
        for(int i=0;i<13;i++) rand_val[i] = 1.0 + randf(-1,1)*1.0e-14;
      }

    printf("Beginning %s test\n", suffix);
//    printf("Beginning %s test for routine %s\n",con2prim_test_names[which_routine]);

      char filename[100];
      sprintf(filename,"unit_test/C2P_Testsuite_%s_%s.asc",con2prim_test_names[which_routine], suffix);
      FILE* outfile = fopen(filename,"w");

      sprintf(filename,"unit_test/C2P_%s_%s_prims.asc",con2prim_test_names[which_routine], suffix);
      FILE* prims_out = fopen(filename,"w");

      sprintf(filename,"unit_test/C2P_%s_%s_cons.asc",con2prim_test_names[which_routine], suffix);
      FILE* cons_out = fopen(filename,"w");

      sprintf(filename,"unit_test/C2P_%s_%s_Tmunu.asc",con2prim_test_names[which_routine], suffix);
      FILE* Tmunu_out = fopen(filename,"w");

      srand(0);

      for(int i=0;i<4/*npoints*/;i++) { // Density loop
        double xrho  = exp(lrmin + dlr*i);
        double P_cold = 0.0;
        double xeps  = 0.0;
        //compute_P_cold(&eos, xrho, &P_cold);
        compute_P_cold_and_eps_cold(&eos, xrho, &P_cold, &xeps);

        // Compute the pressure step size
        const double lpmin        = log(1.0e-30);//-P_cold);
        const double lpmax        = log(10.0*P_cold);
        const double dlp          = (lpmax - lpmin)/(npoints-1);
//      for(int j=0;j<npoints;j++) { // Temperature loop
        for(int j=0;j<npoints;j++) { // Pressure loop
          // Start by setting the prims (rho,Ye,T,P,eps)
//          double xtemp = exp(ltmin + dlt*j);
//          double xye   = Ye_test;
          double xpress  = exp(lpmin + dlp*j);
//          WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xpress,&xeps );

          // Velocity magnitude is derived from the W_test set before
          const double v = sqrt(1.0-1.0/(W_test*W_test));

          double phi, gxx, gyy, gzz, gxy, gxz, gyz, lapse, betax, betay, betaz;
          bool randmet = true;
          if(randmet) {
            gyy = 1.0 + randf(0.0,1.0e-1);
            gzz = 1.0 + randf(0.0,1.0e-1);
            gxy = 0*randf(-1.0e-1,1.0e-1);
            gxz = 0*randf(-1.0e-1,1.0e-1);
            gyz = 0*randf(-1.0e-1,1.0e-1);
            phi = randf(0.0,2.0);
            gxx = ( 12.0*phi - gxy * (gyz*gxz - gxy*gzz) - gxz * (gxy*gyz - gyy*gxz) )/(gyy*gzz - gyz*gyz);
            lapse = 1.0; //randf(1.0e-10,1.0);
            betax = 0.0; //v*randf(-1.0,1.0);
            betay = 0.0; //sqrt(v*v - betax*betax)*randf(0.0,1.0);
            betaz = 0.0; //sqrt(v*v - betax*betax - betay*betay);
          } else {
            gxx = 1.0;
            gyy = 1.0;
            gzz = 1.0;
            gxy = 0.0;
            gxz = 0.0;
            gyz = 0.0;
            lapse = 1.0;
            betax = 0.0;
            betay = 0.0;
            betaz = 0.0;
          }

          // Now set the velocities
          // Then, we set a maximum vy = vx/2. Solving the general expression
          // gxx vx^2 + gyy vy^2 + 2gxy vx vy < 1 in this case gives
          // => vx < ( gxx + gxy + gyy/4 )^(-1/2)
          // which is the bounds for the randf of vx
          const double vbnds = v/sqrt( gxx + gxy + gyy/4.0 );
          const double vx = randf(-vbnds, vbnds);
          const double vy = randf(-vbnds/2.0, vbnds/2.0);
          // Now, we need a vz which is consistent with the metric, vx,
          // and vy. We do this by solving the quadratic equation that
          // arises from the v^2 equation:
          // 0 = gzz vz^2 + vz*(2gxz vx + 2gyz vy)
          //   + gxx vx^2 + gyy vy^2 + 2gxy vx vy - v^2
          double vz=0;

          const double a = gzz;
          const double b = 2.0*( gxz*vx + gyz*vy );
          const double c = gxx*vx*vx + gyy*vy*vy + 2.0*gxy*vx*vy - v*v;
          const double term1 = -b/(2.0*a);
          const double term2 = sqrt(b*b - 4*a*c)/(2.0*a);

          const double s1 = term1 + term2;
          const double s2 = term1 - term2;
          if(s1<1.0 && s1>-1.0) {
            vz = s1;
          } else if (s2<1.0 && s2>-1.0) {
            vz = s2;
          } else {
            printf("Initial velocity has unphysical parameters at loop indices %d %d!", i, j);
          }

          // Finally, the magnetic fields. We'll set them aligned
          // with the velocities, for simplicity.
          const double Bhatx = vx/v;
          const double Bhaty = vy/v;
          const double Bhatz = vz/v;
          const double B     = sqrt(2.0*pow(10.0,logPmoP)*xpress);
          const double Bx    = -Bhatx * B;
          const double By    = -Bhaty * B;
          const double Bz    = -Bhatz * B;

          // Store the metric randomized values into the structs
          metric_quantities metric;
          initialize_metric(lapse,
                            gxx, gxy, gxz,
                            gyy, gyz, gzz,
                            betax, betay, betaz,
                            &metric);

          conservative_quantities cons, cons_orig, cons_undens; // Not initialized because it will be filled based on primitive data.

          primitive_quantities prims, prims_orig, prims_guess;
          initialize_primitives(xrho, xpress, xeps,
                                vx, vy, vz, Bx, By, Bz,
                                poison, poison, poison, // entropy, Y_e=xye, temp=xtemp
                                &prims);
          prims.print = true;

          // Define the stress_energy struct
          stress_energy Tmunu, Tmunu_orig;

          // Then set the conservative variables array
printf("Initial prims %.16e\n %.16e\n %.16e\n %.16e\n %.16e\n", prims.rho, prims.press, prims.vx, prims.vy, prims.vz);
          double tmp_var = (metric.adm_gxx * SQR(prims.vx + metric.betax) +
                            2.0*metric.adm_gxy*(prims.vx + metric.betax)*(prims.vy + metric.betay) +
                            2.0*metric.adm_gxz*(prims.vx + metric.betax)*(prims.vz + metric.betaz) +
                            metric.adm_gyy * SQR(prims.vy + metric.betay) +
                            2.0*metric.adm_gyz*(prims.vy + metric.betay)*(prims.vz + metric.betaz) +
                            metric.adm_gzz * SQR(prims.vz + metric.betaz) );
          double u0 = 1.0/sqrt( SQR(metric.lapse) - tmp_var );

          compute_conservs_and_Tmunu(&params, &eos, &metric, &prims, u0, &cons, &Tmunu);
//printf("initial cons %.16e\n %.16e\n %.16e\n %.16e\n %.16e\n", cons.rho, cons.tau, cons.S_x, cons.S_y, cons.S_z);

          // Store original prims
          prims_orig = prims;
          cons_orig = cons;
          Tmunu_orig = Tmunu;

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
          }

          int check = 0;
          if(cons.rho > 0.0) {

            //This applies the inequality (or "Faber") fixes on the conservatives
            if(eos.eos_type == 0) //Hybrid-only
              apply_inequality_fixes(&params, &eos, &metric, &prims, &cons, &diagnostics);
//printf("inequalities %.16e\n %.16e\n %.16e\n %.16e\n %.16e\n", cons.rho, cons.tau, cons.S_x, cons.S_y, cons.S_z);

            // The Con2Prim routines require the undensitized variables, but IGM evolves the densitized variables.
            undensitize_conservatives(&metric, &cons, &cons_undens);

            // The con2prim routines require primitive guesses in order to perform
            // the recovery. In IllinoisGRMHD, we do not keep track of the primitives
            // in between time steps, and therefore our guesses are *not* the
            // values of the primitives in the previous time level. Instead, we provide
            // guesses based on the conservative variables.

            /************* Conservative-to-primitive recovery ************/

            for(int which_guess=1;which_guess<=2;which_guess++) {
              guess_primitives(&eos, which_guess, &metric, &prims, &cons, &prims_guess);

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
              check = C2P_Select_Hybrid_Method(&params, &eos, con2prim_test_keys[which_routine], &metric, &cons_undens, &prims_guess, &diagnostics);
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
//printf("-rho prims %.16e\n %.16e\n %.16e\n %.16e\n %.16e\n", prims.rho, prims.press, prims.vx, prims.vy, prims.vz);
            diagnostics.failure_checker+=1;
            reset_prims_to_atmosphere(&params, &eos, &metric, &prims, &diagnostics);
 printf("-rho_*: reset to atm\n");
            diagnostics.rho_star_fix_applied++;
          } // if rho_star > 0

          //--------------------------------------------------
          //---------- Primitive recovery completed ----------
          //--------------------------------------------------
          // Enforce limits on primitive variables and recompute conservatives.
//printf("C2P prims %.16e\n %.16e\n %.16e\n %.16e\n %.16e\n", prims.rho, prims.press, prims.vx, prims.vy, prims.vz);
          enforce_primitive_limits_and_output_u0(&params, &eos, &metric, &prims, &u0, &diagnostics);
          compute_conservs_and_Tmunu(&params, &eos, &metric, &prims, u0, &cons, &Tmunu);
//printf("enforced prims %.16e\n %.16e\n %.16e\n %.16e\n %.16e\n", prims.rho, prims.press, prims.vx, prims.vy, prims.vz);

          primitive_quantities prims_error;
          conservative_quantities cons_error;
          stress_energy Tmunu_error;

          double accumulated_error = 0.0;

          conservative_error(&cons_orig, &cons, &cons_error);
          primitive_error(&prims_orig, &prims, &prims_error);
          stress_energy_error(&Tmunu_orig, &Tmunu, &Tmunu_error);
          if( check != 0 ) {
            failures++;
            fprintf(outfile,"Recovery FAILED!\n");
            printf("Recovery FAILED!\n");
            accumulated_error = 1e300;
          } else {
            fprintf(outfile, "Recovery SUCCEEDED FOR POINT %d_%d!\n", i,j);
            printf("Recovery SUCCEEDED FOR POINT %d_%d!\n", i,j);

//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "rho_b", prims_orig.rho, prims.rho, prims_error.rho);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "press", prims_orig.press, prims.press, prims_error.press);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "eps", prims_orig.eps, prims.eps, prims_error.eps);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "vx", prims_orig.vx, prims.vx, prims_error.vx);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "vy", prims_orig.vy, prims.vy, prims_error.vy);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "vz", prims_orig.vz, prims.vz, prims_error.vz);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "Bx", prims_orig.Bx, prims.Bx, prims_error.Bx);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "By", prims_orig.By, prims.By, prims_error.By);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "Bz", prims_orig.Bz, prims.Bz, prims_error.Bz);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "entropy", prims_orig.entropy, prims.entropy, prims_error.entropy);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "Y_e", prims_orig.Y_e, prims.Y_e, prims_error.Y_e);
//            fprintf(outfile, "Relative error for prim %s: %.3e (%e -> %e)\n", "temp", prims_orig.temp, prims.temp, prims_error.temp);
            fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "rho_b", prims_orig.rho, prims.rho, prims_error.rho);
            fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "press", prims_orig.press, prims.press, prims_error.press);
            fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "vx", prims_orig.vx, prims.vx, prims_error.vx);
            fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "vy", prims_orig.vy, prims.vy, prims_error.vy);
            fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "vz", prims_orig.vz, prims.vz, prims_error.vz);
//            fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "Bx", prims_orig.Bx, prims.Bx, prims_error.Bx);
//            fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "By", prims_orig.By, prims.By, prims_error.By);
//            fprintf(outfile, "Relative error for prim %s: (%e -> %e) %.16e\n", "Bz", prims_orig.Bz, prims.Bz, prims_error.Bz);

            fprintf(prims_out, "%d %d %.16e %.16e %.16e %.16e %.16e %.16e "
                               "%.16e %.16e %.16e %.16e %.16e %.16e "
                               "%.16e %.16e %.16e\n",
                               i, j, prims_orig.rho, prims.rho, prims_error.rho,
                               prims_orig.press, prims.press, prims_error.press,
                               prims_orig.vx, prims.vx, prims_error.vx,
                               prims_orig.vy, prims.vy, prims_error.vy,
                               prims_orig.vz, prims.vz, prims_error.vz);

            fprintf(cons_out, "%d %.16e %.16e %.16e %.16e %.16e %.16e "
                              "%.16e %.16e %.16e %.16e %.16e %.16e "
                              "%.16e %.16e %.16e\n",
                               i, cons_orig.rho, cons.rho, cons_error.rho,
                               cons_orig.tau, cons.tau, cons_error.tau,
                               cons_orig.S_x, cons.S_x, cons_error.S_x,
                               cons_orig.S_y, cons.S_y, cons_error.S_y,
                               cons_orig.S_z, cons.S_z, cons_error.S_z);

            fprintf(Tmunu_out, "%d %.16e %.16e %.16e %.16e %.16e %.16e "
                               "%.16e %.16e %.16e %.16e %.16e %.16e "
                               "%.16e %.16e %.16e %.16e %.16e %.16e "
                               "%.16e %.16e %.16e %.16e %.16e %.16e "
                               "%.16e %.16e %.16e %.16e %.16e %.16e\n",
                               i, Tmunu_orig.Ttt, Tmunu.Ttt, Tmunu_error.Ttt,
                               Tmunu_orig.Ttx, Tmunu.Ttx, Tmunu_error.Ttx,
                               Tmunu_orig.Tty, Tmunu.Tty, Tmunu_error.Tty,
                               Tmunu_orig.Ttz, Tmunu.Ttz, Tmunu_error.Ttz,
                               Tmunu_orig.Txx, Tmunu.Txx, Tmunu_error.Txx,
                               Tmunu_orig.Txy, Tmunu.Txy, Tmunu_error.Txy,
                               Tmunu_orig.Txz, Tmunu.Txz, Tmunu_error.Txz,
                               Tmunu_orig.Tyy, Tmunu.Tyy, Tmunu_error.Tyy,
                               Tmunu_orig.Tyz, Tmunu.Tyz, Tmunu_error.Tyz,
                               Tmunu_orig.Tzz, Tmunu.Tzz, Tmunu_error.Tzz);

            accumulated_error = prims_error.rho + prims_error.press /*+ prims_error.eps*/ + prims_error.vx + prims_error.vy + prims_error.vz
                              + prims_error.Bx + prims_error.By + prims_error.Bz /*+ prims_error.entropy + prims_error.Y_e + prims_error.temp*/;
          }

          fprintf(outfile, "Total accumulated error    : %e\n",accumulated_error);
          fprintf(outfile,"%e %e %e\n\n",log10(prims_orig.rho), log10(prims_orig.press), log10(MAX(accumulated_error,1e-16)));
        } // Pressure loop
//        } // Temperature loop
      } // Density loop
    fclose(outfile);
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
