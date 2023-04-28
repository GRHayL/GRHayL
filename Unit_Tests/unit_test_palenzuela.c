#include "unit_tests.h"

int main(int argc, char **argv) {

  unsigned test_type = 1;
  if( argc < 2 || argc > 3 )
    grhayl_error("Correct usage is %s <EOS table path> [test type (1=rho vs T; 2=Pmag vs W]\n", argv[0]);
  else if( argc == 3 )
    test_type = atoi(argv[2]);

  char outfile_suffix[10];
  switch (test_type) {
    case 1:
      grhayl_info("Selected test: rho vs T\n");
      sprintf(outfile_suffix, "rho_vs_T");
      break;
    case 2:
      grhayl_info("Selected test: Pmag vs W\n");
      sprintf(outfile_suffix, "Pmag_vs_W");
      break;
    default:
      grhayl_error("Unknown test key: %d. Supported keys are 1 (rho vs T) and 2 (Pmag vs W).\n");
      break;
  }

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const con2prim_method_t main_routine       = Palenzuela1D_entropy;
  const con2prim_method_t backup_routines[3] = {None,None,None};
  const bool calc_prims_guess                = false;
  const bool evolve_entropy                  = false;
  const bool evolve_temperature              = true;
  const bool update_Tmunu                    = true; //IGM default
  const bool use_cupp_fix                    = true;
  const double Psi6threshold                 = 1e100; //Taken from magnetizedTOV.par
  const double Lorenz_damping_factor         = 0.0;

  const double W_max     = 100.0; //IGM default: 10
  const double rho_b_atm = 1e-12;
  const double rho_b_min = rho_b_atm;
  const double rho_b_max = 1e300; //IGM default
  const double Y_e_atm   = 0.5;
  const double Y_e_min   = 0.05;
  const double Y_e_max   = Y_e_atm;
  const double T_atm     = 1e-2;
  const double T_min     = T_atm;
  const double T_max     = 1e2;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  GRHayL_parameters params;
  initialize_GRHayL(main_routine,
                    backup_routines,
                    evolve_entropy,
                    evolve_temperature,
                    calc_prims_guess,
                    Psi6threshold,
                    update_Tmunu,
                    use_cupp_fix,
                    Lorenz_damping_factor,
                    &params);

  eos_parameters eos;
  initialize_tabulated_eos_functions_and_params(argv[1], W_max,
                                                rho_b_atm, rho_b_min, rho_b_max,
                                                Y_e_atm, Y_e_min, Y_e_max,
                                                T_atm, T_min, T_max, &eos);
  eos.root_finding_precision=1e-12;

  // Number of points in the discretization of rho and T
  const int npoints = 256;

  // const double test_rho_min = 1e-12; // IGM default
  // const double test_rho_max = 1e-3;  // IGM default
  // const double test_T_min   = 1e-2;  // IGM default
  // const double test_T_max   = 1e+2;  // IGM default

  // Compute the density step size
  // const double lrmin        = log(test_rho_min);
  // const double lrmax        = log(test_rho_max);
  // const double dlr          = npoints > 1 ? (lrmax - lrmin)/(npoints-1) : 0;
  const double rho_test = 1.6192159535484857e-07; // 1e11 g/cm3

  // Compute the temperature step size
  // const double ltmin        = log(test_T_min);
  // const double ltmax        = log(test_T_max);
  // const double dlt          = npoints > 1 ? (ltmax - ltmin)/(npoints-1) : 0;
  const double T_test = 5; // in MeV

  // Fix Y_e
  const double Ye_test = 0.1;
  // Fix W
  // const double W_test  = 2.0;
  const double logWm1_min = -5.5;
  const double logWm1_max = +1.5;
  const double dlogWm1    = (logWm1_max-logWm1_min)/(npoints-1);

  // Fix log10(Pmag/P)
  // const double logPmoP = -5.0;
  const double logPmoP_min = -5;
  const double logPmoP_max = +9;
  const double dlogPmoP    = (logPmoP_max - logPmoP_min)/(npoints-1);

  // Set metric quantities to Minkowski
  metric_quantities metric;
  initialize_metric(1, 1, 0, 0, 1, 0, 1, 0, 0, 0, &metric);

  char outfile[256];
  switch(main_routine) {
    case Palenzuela1D:
      sprintf(outfile, "palenzuela_energy_%s.asc", outfile_suffix);
      break;
    case Palenzuela1D_entropy:
      sprintf(outfile, "palenzuela_entropy_energy_%s.asc", outfile_suffix);
      break;
    case Newman1D:
      sprintf(outfile, "newman_energy_%s.asc", outfile_suffix);
      break;
    case Newman1D_entropy:
      sprintf(outfile, "newman_entropy_%s.asc", outfile_suffix);
      break;
    default:
      grhayl_error("Con2Prim routine not supported in this unit test\n");
      break;
  }

  FILE *fp = fopen(outfile, "w");

  srand(0);

  double xrho  = rho_test;
  double xtemp = T_test;
  double xye   = Ye_test;
  double xprs  = 0.0;
  double xeps  = 0.0;
  double xent  = 0.0;
  eos.tabulated_compute_P_eps_S_from_T( &eos, xrho, xye, xtemp, &xprs, &xeps, &xent );

  for(int i=0;i<npoints;i++) {
    for(int j=0;j<npoints;j++) {

      con2prim_diagnostics diagnostics;
      initialize_diagnostics(&diagnostics);

      // double xrho  = exp(lrmin + dlr*i);
      // double xtemp = exp(ltmin + dlt*j);
      // double xye   = Ye_test;
      // double xprs  = 0.0;
      // double xeps  = 0.0;
      // double xent  = 0.0;
      // eos.tabulated_compute_P_eps_S_from_T( &eos, xrho, xye, xtemp, &xprs, &xeps, &xent );

      const double logWm1 = logWm1_min+i*dlogWm1;
      const double W_test = pow(10.0,logWm1)+1;
      const double v      = sqrt(1.0-1.0/(W_test*W_test));
      const double vx     = v*((double)rand())/((double)RAND_MAX);
      const double vy     = sqrt(v*v - vx*vx)*((double)rand())/((double)RAND_MAX);
      const double vz     = sqrt(v*v - vx*vx - vy*vy);

      const double logPmoP = logPmoP_min + j*dlogPmoP;
      const double Bhatx   = vx/v;
      const double Bhaty   = vy/v;
      const double Bhatz   = vz/v;
      const double B       = sqrt(2.0*pow(10.0,logPmoP)*xprs);
      const double Bx      = -Bhatx * B;
      const double By      = -Bhaty * B;
      const double Bz      = -Bhatz * B;

      // Set primitive quantities
      primitive_quantities prims_orig;
      initialize_primitives(xrho, xprs, xeps,
                            vx, vy, vz,
                            Bx, By, Bz,
                            xent, xye, xtemp,
                            &prims_orig);
      limit_v_and_compute_u0(&eos, &metric, &prims_orig, &diagnostics.vel_limited_ptcount);

      // Set prim guesses
      primitive_quantities prims;
      initialize_primitives(0.0/0.0, 0.0/0.0, 0.0/0.0,
                            0.0/0.0, 0.0/0.0, 0.0/0.0,
                            Bx, By, Bz,
                            0.0/0.0, 0.0/0.0, xtemp*0.95,
                            &prims);
      prims.u0 = prims_orig.u0;

      // Compute conserved variables and Tmunu
      conservative_quantities cons;
      __attribute__((unused)) stress_energy dummy;
      compute_conservs_and_Tmunu(&params, &eos, &metric, &prims_orig, &cons, &dummy);

      // Undensitize the conserved variables
      conservative_quantities cons_undens;
      undensitize_conservatives(&metric, &cons, &cons_undens);

      // Now perform the con2prim
      Tabulated_Multi_Method(&params, &eos, &metric, &cons_undens, &prims, &diagnostics);

      prims.vx = prims.vx/prims.u0;
      prims.vy = prims.vy/prims.u0;
      prims.vz = prims.vz/prims.u0;

      double err, err_total=0.0;

      err = relative_error(prims.rho, prims_orig.rho);
      err_total += err;

      err = relative_error(prims.Y_e, prims_orig.Y_e);
      err_total += err;

      err = relative_error(prims.temperature, prims_orig.temperature);
      err_total += err;

      err = relative_error(prims.press, prims_orig.press);
      err_total += err;

      err = relative_error(prims.eps, prims_orig.eps);
      err_total += err;

      err = relative_error(prims.vx, prims_orig.vx);
      err_total += err;

      err = relative_error(prims.vy, prims_orig.vy);
      err_total += err;

      err = relative_error(prims.vz, prims_orig.vz);
      err_total += err;

      err = relative_error(prims.Bx, prims_orig.Bx);
      err_total += err;

      err = relative_error(prims.By, prims_orig.By);
      err_total += err;

      err = relative_error(prims.Bz, prims_orig.Bz);
      err_total += err;

      // fprintf(fp, "%.15e %.15e %.15e\n", log10(prims_orig.rho), log10(prims_orig.temperature), log10(err_total+1e-16));
      fprintf(fp, "%.15e %.15e %.15e\n", logWm1, logPmoP, log10(err_total+1e-16));
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  eos.tabulated_free_memory(&eos);

  return 0;
}
