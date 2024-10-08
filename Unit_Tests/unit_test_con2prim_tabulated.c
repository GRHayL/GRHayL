#include <string.h>
#include "ghl_unit_tests.h"

const char *
get_routine_string(const ghl_con2prim_method_t key) {
  switch(key) {
    case Palenzuela1D:
      return "Palenzuela1D";
    case Palenzuela1D_entropy:
      return "Palenzuela1D_entropy";
    case Newman1D:
      return "Newman1D";
    case Newman1D_entropy:
      return "Newman1D_entropy";
    default:
      ghl_error("Unsupported key %d\n", key);
      return NULL;
  }
}

void
generate_test_data(
    const ghl_parameters *restrict params,
    const ghl_eos_parameters *restrict eos ) {

  const char *routine = get_routine_string(params->main_routine);
  ghl_info("Beginning data generation for %s\n", routine);

  const int npoints = 32; // must be > 1

  for(int perturb=0;perturb<=1;perturb++) {
    for(int vars_key=0;vars_key<=1;vars_key++) {

      const char *perturb_string = perturb ? "perturbed" : "unperturbed";
      const char *vars_string    = vars_key ? "Pmag_vs_Wm1" : "rho_vs_T";

      ghl_info("  Generating %s data for test %s\n", perturb_string, vars_string);
      char filename[256];
      sprintf(filename, "con2prim_tabulated_%s_%s_%s.bin", routine, vars_string, perturb_string);
      FILE *fp = fopen(filename, "wb");

      // Write number of points to file
      fwrite(&npoints, sizeof(int), 1, fp);

      double test_Y_e=0.0/0.0, test_rho_min=0.0/0.0, test_rho_max=0.0/0.0, test_T_min=0.0/0.0, test_T_max=0.0/0.0;
      double test_W=0.0/0.0, test_logPmoP=0.0/0.0, lrmin=0.0/0.0, lrmax=0.0/0.0, dlr=0.0/0.0, ltmin=0.0/0.0, ltmax=0.0/0.0;
      double dlt=0.0/0.0, test_rho=0.0/0.0, test_T=0.0/0.0, test_logWm1_min=0.0/0.0, test_logWm1_max=0.0/0.0;
      double dlogWm1=0.0/0.0, test_logPmoP_min=0.0/0.0, test_logPmoP_max=0.0/0.0, dlogPmoP=0.0/0.0;

      if( vars_key ) {
        // These variables are used in the Pmag vs W test
        test_Y_e         = 0.1;
        test_rho         = 1.6192159535484857e-07;
        test_T           = 5;
        test_logWm1_min  = -5.5;
        test_logWm1_max  = +1.5;
        dlogWm1          = (test_logWm1_max - test_logWm1_min)/(npoints-1);
        test_logPmoP_min = -5;
        test_logPmoP_max = +9;
        dlogPmoP         = (test_logPmoP_max - test_logPmoP_min)/(npoints-1);
        if( perturb ) {
          test_Y_e         *= (1+randf(-1,1)*1e-14);
          test_rho         *= (1+randf(-1,1)*1e-14);
          test_T           *= (1+randf(-1,1)*1e-14);
          test_logWm1_min  *= (1+randf(-1,1)*1e-14);
          test_logWm1_max  *= (1+randf(-1,1)*1e-14);
          dlogWm1          *= (1+randf(-1,1)*1e-14);
          test_logPmoP_min *= (1+randf(-1,1)*1e-14);
          test_logPmoP_max *= (1+randf(-1,1)*1e-14);
          dlogPmoP         *= (1+randf(-1,1)*1e-14);
        }
      }
      else {
        // These variables are used in the P vs T test
        test_Y_e     = 0.1;
        test_rho_min = 1e-12;
        test_rho_max = 1e-3;
        test_T_min   = 1e-2;
        test_T_max   = 1e+2;
        test_W       = 2;
        test_logPmoP = -5.0;
        lrmin        = log10(test_rho_min);
        lrmax        = log10(test_rho_max);
        dlr          = (lrmax - lrmin)/(npoints-1);
        ltmin        = log10(test_T_min);
        ltmax        = log10(test_T_max);
        dlt          = (ltmax - ltmin)/(npoints-1);
        if( perturb ) {
          test_Y_e     *= (1+randf(-1,1)*1e-14);
          test_rho_min *= (1+randf(-1,1)*1e-14);
          test_rho_max *= (1+randf(-1,1)*1e-14);
          test_T_min   *= (1+randf(-1,1)*1e-14);
          test_T_max   *= (1+randf(-1,1)*1e-14);
          test_W       *= (1+randf(-1,1)*1e-14);
          test_logPmoP *= (1+randf(-1,1)*1e-14);
          lrmin        *= (1+randf(-1,1)*1e-14);
          lrmax        *= (1+randf(-1,1)*1e-14);
          dlr          *= (1+randf(-1,1)*1e-14);
          ltmin        *= (1+randf(-1,1)*1e-14);
          ltmax        *= (1+randf(-1,1)*1e-14);
          dlt          *= (1+randf(-1,1)*1e-14);
        }
      }

      // Set metric quantities to Minkowski
      ghl_metric_quantities ADM_metric;
      ghl_initialize_metric(1,
                        0, 0, 0,
                        1, 0, 0,
                        1, 0, 1,
                        &ADM_metric);

      ghl_ADM_aux_quantities metric_aux;
      ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

      srand(0);
      for(int i=0;i<npoints;i++) {
        for(int j=0;j<npoints;j++) {

          ghl_con2prim_diagnostics diagnostics;
          ghl_initialize_diagnostics(&diagnostics);

          double xrho, xtemp, xlogWm1, xW, xlogPmoP;
          if( vars_key ) {
            // Pmag vs W test
            xrho     = test_rho;
            xtemp    = test_T;
            xlogWm1  = test_logWm1_min + i*dlogWm1;
            xW       = pow(10.0, xlogWm1)+1;
            xlogPmoP = test_logPmoP_min + j*dlogPmoP;
          }
          else {
            // rho vs T test
            xrho     = pow(10.0, lrmin + dlr*i);
            xtemp    = pow(10.0, ltmin + dlt*j);
            xW       = test_W;
            xlogWm1  = log10(xW-1);
            xlogPmoP = test_logPmoP;
          }
          double xye   = test_Y_e;
          double xprs  = 0.0;
          double xeps  = 0.0;
          double xent  = 0.0;
          ghl_tabulated_compute_P_eps_S_from_T( eos, xrho, xye, xtemp, &xprs, &xeps, &xent );

          const double v     = sqrt(1.0-1.0/(xW*xW));
          const double vx    = v*((double)rand())/((double)RAND_MAX);
          const double vy    = sqrt(v*v - vx*vx)*((double)rand())/((double)RAND_MAX);
          const double vz    = sqrt(v*v - vx*vx - vy*vy);
          const double Bhatx = vx/v;
          const double Bhaty = vy/v;
          const double Bhatz = vz/v;
          const double B     = sqrt(2.0*pow(10.0,xlogPmoP)*xprs);
          const double Bx    = -Bhatx * B;
          const double By    = -Bhaty * B;
          const double Bz    = -Bhatz * B;

          // Set primitive quantities
          ghl_primitive_quantities prims_orig;
          ghl_initialize_primitives(xrho, xprs, xeps,
                                vx, vy, vz,
                                Bx, By, Bz,
                                xent, xye, xtemp,
                                &prims_orig);
          ghl_error_codes_t error = ghl_limit_v_and_compute_u0(params, &ADM_metric, &prims_orig, &diagnostics.speed_limited);
	  if(error)
            ghl_read_error_codes(error);

          // Set prim guesses
          ghl_primitive_quantities prims;
          ghl_initialize_primitives(0.0/0.0, 0.0/0.0, 0.0/0.0,
                                0.0/0.0, 0.0/0.0, 0.0/0.0,
                                Bx, By, Bz,
                                0.0/0.0, 0.0/0.0, 0.0/0.0,
                                &prims);
          prims.u0 = prims_orig.u0;

          // Compute conserved variables and Tmunu
          ghl_conservative_quantities cons;
          __attribute__((unused)) ghl_stress_energy dummy;
          ghl_compute_conservs_and_Tmunu(&ADM_metric, &metric_aux, &prims_orig, &cons, &dummy);

          // Undensitize the conserved variables
          ghl_conservative_quantities cons_undens;
          ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

          // Now perform the con2prim
          if( ghl_con2prim_tabulated_multi_method(params, eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics) )
            ghl_error("Con2Prim failed\n");

          prims.vU[0] = prims.vU[0]/prims.u0;
          prims.vU[1] = prims.vU[1]/prims.u0;
          prims.vU[2] = prims.vU[2]/prims.u0;

          // Write input and output primitives
          if( !perturb )
            fwrite(&prims_orig, sizeof(ghl_primitive_quantities), 1, fp);
          fwrite(&prims     , sizeof(ghl_primitive_quantities), 1, fp);
        }
      }
      fclose(fp);
    }
  }
}

void
run_unit_test(
    const ghl_parameters *restrict params,
    const ghl_eos_parameters *restrict eos ) {

  const char *routine = get_routine_string(params->main_routine);
  ghl_info("Beginning unit test for %s\n", routine);

  for(int vars_key=0;vars_key<=1;vars_key++) {

    const char *vars_string = vars_key ? "Pmag_vs_Wm1" : "rho_vs_T";
    ghl_info("  Running test %s\n", vars_string);
    char filename[256];
    sprintf(filename, "con2prim_tabulated_%s_%s_unperturbed.bin", routine, vars_string);
    FILE *fp_unpert = fopen(filename, "rb");
    sprintf(filename, "con2prim_tabulated_%s_%s_perturbed.bin", routine, vars_string);
    FILE *fp_pert = fopen(filename, "rb");

    int n1, n2;
    if(fread(&n1, sizeof(int), 1, fp_unpert) != 1 ) ghl_error("Failed to read from file\n");;
    if(fread(&n2, sizeof(int), 1, fp_pert  ) != 1 ) ghl_error("Failed to read from file\n");;
    if( n1 != n2 )
      ghl_error("Problem reading input data files (%d != %d)\n", n1, n2);

    const int npoints = n1;
    ghl_metric_quantities ADM_metric;
    ghl_initialize_metric(1,
                      0, 0, 0,
                      1, 0, 0,
                      1, 0, 1,
                      &ADM_metric);

    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

    for(int i=0;i<npoints;i++) {
      for(int j=0;j<npoints;j++) {

        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

        // Read input primitives from unperturbed data file
        ghl_primitive_quantities prims;
        if( fread(&prims, sizeof(ghl_primitive_quantities), 1, fp_unpert) != 1 )
          ghl_error("Failed to read input primitives from file\n");

        // Compute conserved variables and Tmunu
        ghl_conservative_quantities cons;
        __attribute__((unused)) ghl_stress_energy dummy;
        ghl_compute_conservs_and_Tmunu(&ADM_metric, &metric_aux, &prims, &cons, &dummy);

        // Undensitize the conserved variables
        ghl_conservative_quantities cons_undens;
        ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);

        // Now perform the con2prim
        if( ghl_con2prim_tabulated_multi_method(params, eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics) )
          ghl_error("Con2Prim failed\n");

        prims.vU[0] = prims.vU[0]/prims.u0;
        prims.vU[1] = prims.vU[1]/prims.u0;
        prims.vU[2] = prims.vU[2]/prims.u0;

        // Read unperturbed and perturbed results from file
        ghl_primitive_quantities prims_trusted, prims_pert;
        if( fread(&prims_trusted, sizeof(ghl_primitive_quantities), 1, fp_unpert) != 1 )
          ghl_error("Failed to read trusted primitives from file\n");
        if( fread(&prims_pert   , sizeof(ghl_primitive_quantities), 1, fp_pert) != 1 )
          ghl_error("Failed to read perturbed primitives from file\n");

        // Validate results
        ghl_pert_test_fail(prims_trusted.rho        , prims.rho        , prims_pert.rho        );
        ghl_pert_test_fail(prims_trusted.Y_e        , prims.Y_e        , prims_pert.Y_e        );
        ghl_pert_test_fail(prims_trusted.temperature, prims.temperature, prims_pert.temperature);
        ghl_pert_test_fail(prims_trusted.press      , prims.press      , prims_pert.press      );
        ghl_pert_test_fail(prims_trusted.eps        , prims.eps        , prims_pert.eps        );
        ghl_pert_test_fail(prims_trusted.vU[0]         , prims.vU[0]         , prims_pert.vU[0]);
        ghl_pert_test_fail(prims_trusted.vU[1]         , prims.vU[1]         , prims_pert.vU[1]);
        ghl_pert_test_fail(prims_trusted.vU[2]         , prims.vU[2]         , prims_pert.vU[2]);
      }
    }

    fclose(fp_unpert);
    fclose(fp_pert);
  }
}

int main(int argc, char **argv) {

  if( argc != 3 ) {
    ghl_info("Usage: %s <EOS table path> <test key>\n", argv[0]);
    ghl_info("Available test keys:\n");
    ghl_info("  0 : Generate data\n");
    ghl_info("  1 : Run unit test\n");
    exit(1);
  }

  const char *tablepath = argv[1];
  const int test_key    = atoi(argv[2]);

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const ghl_con2prim_method_t main_routine       = None;
  const ghl_con2prim_method_t backup_routines[3] = {None,None,None};
  const bool calc_prims_guess                = true;
  const bool evolve_entropy                  = false;
  const bool evolve_temperature              = true;
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
  ghl_parameters params;
  ghl_initialize_params(
        main_routine,
        backup_routines,
        evolve_entropy,
        evolve_temperature,
        calc_prims_guess,
        Psi6threshold,
        W_max,
        Lorenz_damping_factor,
        &params);

  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(
        tablepath,
        rho_b_atm, rho_b_min, rho_b_max,
        Y_e_atm, Y_e_min, Y_e_max,
        T_atm, T_min, T_max, &eos);
  eos.root_finding_precision=1e-10;

  if( test_key ) {
    params.main_routine = Palenzuela1D;     run_unit_test(&params, &eos);
    params.main_routine = Newman1D_entropy; run_unit_test(&params, &eos);
    // The routines below fail for high temperatures/magnetizations,
    // respectively. We use the standard Palenzuela routine as a backup.
    // params.backup_routine[0] = Palenzuela1D;
    params.backup_routine[0] = Palenzuela1D;
    params.main_routine = Newman1D;             run_unit_test(&params, &eos);
    params.main_routine = Palenzuela1D_entropy; run_unit_test(&params, &eos);
    ghl_info("All tests succeeded\n");
  }
  else {
    params.main_routine = Palenzuela1D;     generate_test_data(&params, &eos);
    params.main_routine = Newman1D_entropy; generate_test_data(&params, &eos);
    // The routines below fail for high temperatures/magnetizations,
    // respectively. We use the standard Palenzuela routine as a backup.
    params.backup_routine[0] = Palenzuela1D;
    params.main_routine = Newman1D;             generate_test_data(&params, &eos);
    params.main_routine = Palenzuela1D_entropy; generate_test_data(&params, &eos);
  }

  ghl_tabulated_free_memory(&eos);

  return 0;
}
