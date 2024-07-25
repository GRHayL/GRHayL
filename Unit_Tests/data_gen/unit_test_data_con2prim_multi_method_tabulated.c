#include "unit_tests.h"

void write_to_file(
      double **vars,
      const int arraylength,
      const int nvars,
      FILE *restrict outfile) {
  for(int i = 0; i < nvars; i++) {
    fwrite(vars[i], sizeof(double), arraylength, outfile);
  }
}

int main(int argc, char **argv) {

  if(argc != 2) {
    fprintf(stderr, "Usage is: ./data_gen <eos table>\n");
    exit(1);
  }

  const char *tablepath = argv[1];
  {
    FILE *fp = fopen(tablepath, "r");
    if(!fp) {
      fprintf(stderr, "Could not open table file %s\n", tablepath);
      exit(2);
    }
  }

  // These variables set up the tested range of values and number of sampling points.
  // number of sampling points in density and pressure
  const int nrho = 300;
  const int nt = 300;
  const int nye = 50;
  const int sampling = nrho * nt * nye;

  // number of additional points for ensuring we hit all logic branches
  const int ineq_edge_cases = 0;

  const int arraylength = sampling + ineq_edge_cases;

  const double test_rho_min = 1e-12; // Minimum input density
  const double test_rho_max = 1e-3;  // Maximum input density
  const double test_T_min = 1e-2;    // Minimum input temperature
  const double test_T_max = 1e+2;    // Maximum input temperature
  const double test_ye_min = 0.075;  // Minimum input electron fraction
  const double test_ye_max = 0.500;  // Maximum input electron fraction

  // Count number of routines tested
  const int num_routines = 6;
  int methods[num_routines];
  bool uses_entropy[num_routines];

  methods[0] = Palenzuela1D_entropy;
  uses_entropy[0] = true;
  methods[1] = Newman1D;
  uses_entropy[1] = false;
  methods[2] = Newman1D_entropy;
  uses_entropy[2] = true;

  // To ensure the behavior remains the same for functions after Con2Prim,
  // we explicitly set a routine to always be last.
  methods[num_routines - 1] = Palenzuela1D;
  uses_entropy[num_routines - 1] = false;

  double dummy2, dummy3;
  const double poison = 0.0 / 0.0;

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int main_routine = None;
  const int backup_routine[3] = { None, None, None };
  const bool evolve_entropy = true;
  const bool evolve_temperature = true;
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;
  const double W_max = 10.0;
  const double Lorenz_damping_factor = 0.0;

  const double rho_b_atm = 1e-12;
  const double rho_b_min = rho_b_atm;
  const double rho_b_max = 1e300; // IGM default
  const double Y_e_atm = 0.5;
  const double Y_e_min = 0.05;
  const double Y_e_max = Y_e_atm;
  const double T_atm = 1e-2;
  const double T_min = T_atm;
  const double T_max = 1e2;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(
        main_routine, backup_routine, evolve_entropy, evolve_temperature,
        calc_prims_guess, Psi6threshold, W_max, Lorenz_damping_factor, &params);

  ghl_eos_parameters eos;
  ghl_initialize_tabulated_eos_functions_and_params(
        tablepath, rho_b_atm, rho_b_min, rho_b_max, Y_e_atm, Y_e_min, Y_e_max, T_atm,
        T_min, T_max, &eos);
  eos.root_finding_precision = 1e-10;

  // rho will vary between rho_min and rho_max (uniformly in log space)
  // pressure will vary between 1.0e-30 and 10*P_cold (uniformly in log space)

  // Compute the density step size
  const double lrmin = log10(test_rho_min);
  const double lrmax = log10(test_rho_max);
  const double dlr = (lrmax - lrmin) / (nrho - 1);

  const double ltmin = log10(test_T_min);
  const double ltmax = log10(test_T_max);
  const double dlt = (ltmax - ltmin) / (nt - 1);

  const double yemin = test_ye_min;
  const double yemax = test_ye_max;
  const double dye = (yemax - yemin) / (nye - 1);

  double *lapse = (double *)malloc(sizeof(double) * arraylength);
  double *betax = (double *)malloc(sizeof(double) * arraylength);
  double *betay = (double *)malloc(sizeof(double) * arraylength);
  double *betaz = (double *)malloc(sizeof(double) * arraylength);

  double *gxx = (double *)malloc(sizeof(double) * arraylength);
  double *gxy = (double *)malloc(sizeof(double) * arraylength);
  double *gxz = (double *)malloc(sizeof(double) * arraylength);
  double *gyy = (double *)malloc(sizeof(double) * arraylength);
  double *gyz = (double *)malloc(sizeof(double) * arraylength);
  double *gzz = (double *)malloc(sizeof(double) * arraylength);

  double *Ttt = (double *)malloc(sizeof(double) * arraylength);
  double *Ttx = (double *)malloc(sizeof(double) * arraylength);
  double *Tty = (double *)malloc(sizeof(double) * arraylength);
  double *Ttz = (double *)malloc(sizeof(double) * arraylength);
  double *Txx = (double *)malloc(sizeof(double) * arraylength);
  double *Txy = (double *)malloc(sizeof(double) * arraylength);
  double *Txz = (double *)malloc(sizeof(double) * arraylength);
  double *Tyy = (double *)malloc(sizeof(double) * arraylength);
  double *Tyz = (double *)malloc(sizeof(double) * arraylength);
  double *Tzz = (double *)malloc(sizeof(double) * arraylength);

  double *rho_b = (double *)malloc(sizeof(double) * arraylength);
  double *press = (double *)malloc(sizeof(double) * arraylength);
  double *eps = (double *)malloc(sizeof(double) * arraylength);
  double *vx = (double *)malloc(sizeof(double) * arraylength);
  double *vy = (double *)malloc(sizeof(double) * arraylength);
  double *vz = (double *)malloc(sizeof(double) * arraylength);
  double *Bx = (double *)malloc(sizeof(double) * arraylength);
  double *By = (double *)malloc(sizeof(double) * arraylength);
  double *Bz = (double *)malloc(sizeof(double) * arraylength);
  double *entropy = (double *)malloc(sizeof(double) * arraylength);
  double *ye = (double *)malloc(sizeof(double) * arraylength);
  double *temperature = (double *)malloc(sizeof(double) * arraylength);

  double *rho_star = (double *)malloc(sizeof(double) * arraylength);
  double *tau = (double *)malloc(sizeof(double) * arraylength);
  double *S_x = (double *)malloc(sizeof(double) * arraylength);
  double *S_y = (double *)malloc(sizeof(double) * arraylength);
  double *S_z = (double *)malloc(sizeof(double) * arraylength);
  double *ent_star = (double *)malloc(sizeof(double) * arraylength);
  double *ye_star = (double *)malloc(sizeof(double) * arraylength);

  double *rho_b_orig = (double *)malloc(sizeof(double) * arraylength);
  double *press_orig = (double *)malloc(sizeof(double) * arraylength);
  double *eps_orig = (double *)malloc(sizeof(double) * arraylength);
  double *vx_orig = (double *)malloc(sizeof(double) * arraylength);
  double *vy_orig = (double *)malloc(sizeof(double) * arraylength);
  double *vz_orig = (double *)malloc(sizeof(double) * arraylength);
  double *ent_orig = (double *)malloc(sizeof(double) * arraylength);
  double *ye_orig = (double *)malloc(sizeof(double) * arraylength);
  double *temp_orig = (double *)malloc(sizeof(double) * arraylength);

  double *rho_star_orig = (double *)malloc(sizeof(double) * arraylength);
  double *tau_orig = (double *)malloc(sizeof(double) * arraylength);
  double *S_x_orig = (double *)malloc(sizeof(double) * arraylength);
  double *S_y_orig = (double *)malloc(sizeof(double) * arraylength);
  double *S_z_orig = (double *)malloc(sizeof(double) * arraylength);
  double *ent_star_orig = (double *)malloc(sizeof(double) * arraylength);
  double *ye_star_orig = (double *)malloc(sizeof(double) * arraylength);

  double *rho_b_pert = (double *)malloc(sizeof(double) * arraylength);
  double *press_pert = (double *)malloc(sizeof(double) * arraylength);
  double *eps_pert = (double *)malloc(sizeof(double) * arraylength);
  double *vx_pert = (double *)malloc(sizeof(double) * arraylength);
  double *vy_pert = (double *)malloc(sizeof(double) * arraylength);
  double *vz_pert = (double *)malloc(sizeof(double) * arraylength);
  double *Bx_pert = (double *)malloc(sizeof(double) * arraylength);
  double *By_pert = (double *)malloc(sizeof(double) * arraylength);
  double *Bz_pert = (double *)malloc(sizeof(double) * arraylength);
  double *ent_pert = (double *)malloc(sizeof(double) * arraylength);
  double *ye_pert = (double *)malloc(sizeof(double) * arraylength);
  double *temp_pert = (double *)malloc(sizeof(double) * arraylength);

  double *rho_star_pert = (double *)malloc(sizeof(double) * arraylength);
  double *tau_pert = (double *)malloc(sizeof(double) * arraylength);
  double *S_x_pert = (double *)malloc(sizeof(double) * arraylength);
  double *S_y_pert = (double *)malloc(sizeof(double) * arraylength);
  double *S_z_pert = (double *)malloc(sizeof(double) * arraylength);
  double *ent_star_pert = (double *)malloc(sizeof(double) * arraylength);
  double *ye_star_pert = (double *)malloc(sizeof(double) * arraylength);

  // Some functions produce/need u^0, so we need this to
  // ouput the data
  double *u0 = (double *)malloc(sizeof(double) * arraylength);

  // We also want to capture the success/failure value of
  // the C2P method
  int *c2p_check = (int *)malloc(sizeof(int) * arraylength);

  for(int k = 0; k < nye; k++) { // Electron fraction loop
    const double xye = yemin + dye * k;
    for(int j = 0; j < nt; j++) { // Temperature loop
      const double xtemp = pow(10.0, ltmin + dlt * j);
      for(int i = 0; i < nrho; i++) { // Density loop
        const double xrho = pow(10.0, lrmin + dlr * i);

        // Compute remaining primitives
        double xpress = 0.0;
        double xeps = 0.0;
        double xent = 0.0;
        ghl_tabulated_compute_P_eps_S_from_T(
              &eos, xrho, xye, xtemp, &xpress, &xeps, &xent);

        const int index = i + nrho * (j + nt * k);

        // Define the various GRHayL structs for the unit tests
        ghl_metric_quantities ADM_metric;
        ghl_primitive_quantities prims;
        ghl_conservative_quantities cons;

        // Generate random data to serve as the 'true' primitive values
        // and a randomized metric
        double local_lapse, local_betax, local_betay, local_betaz, local_gxx, local_gxy,
              local_gxz, local_gyy, local_gyz, local_gzz;
        ghl_randomize_metric(
              &local_lapse, &local_betax, &local_betay, &local_betaz, &local_gxx,
              &local_gxy, &local_gxz, &local_gyy, &local_gyz, &local_gzz);

        double local_vx, local_vy, local_vz, local_Bx, local_By, local_Bz;
        ghl_randomize_primitives(
              &eos, xrho, xpress, &local_vx, &local_vy, &local_vz, &local_Bx, &local_By,
              &local_Bz);

        ghl_initialize_metric(
              local_lapse, local_betax, local_betay, local_betaz, local_gxx, local_gxy,
              local_gxz, local_gyy, local_gyz, local_gzz, &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        ghl_initialize_primitives(
              xrho, xpress, poison, local_vx, local_vy, local_vz, local_Bx, local_By,
              local_Bz, xent, xye, xtemp, &prims);

        const int speed_limited __attribute__((unused))
        = ghl_limit_v_and_compute_u0(&params, &ADM_metric, &prims);

        // Compute conservatives based on these primitives
        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims, &cons);

        lapse[index] = ADM_metric.lapse;
        betax[index] = ADM_metric.betaU[0];
        betay[index] = ADM_metric.betaU[1];
        betaz[index] = ADM_metric.betaU[2];

        gxx[index] = ADM_metric.gammaDD[0][0];
        gxy[index] = ADM_metric.gammaDD[0][1];
        gxz[index] = ADM_metric.gammaDD[0][2];
        gyy[index] = ADM_metric.gammaDD[1][1];
        gyz[index] = ADM_metric.gammaDD[1][2];
        gzz[index] = ADM_metric.gammaDD[2][2];

        ghl_return_primitives(
              &prims, &rho_b[index], &press[index], &eps[index], &vx[index], &vy[index],
              &vz[index], &Bx[index], &By[index], &Bz[index], &entropy[index],
              &ye[index], &temperature[index]);

        ghl_return_conservatives(
              &cons, &rho_star[index], &tau[index], &S_x[index], &S_y[index],
              &S_z[index], &ent_star[index], &ye_star[index]);

        rho_b_orig[index] = rho_b[index];
        press_orig[index] = press[index];
        eps_orig[index] = eps[index];
        vx_orig[index] = vx[index];
        vy_orig[index] = vy[index];
        vz_orig[index] = vz[index];
        ent_orig[index] = entropy[index];
        ye_orig[index] = ye[index];
        temp_orig[index] = temperature[index];

        rho_star_orig[index] = rho_star[index];
        tau_orig[index] = tau[index];
        S_x_orig[index] = S_x[index];
        S_y_orig[index] = S_y[index];
        S_z_orig[index] = S_z[index];
        ent_star_orig[index] = ent_star[index];
        ye_star_orig[index] = ye_star[index];
      }
    }
  }

  /*
     Now we set up some tests for edge cases. Note
     that since the ghl_apply_conservative_limits
     edge cases require changing params, they need
     to be last.

     ghl_apply_conservative_limits branches
     Edge case 1) small B: S fix 1
     Edge case 2) \psi^6 too large
       Subcase a) \tau fix 1 & 2, don't fix S
       Subcase b) Don't fix \tau, S fix 2
  */

  ghl_metric_quantities flat_metric;
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, &flat_metric);

  if(ineq_edge_cases > 0) {
    for(int i = sampling; i < sampling + ineq_edge_cases; i++) {
      lapse[i] = 1.0;
      betax[i] = 0.0;
      betay[i] = 0.0;
      betaz[i] = 0.0;
      gxx[i] = 1.0;
      gxy[i] = 0.0;
      gxz[i] = 0.0;
      gyy[i] = 1.0;
      gyz[i] = 0.0;
      gzz[i] = 1.0;

      rho_star[i] = 1e-2;
    }

    Bx[sampling] = By[sampling] = Bz[sampling] = 1e-160;
    tau[sampling] = 1e-2;
    S_x[sampling] = S_y[sampling] = S_z[sampling]
          = 1000.0 * tau[sampling] * (tau[sampling] + 2.0 * rho_star[sampling]);

    Bx[sampling + 1] = By[sampling + 1] = Bz[sampling + 1] = 1e-2;
    // Flat space B^2 with bar rescaling
    double Bbar2
          = (Bx[sampling + 1] * Bx[sampling + 1] + By[sampling + 1] * By[sampling + 1]
             + Bz[sampling + 1] * Bz[sampling + 1])
            * SQR(ONE_OVER_SQRT_4PI);
    tau[sampling + 1] = flat_metric.sqrt_detgamma * Bbar2 / 4.0;
    S_x[sampling + 1] = S_y[sampling + 1] = S_z[sampling + 1]
          = eos.tau_atm * (eos.tau_atm + 2.0 * rho_star[sampling + 1]);

    Bx[sampling + 2] = By[sampling + 2] = Bz[sampling + 2] = 1e-2;
    // Flat space B^2 with bar rescaling
    Bbar2 = (Bx[sampling + 2] * Bx[sampling + 2] + By[sampling + 2] * By[sampling + 2]
             + Bz[sampling + 2] * Bz[sampling + 2])
            * SQR(ONE_OVER_SQRT_4PI);
    tau[sampling + 2] = 2.0 * flat_metric.sqrt_detgamma * Bbar2;
    S_x[sampling + 2] = S_y[sampling + 2] = S_z[sampling + 2]
          = 1000 * tau[sampling + 2]
            * (tau[sampling + 2] + 2.0 * rho_star[sampling + 2]);

    ent_star[sampling] = ent_star[sampling + 1] = ent_star[sampling + 2]
          = 1e-8; // no idea what it should be...

    for(int i = sampling; i < arraylength; i++) {
      rho_b_orig[i] = rho_b[i];
      press_orig[i] = press[i];
      eps_orig[i] = eps[i];
      vx_orig[i] = vx[i];
      vy_orig[i] = vy[i];
      vz_orig[i] = vz[i];
      ent_orig[i] = entropy[i];

      rho_star_orig[i] = rho_star[i];
      tau_orig[i] = tau[i];
      S_x_orig[i] = S_x[i];
      S_y_orig[i] = S_y[i];
      S_z_orig[i] = S_z[i];
      ent_star_orig[i] = ent_star[i];
    }
  }

  // Generate perturbed initial data
  for(int i = 0; i < arraylength; i++) {
    rho_b_pert[i] = rho_b[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    press_pert[i] = press[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    eps_pert[i] = eps[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    vx_pert[i] = vx[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    vy_pert[i] = vy[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    vz_pert[i] = vz[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    ent_pert[i] = entropy[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    ye_pert[i] = ye[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    temp_pert[i] = temperature[i] * (1.0 + randf(-1, 1) * 1.0e-14);

    rho_star_pert[i] = rho_star[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    tau_pert[i] = tau[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    S_x_pert[i] = S_x[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    S_y_pert[i] = S_y[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    S_z_pert[i] = S_z[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    ent_star_pert[i] = ent_star[i] * (1.0 + randf(-1, 1) * 1.0e-14);
    ye_star_pert[i] = ye_star[i] * (1.0 + randf(-1, 1) * 1.0e-14);
  }

  const int metric_length = 10;
  const int Tmunu_length = 10;
  const int prims_length = 6;
  const int cons_length = 5;

  double *metric_array[10]
        = { lapse, betax, betay, betaz, gxx, gxy, gxz, gyy, gyz, gzz };
  double *Tmunu_array[10] = { Ttt, Ttx, Tty, Ttz, Txx, Txy, Txz, Tyy, Tyz, Tzz };
  double *prims_array[6] = { rho_b, press, eps, vx, vy, vz };
  double *cons_array[5] = { rho_star, tau, S_x, S_y, S_z };

  FILE *outfile = fopen_with_check("metric_Bfield_initial_data_tabulated.bin", "wb");
  fwrite(&arraylength, sizeof(int), 1, outfile);
  write_to_file(metric_array, arraylength, metric_length, outfile);
  fwrite(Bx, sizeof(double), arraylength, outfile);
  fwrite(By, sizeof(double), arraylength, outfile);
  fwrite(Bz, sizeof(double), arraylength, outfile);
  fclose(outfile);

  char filename[100];
  char pert_suffix[7] = "";

  for(int perturb = 0; perturb < 2; perturb++) {
    for(int i = 0; i < arraylength; i++) {

    if(perturb) {
      sprintf(pert_suffix, "_pert");
      for(int i = 0; i < arraylength; i++) {
        // Cycle data
        rho_b_orig[i] = rho_b_pert[i];
        press_orig[i] = press_pert[i];
        eps_orig[i] = eps_pert[i];
        vx_orig[i] = vx_pert[i];
        vy_orig[i] = vy_pert[i];
        vz_orig[i] = vz_pert[i];
        ent_orig[i] = ent_pert[i];
        ye_orig[i] = ye_pert[i];
        temp_orig[i] = temp_pert[i];

        rho_star_orig[i] = rho_star_pert[i];
        tau_orig[i] = tau_pert[i];
        S_x_orig[i] = S_x_pert[i];
        S_y_orig[i] = S_y_pert[i];
        S_z_orig[i] = S_z_pert[i];
        ent_star_orig[i] = ent_star_pert[i];
        ye_star_orig[i] = ye_star_pert[i];
      }
    }

    sprintf(filename, "con2prim_multi_method_tabulated_output%.5s.bin", pert_suffix);
    outfile = fopen_with_check(filename, "wb");
    for(int routine = 0; routine < num_routines; routine++) {
      params.main_routine = methods[routine];
      params.evolve_entropy = uses_entropy[routine];
      int fcnt = 0;
      for(int i = 0; i < arraylength; i++) {
        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);
        ghl_metric_quantities ADM_metric;
        ghl_primitive_quantities prims;
        ghl_conservative_quantities cons, cons_undens;

        ghl_initialize_metric(
              lapse[i], betax[i], betay[i], betaz[i], gxx[i], gxy[i], gxz[i], gyy[i],
              gyz[i], gzz[i], &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        ghl_initialize_primitives(
              poison, poison, poison, poison, poison, poison, Bx[i], By[i], Bz[i],
              poison, poison, poison, &prims);

        ghl_initialize_conservatives(
              rho_star_orig[i], tau_orig[i], S_x_orig[i], S_y_orig[i], S_z_orig[i],
              ent_star[i], poison, &cons);

        ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);
        if(params.calc_prim_guess) {
          ghl_guess_primitives(&eos, &ADM_metric, &cons, &prims);
        }
        c2p_check[i] = ghl_con2prim_hybrid_select_method(
              methods[routine], &params, &eos, &ADM_metric, &metric_aux, &cons_undens,
              &prims, &diagnostics);

        ghl_return_primitives(
              &prims, &rho_b[i], &press[i], &eps[i], &vx[i], &vy[i], &vz[i], &Bx[i],
              &By[i], &Bz[i], &entropy[i], &dummy2, &dummy3);
        if(c2p_check[i]) {
          fcnt++;
        }
      }
      write_to_file(prims_array, arraylength, prims_length, outfile);
      if(params.evolve_entropy) {
        fwrite(entropy, sizeof(double), arraylength, outfile);
      }
      fwrite(c2p_check, sizeof(int), arraylength, outfile);
      printf(
            "Routine %s had %d failures out of %d points.\n",
            ghl_get_con2prim_routine_name(methods[routine]), fcnt, arraylength);
    }
    fclose(outfile);

    if(!perturb) {
      outfile = fopen_with_check(
            "enforce_primitive_limits_and_compute_u0_input.bin", "wb");
      write_to_file(prims_array, arraylength, prims_length, outfile);
      fwrite(entropy, sizeof(double), arraylength, outfile);
      fclose(outfile);
    }

    for(int i = 0; i < arraylength; i++) {
      // Cycle data
      rho_b_orig[i] = rho_b[i];
      press_orig[i] = press[i];
      eps_orig[i] = eps[i];
      vx_orig[i] = vx[i];
      vy_orig[i] = vy[i];
      vz_orig[i] = vz[i];
      ent_orig[i] = entropy[i];

      rho_star_orig[i] = rho_star[i];
      tau_orig[i] = tau[i];
      S_x_orig[i] = S_x[i];
      S_y_orig[i] = S_y[i];
      S_z_orig[i] = S_z[i];
      ent_star_orig[i] = ent_star[i];

      ghl_metric_quantities ADM_metric;
      ghl_primitive_quantities prims;

      ghl_initialize_metric(
            lapse[i], betax[i], betay[i], betaz[i], gxx[i], gxy[i], gxz[i], gyy[i],
            gyz[i], gzz[i], &ADM_metric);

      ghl_ADM_aux_quantities metric_aux;
      ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

      ghl_initialize_primitives(
            rho_b_orig[i], press_orig[i], eps_orig[i], vx_orig[i], vy_orig[i],
            vz_orig[i], Bx[i], By[i], Bz[i], ent_orig[i], poison, poison, &prims);

      const int speed_limited __attribute__((unused))
      = ghl_enforce_primitive_limits_and_compute_u0(&params, &eos, &ADM_metric, &prims);

      ghl_return_primitives(
            &prims, &rho_b[i], &press[i], &eps[i], &vx[i], &vy[i], &vz[i], &Bx[i],
            &By[i], &Bz[i], &entropy[i], &dummy2, &dummy3);
      u0[i] = prims.u0;
    }

    sprintf(
          filename, "enforce_primitive_limits_and_compute_u0_output%.5s.bin",
          pert_suffix);
    outfile = fopen_with_check(filename, "wb");
    write_to_file(prims_array, arraylength, prims_length, outfile);
    fwrite(entropy, sizeof(double), arraylength, outfile);
    fwrite(u0, sizeof(double), arraylength, outfile);
    fclose(outfile);

    for(int i = 0; i < arraylength; i++) {
      // Cycle data
      rho_b_orig[i] = rho_b[i];
      press_orig[i] = press[i];
      eps_orig[i] = eps[i];
      vx_orig[i] = vx[i];
      vy_orig[i] = vy[i];
      vz_orig[i] = vz[i];
      ent_orig[i] = entropy[i];

      rho_star_orig[i] = rho_star[i];
      tau_orig[i] = tau[i];
      S_x_orig[i] = S_x[i];
      S_y_orig[i] = S_y[i];
      S_z_orig[i] = S_z[i];
      ent_star_orig[i] = ent_star[i];

      ghl_metric_quantities ADM_metric;
      ghl_primitive_quantities prims;
      ghl_conservative_quantities cons;
      ghl_stress_energy Tmunu;

      ghl_initialize_metric(
            lapse[i], betax[i], betay[i], betaz[i], gxx[i], gxy[i], gxz[i], gyy[i],
            gyz[i], gzz[i], &ADM_metric);

      ghl_ADM_aux_quantities metric_aux;
      ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

      ghl_initialize_primitives(
            rho_b_orig[i], press_orig[i], eps_orig[i], vx_orig[i], vy_orig[i],
            vz_orig[i], Bx[i], By[i], Bz[i], ent_orig[i], poison, poison, &prims);
      prims.u0 = u0[i];

      ghl_compute_conservs_and_Tmunu(&ADM_metric, &metric_aux, &prims, &cons, &Tmunu);

      ghl_return_conservatives(
            &cons, &rho_star[i], &tau[i], &S_x[i], &S_y[i], &S_z[i], &ent_star[i],
            &dummy2);

      ghl_return_stress_energy(
            &Tmunu, &Ttt[i], &Ttx[i], &Tty[i], &Ttz[i], &Txx[i], &Txy[i], &Txz[i],
            &Tyy[i], &Tyz[i], &Tzz[i]);
    }

    sprintf(filename, "compute_conservs_and_Tmunu_output%.5s.bin", pert_suffix);
    outfile = fopen_with_check(filename, "wb");
    write_to_file(cons_array, arraylength, cons_length, outfile);
    fwrite(ent_star, sizeof(double), arraylength, outfile);
    write_to_file(Tmunu_array, arraylength, Tmunu_length, outfile);
    fclose(outfile);
  }
  printf("The file con2prim_multi_method_hybrid_input.bin contains the same data as"
         " apply_conservative_limits_output.bin. We should try to make this a symlink "
         "or something..\n");
  int ignore __attribute__((unused)) = system("cp apply_conservative_limits_output.bin "
                                              "con2prim_multi_method_hybrid_input.bin");

  printf("The file compute_conservs_and_Tmunu_input.bin contains the same data as"
         " enforce_primitive_limits_and_compute_u0_output.bin. We should try to make "
         "this a symlink or something.\n");
  ignore = system("cp enforce_primitive_limits_and_compute_u0_output.bin "
                  "compute_conservs_and_Tmunu_input.bin");

  return 0;
}
