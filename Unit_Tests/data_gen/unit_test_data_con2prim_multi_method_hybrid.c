#include "unit_tests.h"

void write_to_file(double **vars, const int arraylength, const int nvars, FILE *restrict outfile) {
  for(int i=0; i<nvars; i++)
    fwrite(vars[i], sizeof(double), arraylength, outfile);
}

int main(int argc, char **argv) {

  // These variables set up the tested range of values and number of sampling points.
  // number of sampling points in density and pressure
  const int npoints = 300;
  const int sampling = npoints*npoints;

  // number of additional points for ensuring we hit all logic branches
  const int ineq_edge_cases = 3;

  const int arraylength = sampling + ineq_edge_cases;

  double test_rho_min = 1e-12; //Minimum input density
  double test_rho_max = 1e-3; //Maximum input density

  // Count number of routines tested
  int num_routines = 2;
  int c2p_keys[num_routines];

  c2p_keys[0] = Font1D;
  //c2p_keys[1] = Palenzuela1D;

  // To ensure the behavior remains the same for functions after Con2Prim,
  // we explicitly set a routine to always be last.
  c2p_keys[num_routines-1] = Noble2D;

  double dummy1, dummy2, dummy3;
  const double poison = 0.0/0.0;

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {None,None,None};
  const bool evolve_entropy = false;
  const bool evolve_temperature = false;
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;
  const bool Cupp_fix = true;

  const int neos = 1;
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(
        None, backup_routine, evolve_entropy, evolve_temperature,
        calc_prims_guess, Psi6threshold, Cupp_fix, 0.0, &params);

  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(
        W_max, rho_b_min, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &eos);

  // rho will vary between rho_min and rho_max (uniformly in log space)
  // pressure will vary between 1.0e-30 and 10*P_cold (uniformly in log space)

  // Compute the density step size
  const double lrmin        = log(test_rho_min);
  const double lrmax        = log(test_rho_max);
  const double dlr          = (lrmax - lrmin)/(npoints-1);

  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  double *gxx = (double*) malloc(sizeof(double)*arraylength);
  double *gxy = (double*) malloc(sizeof(double)*arraylength);
  double *gxz = (double*) malloc(sizeof(double)*arraylength);
  double *gyy = (double*) malloc(sizeof(double)*arraylength);
  double *gyz = (double*) malloc(sizeof(double)*arraylength);
  double *gzz = (double*) malloc(sizeof(double)*arraylength);

  double *Ttt = (double*) malloc(sizeof(double)*arraylength);
  double *Ttx = (double*) malloc(sizeof(double)*arraylength);
  double *Tty = (double*) malloc(sizeof(double)*arraylength);
  double *Ttz = (double*) malloc(sizeof(double)*arraylength);
  double *Txx = (double*) malloc(sizeof(double)*arraylength);
  double *Txy = (double*) malloc(sizeof(double)*arraylength);
  double *Txz = (double*) malloc(sizeof(double)*arraylength);
  double *Tyy = (double*) malloc(sizeof(double)*arraylength);
  double *Tyz = (double*) malloc(sizeof(double)*arraylength);
  double *Tzz = (double*) malloc(sizeof(double)*arraylength);

  double *rho_b = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *eps = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vy = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);
  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  double *rho_star = (double*) malloc(sizeof(double)*arraylength);
  double *tau = (double*) malloc(sizeof(double)*arraylength);
  double *S_x = (double*) malloc(sizeof(double)*arraylength);
  double *S_y = (double*) malloc(sizeof(double)*arraylength);
  double *S_z = (double*) malloc(sizeof(double)*arraylength);

  double *rho_b_orig = (double*) malloc(sizeof(double)*arraylength);
  double *press_orig = (double*) malloc(sizeof(double)*arraylength);
  double *eps_orig = (double*) malloc(sizeof(double)*arraylength);
  double *vx_orig = (double*) malloc(sizeof(double)*arraylength);
  double *vy_orig = (double*) malloc(sizeof(double)*arraylength);
  double *vz_orig = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_orig = (double*) malloc(sizeof(double)*arraylength);
  double *By_orig = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_orig = (double*) malloc(sizeof(double)*arraylength);

  double *rho_star_orig = (double*) malloc(sizeof(double)*arraylength);
  double *tau_orig = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_orig = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_orig = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_orig = (double*) malloc(sizeof(double)*arraylength);

  double *rho_b_pert = (double*) malloc(sizeof(double)*arraylength);
  double *press_pert = (double*) malloc(sizeof(double)*arraylength);
  double *eps_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vy_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vz_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *By_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_pert = (double*) malloc(sizeof(double)*arraylength);

  double *rho_star_pert = (double*) malloc(sizeof(double)*arraylength);
  double *tau_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_pert = (double*) malloc(sizeof(double)*arraylength);

  // Some functions produce/need u^0, so we need this to
  // ouput the data
  double *u0 = (double*) malloc(sizeof(double)*arraylength);

  // We also want to capture the success/failure value of
  // the C2P method
  int *c2p_check = (int*) malloc(sizeof(int)*arraylength);

  for(int j=0;j<npoints;j++) { // Density loop
    const double xrho  = exp(lrmin + dlr*j);
    double P_cold = 0.0;
    double eps_cold = 0.0;
    ghl_hybrid_compute_P_cold_and_eps_cold(&eos, xrho, &P_cold, &eps_cold);

    // Compute the pressure step size
    const double lpmin        = log(1.0e-30);//-P_cold);
    const double lpmax        = log(10.0*P_cold);
    const double dlp          = (lpmax - lpmin)/(npoints-1);
    for(int i=0;i<npoints;i++) { // Pressure loop
      const int index = indexf(npoints,i,j,0);

      // Start by setting the prims (rho,Ye,T,P,eps)
      //double xtemp = exp(ltmin + dlt*i);
      //double xye   = Ye_test;
      const double xpress  = exp(lpmin + dlp*i);
      //WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xpress,&xeps );

      // Define the various GRHayL structs for the unit tests
      ghl_metric_quantities ADM_metric;
      ghl_primitive_quantities prims;
      ghl_conservative_quantities cons;

      // Generate random data to serve as the 'true' primitive values
      // and a randomized metric
      double local_lapse, local_betax, local_betay, local_betaz, local_gxx, local_gxy, local_gxz, local_gyy, local_gyz, local_gzz;
      ghl_randomize_metric(
          &local_lapse, &local_betax, &local_betay, &local_betaz,
          &local_gxx, &local_gxy, &local_gxz,
          &local_gyy, &local_gyz, &local_gzz);

      double local_eps, local_vx, local_vy, local_vz, local_Bx, local_By, local_Bz;
      ghl_randomize_primitives(
          &eos, xrho, xpress, &local_eps,
          &local_vx, &local_vy, &local_vz,
          &local_Bx, &local_By, &local_Bz);

      ghl_initialize_metric(
          local_lapse, local_betax, local_betay, local_betaz,
          local_gxx, local_gxy, local_gxz,
          local_gyy, local_gyz, local_gzz,
          &ADM_metric);

      ghl_ADM_aux_quantities metric_aux;
      ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

      ghl_initialize_primitives(
            xrho, xpress, local_eps,
            local_vx, local_vy, local_vz,
            local_Bx, local_By, local_Bz,
            poison, poison, poison, // entropy, Y_e, temp
            &prims);

      const int speed_limited __attribute__((unused)) = ghl_limit_v_and_compute_u0(&eos, &ADM_metric, &prims);

      // We need epsilon to compute the enthalpy in ghl_compute_conservs_and_Tmunu;
      // This normally happens in the ghl_enforce_primitive_limits_and_compute_u0 function
      prims.eps = eps_cold + (prims.press-P_cold)/(eos.Gamma_th-1.0)/prims.rho;

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

      ghl_return_primitives(&prims,
            &rho_b[index], &press[index], &eps[index],
            &vx[index], &vy[index], &vz[index],
            &Bx[index], &By[index], &Bz[index],
            &dummy1, &dummy2, &dummy3);

      ghl_return_conservatives(
            &cons, &rho_star[index], &tau[index],
            &S_x[index], &S_y[index], &S_z[index],
            &dummy1, &dummy2);

      rho_b_orig[index] = rho_b[index];
      press_orig[index] = press[index];
      eps_orig[index] = eps[index];
      vx_orig[index] = vx[index];
      vy_orig[index] = vy[index];
      vz_orig[index] = vz[index];
      Bx_orig[index] = Bx[index];
      By_orig[index] = By[index];
      Bz_orig[index] = Bz[index];

      rho_star_orig[index] = rho_star[index];
      tau_orig[index] = tau[index];
      S_x_orig[index] = S_x[index];
      S_y_orig[index] = S_y[index];
      S_z_orig[index] = S_z[index];
    }
  }

  // Now we set up some tests for edge cases
  // Note that since the ghl_apply_conservative_limits
  // edge cases require changing params, they need to
  // be last.

  // ghl_apply_conservative_limits branches
  // Edge case 1) small B: S fix 1
  // Edge case 2) \psi^6 too large
  //   Subcase a) \tau fix 1 & 2, don't fix S
  //   Subcase b) Don't fix \tau, S fix 2

  ghl_metric_quantities flat_metric;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        &flat_metric);

  for(int i=sampling; i<sampling+ineq_edge_cases; i++) {
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
  S_x[sampling] = S_y[sampling] = S_z[sampling] = 1000.0*tau[sampling]*(tau[sampling] + 2.0*rho_star[sampling]);

  Bx[sampling+1] = By[sampling+1] = Bz[sampling+1] = 1e-2;
  // Flat space B^2 with bar rescaling
  double Bbar2 = (Bx[sampling+1]*Bx[sampling+1] + By[sampling+1]*By[sampling+1] + Bz[sampling+1]*Bz[sampling+1])*SQR(ONE_OVER_SQRT_4PI);
  tau[sampling+1] = flat_metric.sqrt_detgamma*Bbar2/4.0;
  S_x[sampling+1] = S_y[sampling+1] = S_z[sampling+1] = eos.tau_atm*(eos.tau_atm + 2.0*rho_star[sampling+1]);

  Bx[sampling+2] = By[sampling+2] = Bz[sampling+2] = 1e-2;
  // Flat space B^2 with bar rescaling
  Bbar2 = (Bx[sampling+2]*Bx[sampling+2] + By[sampling+2]*By[sampling+2] + Bz[sampling+2]*Bz[sampling+2])*SQR(ONE_OVER_SQRT_4PI);
  tau[sampling+2] = 2.0*flat_metric.sqrt_detgamma*Bbar2;
  S_x[sampling+2] = S_y[sampling+2] = S_z[sampling+2] = 1000*tau[sampling+2]*(tau[sampling+2] + 2.0*rho_star[sampling+2]);

  for(int i=sampling; i<arraylength; i++) {
    rho_b_orig[i] = rho_b[i];
    press_orig[i] = press[i];
    eps_orig[i] = eps[i];
    vx_orig[i] = vx[i];
    vy_orig[i] = vy[i];
    vz_orig[i] = vz[i];
    Bx_orig[i] = Bx[i];
    By_orig[i] = By[i];
    Bz_orig[i] = Bz[i];

    rho_star_orig[i] = rho_star[i];
    tau_orig[i] = tau[i];
    S_x_orig[i] = S_x[i];
    S_y_orig[i] = S_y[i];
    S_z_orig[i] = S_z[i];
  }

  // Generate perturbed initial data
  for(int i=0; i<arraylength; i++) {
    rho_b_pert[i] = rho_b[i]*(1.0 + randf(-1,1)*1.0e-14);
    press_pert[i] = press[i]*(1.0 + randf(-1,1)*1.0e-14);
    eps_pert[i] = eps[i]*(1.0 + randf(-1,1)*1.0e-14);
    vx_pert[i] = vx[i]*(1.0 + randf(-1,1)*1.0e-14);
    vy_pert[i] = vy[i]*(1.0 + randf(-1,1)*1.0e-14);
    vz_pert[i] = vz[i]*(1.0 + randf(-1,1)*1.0e-14);

    rho_star_pert[i] = rho_star[i]*(1.0 + randf(-1,1)*1.0e-14);
    tau_pert[i] = tau[i]*(1.0 + randf(-1,1)*1.0e-14);
    S_x_pert[i] = S_x[i]*(1.0 + randf(-1,1)*1.0e-14);
    S_y_pert[i] = S_y[i]*(1.0 + randf(-1,1)*1.0e-14);
    S_z_pert[i] = S_z[i]*(1.0 + randf(-1,1)*1.0e-14);
  }

  const int metric_length = 10;
  const int Tmunu_length = 10;
  const int prims_length = 6;
  const int cons_length = 5;

  double *metric_array[10] = {lapse, betax, betay, betaz, gxx, gxy, gxz, gyy, gyz, gzz};
  double *Tmunu_array[10] = {Ttt, Ttx, Tty, Ttz, Txx, Txy, Txz, Tyy, Tyz, Tzz};
  double *prims_array[6] = {rho_b, press, eps, vx, vy, vz};
  double *cons_array[5] = {rho_star, tau, S_x, S_y, S_z};

  FILE* outfile = fopen_with_check("metric_Bfield_initial_data.bin", "wb");
  fwrite(&arraylength, sizeof(int), 1, outfile);
  write_to_file(metric_array, arraylength, metric_length, outfile);
  fwrite(Bx, sizeof(double), arraylength, outfile);
  fwrite(By, sizeof(double), arraylength, outfile);
  fwrite(Bz, sizeof(double), arraylength, outfile);
  fclose(outfile);

  outfile = fopen_with_check("apply_conservative_limits_input.bin", "wb");
  write_to_file(cons_array,  arraylength, cons_length,  outfile);
  fclose(outfile);

  char filename[100];
  char pert_suffix[7] = "";

  for(int perturb=0; perturb<2; perturb++) {
    for(int i=0; i<arraylength; i++) {

      // Loop over data generation again with perturbed data
      if(perturb) {
        sprintf(pert_suffix, "_pert");

        rho_b_orig[i] = rho_b_pert[i];
        press_orig[i] = press_pert[i];
        eps_orig[i] = eps_pert[i];
        vx_orig[i] = vx_pert[i];
        vy_orig[i] = vy_pert[i];
        vz_orig[i] = vz_pert[i];
        Bx_orig[i] = Bx_pert[i];
        By_orig[i] = By_pert[i];
        Bz_orig[i] = Bz_pert[i];

        rho_star_orig[i] = rho_star_pert[i];
        tau_orig[i] = tau_pert[i];
        S_x_orig[i] = S_x_pert[i];
        S_y_orig[i] = S_y_pert[i];
        S_z_orig[i] = S_z_pert[i];
      }

      ghl_con2prim_diagnostics diagnostics;
      ghl_initialize_diagnostics(&diagnostics);
      ghl_metric_quantities ADM_metric;
      ghl_primitive_quantities prims;
      ghl_conservative_quantities cons;

      ghl_initialize_metric(
            lapse[i], betax[i], betay[i], betaz[i],
            gxx[i], gxy[i], gxz[i],
            gyy[i], gyz[i], gzz[i],
            &ADM_metric);

      ghl_ADM_aux_quantities metric_aux;
      ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

      ghl_initialize_primitives(
            rho_b_orig[i], press_orig[i], eps_orig[i],
            vx_orig[i], vy_orig[i], vz_orig[i],
            Bx_orig[i], By_orig[i], Bz_orig[i],
            poison, poison, poison, &prims);

      ghl_initialize_conservatives(
            rho_star_orig[i], tau_orig[i],
            S_x_orig[i], S_y_orig[i], S_z_orig[i],
            poison, poison, &cons);

      if(i == arraylength-1 || i == arraylength-2)
        params.psi6threshold = 0.0;
      ghl_apply_conservative_limits(
            &params, &eos, &ADM_metric,
            &prims, &cons, &diagnostics);
      if(i == arraylength-1 || i == arraylength-2)
        params.psi6threshold = Psi6threshold;

      ghl_return_conservatives(
            &cons, &rho_star[i], &tau[i],
            &S_x[i], &S_y[i], &S_z[i],
            &dummy1, &dummy2);
    }

    sprintf(filename,"apply_conservative_limits_output%.5s.bin", pert_suffix);
    outfile = fopen_with_check(filename, "wb");
    write_to_file(cons_array, arraylength, cons_length, outfile);
    fclose(outfile);

    if(!perturb) {
      printf("The file con2prim_multi_method_hybrid_input.bin contains the same data as"
             " apply_conservative_limits_output.bin, so a symlink handles this duplication.\n");
    }

    for(int i=0; i<arraylength; i++) {
      // Cycle data
      rho_b_orig[i] = rho_b[i];
      press_orig[i] = press[i];
      eps_orig[i] = eps[i];
      vx_orig[i] = vx[i];
      vy_orig[i] = vy[i];
      vz_orig[i] = vz[i];
      Bx_orig[i] = Bx[i];
      By_orig[i] = By[i];
      Bz_orig[i] = Bz[i];

      rho_star_orig[i] = rho_star[i];
      tau_orig[i] = tau[i];
      S_x_orig[i] = S_x[i];
      S_y_orig[i] = S_y[i];
      S_z_orig[i] = S_z[i];
    }

    sprintf(filename,"con2prim_multi_method_hybrid_output%.5s.bin", pert_suffix);
    outfile = fopen_with_check(filename, "wb");
    for(int routine=0; routine<num_routines; routine++) {
      params.main_routine = c2p_keys[routine];
      for(int i=0; i<arraylength; i++) {
        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);
        ghl_metric_quantities ADM_metric;
        ghl_primitive_quantities prims;
        ghl_conservative_quantities cons, cons_undens;

        ghl_initialize_metric(
              lapse[i], betax[i], betay[i], betaz[i],
              gxx[i], gxy[i], gxz[i],
              gyy[i], gyz[i], gzz[i],
              &ADM_metric);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        ghl_initialize_primitives(
              rho_b_orig[i], press_orig[i], eps_orig[i],
              vx_orig[i], vy_orig[i], vz_orig[i],
              Bx_orig[i], By_orig[i], Bz_orig[i],
              poison, poison, poison, &prims);

        ghl_initialize_conservatives(
              rho_star_orig[i], tau_orig[i],
              S_x_orig[i], S_y_orig[i], S_z_orig[i],
              poison, poison, &cons);

        ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens);
        c2p_check[i] = ghl_con2prim_hybrid_multi_method(
              &params, &eos, &ADM_metric, &metric_aux,
              &cons_undens, &prims, &diagnostics);

        ghl_return_primitives(
              &prims, &rho_b[i], &press[i], &eps[i],
              &vx[i], &vy[i], &vz[i],
              &Bx[i], &By[i], &Bz[i],
              &dummy1, &dummy2, &dummy3);
      }

      write_to_file(prims_array, arraylength, prims_length, outfile);
      fwrite(c2p_check, sizeof(int), arraylength, outfile);
    }
    fclose(outfile);

    if(!perturb) {
      outfile = fopen_with_check("enforce_primitive_limits_and_compute_u0_input.bin", "wb");
      write_to_file(prims_array, arraylength, prims_length, outfile);
      fclose(outfile);
    }

    for(int i=0; i<arraylength; i++) {
      // Cycle data
      rho_b_orig[i] = rho_b[i];
      press_orig[i] = press[i];
      eps_orig[i] = eps[i];
      vx_orig[i] = vx[i];
      vy_orig[i] = vy[i];
      vz_orig[i] = vz[i];
      Bx_orig[i] = Bx[i];
      By_orig[i] = By[i];
      Bz_orig[i] = Bz[i];

      rho_star_orig[i] = rho_star[i];
      tau_orig[i] = tau[i];
      S_x_orig[i] = S_x[i];
      S_y_orig[i] = S_y[i];
      S_z_orig[i] = S_z[i];

      ghl_metric_quantities ADM_metric;
      ghl_primitive_quantities prims;

      ghl_initialize_metric(
            lapse[i], betax[i], betay[i], betaz[i],
            gxx[i], gxy[i], gxz[i],
            gyy[i], gyz[i], gzz[i],
            &ADM_metric);

      ghl_ADM_aux_quantities metric_aux;
      ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

      ghl_initialize_primitives(
            rho_b_orig[i], press_orig[i], eps_orig[i],
            vx_orig[i], vy_orig[i], vz_orig[i],
            Bx_orig[i], By_orig[i], Bz_orig[i],
            poison, poison, poison, &prims);

      const int speed_limited __attribute__((unused)) = ghl_enforce_primitive_limits_and_compute_u0(
            &params, &eos, &ADM_metric, &prims);

      ghl_return_primitives(
            &prims, &rho_b[i], &press[i], &eps[i],
            &vx[i], &vy[i], &vz[i],
            &Bx[i], &By[i], &Bz[i],
            &dummy1, &dummy2, &dummy3);
      u0[i] = prims.u0;
    }

    sprintf(filename,"enforce_primitive_limits_and_compute_u0_output%.5s.bin", pert_suffix);
    outfile = fopen_with_check(filename, "wb");
    write_to_file(prims_array, arraylength, prims_length, outfile);
    fwrite(u0, sizeof(double), arraylength, outfile);
    fclose(outfile);

    if(!perturb) {
      printf("The file compute_conservs_and_Tmunu_input.bin contains the same data as"
             " enforce_primitive_limits_and_compute_u0_output.bin, so a symlink handles this duplication.\n");
    }

    for(int i=0; i<arraylength; i++) {
      // Cycle data
      rho_b_orig[i] = rho_b[i];
      press_orig[i] = press[i];
      eps_orig[i] = eps[i];
      vx_orig[i] = vx[i];
      vy_orig[i] = vy[i];
      vz_orig[i] = vz[i];
      Bx_orig[i] = Bx[i];
      By_orig[i] = By[i];
      Bz_orig[i] = Bz[i];

      rho_star_orig[i] = rho_star[i];
      tau_orig[i] = tau[i];
      S_x_orig[i] = S_x[i];
      S_y_orig[i] = S_y[i];
      S_z_orig[i] = S_z[i];

      ghl_metric_quantities ADM_metric;
      ghl_primitive_quantities prims;
      ghl_conservative_quantities cons;
      ghl_stress_energy Tmunu;

      ghl_initialize_metric(
            lapse[i], betax[i], betay[i], betaz[i],
            gxx[i], gxy[i], gxz[i],
            gyy[i], gyz[i], gzz[i],
            &ADM_metric);

      ghl_ADM_aux_quantities metric_aux;
      ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

      ghl_initialize_primitives(
            rho_b_orig[i], press_orig[i], eps_orig[i],
            vx_orig[i], vy_orig[i], vz_orig[i],
            Bx_orig[i], By_orig[i], Bz_orig[i],
            poison, poison, poison, &prims);
      prims.u0 = u0[i];

      ghl_compute_conservs_and_Tmunu(
            &ADM_metric, &metric_aux, &prims, &cons, &Tmunu);

      ghl_return_conservatives(
            &cons, &rho_star[i], &tau[i],
            &S_x[i], &S_y[i], &S_z[i],
            &dummy1, &dummy2);

      ghl_return_stress_energy(&Tmunu,
             &Ttt[i], &Ttx[i], &Tty[i],
             &Ttz[i], &Txx[i], &Txy[i],
             &Txz[i], &Tyy[i], &Tyz[i],
             &Tzz[i]);

    }

    sprintf(filename,"compute_conservs_and_Tmunu_output%.5s.bin", pert_suffix);
    outfile = fopen_with_check(filename, "wb");
    write_to_file(cons_array, arraylength, cons_length, outfile);
    write_to_file(Tmunu_array, arraylength, Tmunu_length, outfile);
    fclose(outfile);
  }
  return 0;
}
