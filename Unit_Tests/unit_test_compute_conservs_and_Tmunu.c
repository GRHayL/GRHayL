#include "unit_tests.h"

int main(int argc, char **argv) {

  FILE* infile = fopen_with_check("metric_Bfield_initial_data.bin","rb");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, infile);
  if( key != 1 || arraylength < 1)
    ghl_error("An error has occured with reading the grid size. "
                 "Please check that metric_initial_data.bin"
                 "is up-to-date with current test version.\n");

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {None,None,None};
  const bool evolve_entropy = true;
  const bool evolve_temperature = false;
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;
  const bool ignore_negative_pressure = true;
  const double W_max = 10.0;
  const double Lorenz_damping_factor = 0.0;

  const int neos = 1;
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
        None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
        Psi6threshold, ignore_negative_pressure, W_max, Lorenz_damping_factor, &params);

  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_min, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &eos);

  // Allocate memory for the metric data
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

  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  key  = fread(lapse, sizeof(double), arraylength, infile);
  key += fread(betax, sizeof(double), arraylength, infile);
  key += fread(betay, sizeof(double), arraylength, infile);
  key += fread(betaz, sizeof(double), arraylength, infile);

  key += fread(gxx, sizeof(double), arraylength, infile);
  key += fread(gxy, sizeof(double), arraylength, infile);
  key += fread(gxz, sizeof(double), arraylength, infile);
  key += fread(gyy, sizeof(double), arraylength, infile);
  key += fread(gyz, sizeof(double), arraylength, infile);
  key += fread(gzz, sizeof(double), arraylength, infile);

  key += fread(Bx, sizeof(double), arraylength, infile);
  key += fread(By, sizeof(double), arraylength, infile);
  key += fread(Bz, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*13)
    ghl_error("An error has occured with reading in metric data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the initial primitive data
  double *rho_b = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *eps = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vy = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);
  double *entropy = (double*) malloc(sizeof(double)*arraylength);

  // This function uses u^0, so we need to store it as well
  double *u0 = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("compute_conservs_and_Tmunu_input.bin","rb");
  key  = fread(rho_b, sizeof(double), arraylength, infile);
  key += fread(press, sizeof(double), arraylength, infile);
  key += fread(eps, sizeof(double), arraylength, infile);
  key += fread(vx, sizeof(double), arraylength, infile);
  key += fread(vy, sizeof(double), arraylength, infile);
  key += fread(vz, sizeof(double), arraylength, infile);
  key += fread(entropy, sizeof(double), arraylength, infile);
  key += fread(u0, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*8)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the trusted conservative and Tmunu data
  double *rho_star_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *tau_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *ent_trusted = (double*) malloc(sizeof(double)*arraylength);

  double *Ttt_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Ttx_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Tty_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Ttz_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Txx_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Txy_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Txz_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Tyy_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Tyz_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Tzz_trusted = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("compute_conservs_and_Tmunu_output.bin","rb");
  key  = fread(rho_star_trusted, sizeof(double), arraylength, infile);
  key += fread(tau_trusted, sizeof(double), arraylength, infile);
  key += fread(S_x_trusted, sizeof(double), arraylength, infile);
  key += fread(S_y_trusted, sizeof(double), arraylength, infile);
  key += fread(S_z_trusted, sizeof(double), arraylength, infile);
  key += fread(ent_trusted, sizeof(double), arraylength, infile);

  key += fread(Ttt_trusted, sizeof(double), arraylength, infile);
  key += fread(Ttx_trusted, sizeof(double), arraylength, infile);
  key += fread(Tty_trusted, sizeof(double), arraylength, infile);
  key += fread(Ttz_trusted, sizeof(double), arraylength, infile);
  key += fread(Txx_trusted, sizeof(double), arraylength, infile);
  key += fread(Txy_trusted, sizeof(double), arraylength, infile);
  key += fread(Txz_trusted, sizeof(double), arraylength, infile);
  key += fread(Tyy_trusted, sizeof(double), arraylength, infile);
  key += fread(Tyz_trusted, sizeof(double), arraylength, infile);
  key += fread(Tzz_trusted, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*16)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the perturbed conservative and Tmunu data
  double *rho_star_pert = (double*) malloc(sizeof(double)*arraylength);
  double *tau_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_pert = (double*) malloc(sizeof(double)*arraylength);
  double *ent_pert = (double*) malloc(sizeof(double)*arraylength);

  double *Ttt_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Ttx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Tty_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Ttz_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Txx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Txy_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Txz_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Tyy_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Tyz_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Tzz_pert = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("compute_conservs_and_Tmunu_output_pert.bin","rb");
  key  = fread(rho_star_pert, sizeof(double), arraylength, infile);
  key += fread(tau_pert, sizeof(double), arraylength, infile);
  key += fread(S_x_pert, sizeof(double), arraylength, infile);
  key += fread(S_y_pert, sizeof(double), arraylength, infile);
  key += fread(S_z_pert, sizeof(double), arraylength, infile);
  key += fread(ent_pert, sizeof(double), arraylength, infile);

  key += fread(Ttt_pert, sizeof(double), arraylength, infile);
  key += fread(Ttx_pert, sizeof(double), arraylength, infile);
  key += fread(Tty_pert, sizeof(double), arraylength, infile);
  key += fread(Ttz_pert, sizeof(double), arraylength, infile);
  key += fread(Txx_pert, sizeof(double), arraylength, infile);
  key += fread(Txy_pert, sizeof(double), arraylength, infile);
  key += fread(Txz_pert, sizeof(double), arraylength, infile);
  key += fread(Tyy_pert, sizeof(double), arraylength, infile);
  key += fread(Tyz_pert, sizeof(double), arraylength, infile);
  key += fread(Tzz_pert, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*16)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  const double poison = 0.0/0.0;

  for(int i=0;i<arraylength;i++) {

    // Define the various GRHayL structs for the unit tests
    ghl_con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);
    ghl_metric_quantities ADM_metric;
    ghl_primitive_quantities prims;
    ghl_conservative_quantities cons;
    ghl_stress_energy Tmunu;

    // Read initial data accompanying trusted output
    ghl_initialize_metric(
          lapse[i],
          betax[i], betay[i], betaz[i],
          gxx[i], gxy[i], gxz[i],
          gyy[i], gyz[i], gzz[i],
          &ADM_metric);

    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

    ghl_initialize_primitives(
          rho_b[i], press[i], eps[i],
          vx[i], vy[i], vz[i],
          Bx[i], By[i], Bz[i],
          entropy[i], poison, poison,
          &prims);
    prims.u0 = u0[i];

    // This computes the conservatives and stress-energy tensor from the new primitives
    ghl_compute_conservs_and_Tmunu(&ADM_metric, &metric_aux, &prims, &cons, &Tmunu);

    ghl_conservative_quantities cons_trusted, cons_pert;
    ghl_stress_energy Tmunu_trusted, Tmunu_pert;

    ghl_initialize_conservatives(
          rho_star_trusted[i], tau_trusted[i],
          S_x_trusted[i], S_y_trusted[i], S_z_trusted[i],
          ent_trusted[i], poison, &cons_trusted);

    ghl_initialize_conservatives(
          rho_star_pert[i], tau_pert[i],
          S_x_pert[i], S_y_pert[i], S_z_pert[i],
          ent_pert[i], poison, &cons_pert);

    ghl_initialize_stress_energy(
          Ttt_trusted[i], Ttx_trusted[i],
          Tty_trusted[i], Ttz_trusted[i],
          Txx_trusted[i], Txy_trusted[i],
          Txz_trusted[i], Tyy_trusted[i],
          Tyz_trusted[i], Tzz_trusted[i],
          &Tmunu_trusted);

    ghl_initialize_stress_energy(
          Ttt_pert[i], Ttx_pert[i],
          Tty_pert[i], Ttz_pert[i],
          Txx_pert[i], Txy_pert[i],
          Txz_pert[i], Tyy_pert[i],
          Tyz_pert[i], Tzz_pert[i],
          &Tmunu_pert);

    ghl_pert_test_fail_conservatives(params.evolve_entropy, &cons_trusted, &cons, &cons_pert);
    ghl_pert_test_fail_stress_energy(&Tmunu_trusted, &Tmunu, &Tmunu_pert);

    /*
       GRHayL also has standalone functions ghl_compute_conservs() and ghl_compute_TDNmunu().
       These exist for several reasons:
         1) the user only wants to compute the conservs
         2) the Tmunu variable is at a different centering than the conservatives,
            requiring some sort of interpolation which the user handles separately.
       We can easily use this test to also check these functions.
    */
    ghl_compute_conservs(&ADM_metric, &metric_aux, &prims, &cons);
    ghl_compute_TDNmunu(&ADM_metric, &metric_aux, &prims, &Tmunu);
    ghl_pert_test_fail_conservatives(params.evolve_entropy, &cons_trusted, &cons, &cons_pert);
    ghl_pert_test_fail_stress_energy(&Tmunu_trusted, &Tmunu, &Tmunu_pert);
  }
  ghl_info("ghl_compute_conservs_and_Tmunu function test has passed!\n");
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(rho_b); free(press); free(eps);
  free(vx); free(vy); free(vz);
  free(Bx); free(By); free(Bz);
  free(u0);
  free(rho_star_trusted); free(tau_trusted);
  free(S_x_trusted); free(S_y_trusted); free(S_z_trusted);
  free(rho_star_pert); free(tau_pert);
  free(S_x_pert); free(S_y_pert); free(S_z_pert);
  free(Ttt_trusted); free(Ttx_trusted); free(Tty_trusted);
  free(Ttz_trusted); free(Txx_trusted); free(Txy_trusted);
  free(Txz_trusted); free(Tyy_trusted); free(Tyz_trusted);
  free(Tzz_trusted);
  free(Ttt_pert); free(Ttx_pert); free(Tty_pert);
  free(Ttz_pert); free(Txx_pert); free(Txy_pert);
  free(Txz_pert); free(Tyy_pert); free(Tyz_pert);
  free(Tzz_pert);
}
