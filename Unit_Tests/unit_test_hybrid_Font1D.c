#include "unit_tests.h"

int main(int argc, char **argv) {

  FILE* infile = fopen_with_check("metric_initial_data.bin","rb");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, infile);
  if( key != 1 || arraylength < 1 )
    ghl_error("An error has occured with reading the grid size. "
                 "Please check that metric_initial_data.bin"
                 "is up-to-date with current test version.\n");

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int main_routine = Font1D;
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
  ghl_initialize_params(main_routine, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
                    Psi6threshold, Cupp_fix, 0.0 /*Lorenz damping factor*/, &params);

  eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(W_max,
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

  fclose(infile);
  if(key != arraylength*10)
    ghl_error("An error has occured with reading in metric data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the initial primitive data
  double *rho_b = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *eps = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vy = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);
  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for the initial conservative data
  double *rho_star = (double*) malloc(sizeof(double)*arraylength);
  double *tau = (double*) malloc(sizeof(double)*arraylength);
  double *S_x = (double*) malloc(sizeof(double)*arraylength);
  double *S_y = (double*) malloc(sizeof(double)*arraylength);
  double *S_z = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("con2prim_multi_method_hybrid_input.bin","rb");
  key  = fread(rho_b, sizeof(double), arraylength, infile);
  key += fread(press, sizeof(double), arraylength, infile);
  key += fread(eps, sizeof(double), arraylength, infile);
  key += fread(vx, sizeof(double), arraylength, infile);
  key += fread(vy, sizeof(double), arraylength, infile);
  key += fread(vz, sizeof(double), arraylength, infile);
  key += fread(Bx, sizeof(double), arraylength, infile);
  key += fread(By, sizeof(double), arraylength, infile);
  key += fread(Bz, sizeof(double), arraylength, infile);

  key += fread(rho_star, sizeof(double), arraylength, infile);
  key += fread(tau, sizeof(double), arraylength, infile);
  key += fread(S_x, sizeof(double), arraylength, infile);
  key += fread(S_y, sizeof(double), arraylength, infile);
  key += fread(S_z, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*14)
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the trusted primitive data
  double *rho_b_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *press_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vx_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vy_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vz_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *eps_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *By_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_trusted = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for the returned value of C2P routine
  int *c2p_check = (int*) malloc(sizeof(int)*arraylength);

  infile = fopen_with_check("con2prim_multi_method_hybrid_Font1D_output.bin","rb");
  key  = fread(rho_b_trusted, sizeof(double), arraylength, infile);
  key += fread(press_trusted, sizeof(double), arraylength, infile);
  key += fread(eps_trusted, sizeof(double), arraylength, infile);
  key += fread(vx_trusted, sizeof(double), arraylength, infile);
  key += fread(vy_trusted, sizeof(double), arraylength, infile);
  key += fread(vz_trusted, sizeof(double), arraylength, infile);
  key += fread(Bx_trusted, sizeof(double), arraylength, infile);
  key += fread(By_trusted, sizeof(double), arraylength, infile);
  key += fread(Bz_trusted, sizeof(double), arraylength, infile);
  key += fread(c2p_check, sizeof(int), arraylength, infile);

  fclose(infile);
  if(key != arraylength*10)
    ghl_error("An error has occured with reading in trusted data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the perturbed primitive data
  double *rho_b_pert = (double*) malloc(sizeof(double)*arraylength);
  double *press_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vy_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vz_pert = (double*) malloc(sizeof(double)*arraylength);
  double *eps_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *By_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_pert = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("con2prim_multi_method_hybrid_Font1D_output_pert.bin","rb");
  key  = fread(rho_b_pert, sizeof(double), arraylength, infile);
  key += fread(press_pert, sizeof(double), arraylength, infile);
  key += fread(eps_pert, sizeof(double), arraylength, infile);
  key += fread(vx_pert, sizeof(double), arraylength, infile);
  key += fread(vy_pert, sizeof(double), arraylength, infile);
  key += fread(vz_pert, sizeof(double), arraylength, infile);
  key += fread(Bx_pert, sizeof(double), arraylength, infile);
  key += fread(By_pert, sizeof(double), arraylength, infile);
  key += fread(Bz_pert, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*9)
    ghl_error("An error has occured with reading in perturbed data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  const double poison = 0.0/0.0;

  for(int i=0;i<arraylength;i++) {

    // Define the various GRHayL structs for the unit tests
    con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);
    metric_quantities ADM_metric;
    primitive_quantities prims;
    conservative_quantities cons, cons_undens;

    // Read initial data accompanying trusted output
    ghl_initialize_metric(lapse[i],
                      betax[i], betay[i], betaz[i],
                      gxx[i], gxy[i], gxz[i],
                      gyy[i], gyz[i], gzz[i],
                      &ADM_metric);

    ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

    ghl_initialize_primitives(
                      rho_b[i], press[i], eps[i],
                      vx[i], vy[i], vz[i],
                      Bx[i], By[i], Bz[i],
                      poison, poison, poison,
                      &prims);

    ghl_initialize_conservatives(rho_star[i], tau[i],
                             S_x[i], S_y[i], S_z[i],
                             poison, poison, &cons);

    //This uses the Font fix method to compute primitives from conservatives.
    ghl_undensitize_conservatives(metric_aux.psi6, &cons, &cons_undens);
    ghl_guess_primitives(&eos, &ADM_metric, &metric_aux, &cons, &prims);
    const int check = ghl_hybrid_Font_fix(&params, &eos, &ADM_metric, &metric_aux, &cons_undens, &prims, &diagnostics);
    if( check != c2p_check[i] )
      ghl_error("ghl_hybrid_Font1D has different return value: %d vs %d\n", check, c2p_check[i]);

    primitive_quantities prims_trusted, prims_pert;
    ghl_initialize_primitives(
                      rho_b_trusted[i], press_trusted[i], eps_trusted[i],
                      vx_trusted[i], vy_trusted[i], vz_trusted[i],
                      Bx_trusted[i], By_trusted[i], Bz_trusted[i],
                      poison, poison, poison,
                      &prims_trusted);

    ghl_initialize_primitives(
                      rho_b_pert[i], press_pert[i], eps_pert[i],
                      vx_pert[i], vy_pert[i], vz_pert[i],
                      Bx_pert[i], By_pert[i], Bz_pert[i],
                      poison, poison, poison,
                      &prims_pert);

    validate_primitives(params.evolve_entropy, &eos, &prims_trusted, &prims, &prims_pert);
  }
  ghl_info("ghl_hybrid_Font1D function test has passed!\n");
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(rho_b); free(press); free(eps);
  free(vx); free(vy); free(vz);
  free(Bx); free(By); free(Bz);
  free(rho_star); free(tau);
  free(S_x); free(S_y); free(S_z);
  free(rho_b_trusted); free(press_trusted); free(eps_trusted);
  free(vx_trusted); free(vy_trusted); free(vz_trusted);
  free(Bx_trusted); free(By_trusted); free(Bz_trusted);
  free(rho_b_pert); free(press_pert); free(eps_pert);
  free(vx_pert); free(vy_pert); free(vz_pert);
  free(Bx_pert); free(By_pert); free(Bz_pert);
}
