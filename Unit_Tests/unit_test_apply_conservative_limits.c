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
  ghl_initialize_params(None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
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

  infile = fopen_with_check("apply_conservative_limits_input.bin","rb");
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

  // Allocate memory for the trusted conservative data
  double *rho_star_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *tau_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_trusted = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("apply_conservative_limits_output.bin","rb");
  key  = fread(rho_star_trusted, sizeof(double), arraylength, infile);
  key += fread(tau_trusted, sizeof(double), arraylength, infile);
  key += fread(S_x_trusted, sizeof(double), arraylength, infile);
  key += fread(S_y_trusted, sizeof(double), arraylength, infile);
  key += fread(S_z_trusted, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*5)
    ghl_error("An error has occured with reading in trusted data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  // Allocate memory for the perturbed conservative data
  double *rho_star_pert = (double*) malloc(sizeof(double)*arraylength);
  double *tau_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_pert = (double*) malloc(sizeof(double)*arraylength);

  infile = fopen_with_check("apply_conservative_limits_output_pert.bin","rb");
  key  = fread(rho_star_pert, sizeof(double), arraylength, infile);
  key += fread(tau_pert, sizeof(double), arraylength, infile);
  key += fread(S_x_pert, sizeof(double), arraylength, infile);
  key += fread(S_y_pert, sizeof(double), arraylength, infile);
  key += fread(S_z_pert, sizeof(double), arraylength, infile);

  fclose(infile);
  if(key != arraylength*5)
    ghl_error("An error has occured with reading in perturbed data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  const double poison = 0.0/0.0;

  for(int i=0;i<arraylength;i++) {
    // Define the various GRHayL structs for the unit tests
    con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);
    metric_quantities ADM_metric;
    primitive_quantities prims;
    conservative_quantities cons;

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

    //This applies inequality fixes on the conservatives
    if(i == arraylength-1 || i == arraylength-2)
      params.psi6threshold = 0.0;
    ghl_apply_conservative_limits(&params, &eos, &ADM_metric, &prims, &cons, &diagnostics);
    if(i == arraylength-1 || i == arraylength-2)
      params.psi6threshold = Psi6threshold;

    conservative_quantities cons_trusted, cons_pert;
    ghl_initialize_conservatives(rho_star_trusted[i], tau_trusted[i],
                             S_x_trusted[i], S_y_trusted[i], S_z_trusted[i],
                             poison, poison, &cons_trusted);

    ghl_initialize_conservatives(rho_star_pert[i], tau_pert[i],
                             S_x_pert[i], S_y_pert[i], S_z_pert[i],
                             poison, poison, &cons_pert);


    ghl_validate_conservatives(params.evolve_entropy, &cons_trusted, &cons, &cons_pert);
  }
  ghl_info("ghl_apply_conservative_limits function test has passed!\n");
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(rho_b); free(press); free(eps);
  free(vx); free(vy); free(vz);
  free(Bx); free(By); free(Bz);
  free(rho_star); free(tau);
  free(S_x); free(S_y); free(S_z);
  free(rho_star_trusted); free(tau_trusted);
  free(S_x_trusted); free(S_y_trusted); free(S_z_trusted);
  free(rho_star_pert); free(tau_pert);
  free(S_x_pert); free(S_y_pert); free(S_z_pert);
}
