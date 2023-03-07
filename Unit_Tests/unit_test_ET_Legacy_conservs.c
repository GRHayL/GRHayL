// Thorn      : GRHayL
// File       : unit_test_data_con2prim.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.
#include "unit_tests.h"


int main(int argc, char **argv) {

  FILE* input = fopen("ET_Legacy_conservs_input.bin", "rb");
  check_file_was_successfully_open(input, "ET_Legacy_conservs_input.bin");

  int npoints;
  int key = fread(&npoints, sizeof(int), 1, input);
  const int arraylength = npoints*npoints;

  double poison = 1e200;
  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  int backup_routine[3] = {None,None,None};
  bool calc_prims_guess = true;
  double Psi6threshold = 1e100; //Taken from magnetizedTOV.par
  int update_Tmunu = 1; //IGM default

  int neos = 1;
  double W_max = 10.0; //IGM default
  double rho_b_min = 1e-12;
  double rho_b_max = 1e300; //IGM default
  double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  double rho_ppoly[1] = {0.0};
  double Gamma_ppoly[1] = {2.0};
  double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  GRHayL_parameters params;
  initialize_GRHayL(Noble2D, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess, Psi6threshold, update_Tmunu, 1 /*Cupp Fix*/, &params);

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

  // Allocate memory for the metrics
  double *gxx = (double*) malloc(sizeof(double)*arraylength);
  double *gxy = (double*) malloc(sizeof(double)*arraylength);
  double *gxz = (double*) malloc(sizeof(double)*arraylength);
  double *gyy = (double*) malloc(sizeof(double)*arraylength);
  double *gyz = (double*) malloc(sizeof(double)*arraylength);
  double *gzz = (double*) malloc(sizeof(double)*arraylength);

  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for the primitives
  double *rho_b = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *eps = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vy = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);
  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for the conservatives 
  double *rho_star = (double*) malloc(sizeof(double)*arraylength);
  double *tau = (double*) malloc(sizeof(double)*arraylength);
  double *S_x = (double*) malloc(sizeof(double)*arraylength);
  double *S_y = (double*) malloc(sizeof(double)*arraylength);
  double *S_z = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for for trusted output
  double *rho_b_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *press_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vx_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vy_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vz_trusted = (double*) malloc(sizeof(double)*arraylength);

  double *rho_star_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *tau_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_trusted = (double*) malloc(sizeof(double)*arraylength);

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

  // Allocate memory for for perturbed output
  double *rho_b_pert = (double*) malloc(sizeof(double)*arraylength);
  double *press_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vy_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vz_pert = (double*) malloc(sizeof(double)*arraylength);

  double *rho_star_pert = (double*) malloc(sizeof(double)*arraylength);
  double *tau_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_pert = (double*) malloc(sizeof(double)*arraylength);

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

  // Read in data from file to ensure portability
  key  = fread(gxx,   sizeof(double), arraylength, input);
  key += fread(gxy,   sizeof(double), arraylength, input);
  key += fread(gxz,   sizeof(double), arraylength, input);
  key += fread(gyy,   sizeof(double), arraylength, input);
  key += fread(gyz,   sizeof(double), arraylength, input);
  key += fread(gzz,   sizeof(double), arraylength, input);
  key += fread(lapse, sizeof(double), arraylength, input);
  key += fread(betax, sizeof(double), arraylength, input);
  key += fread(betay, sizeof(double), arraylength, input);
  key += fread(betaz, sizeof(double), arraylength, input);

  key += fread(rho_b, sizeof(double), arraylength, input);
  key += fread(press, sizeof(double), arraylength, input);
  key += fread(eps, sizeof(double), arraylength, input);
  key += fread(vx, sizeof(double), arraylength, input);
  key += fread(vy, sizeof(double), arraylength, input);
  key += fread(vz, sizeof(double), arraylength, input);
  key += fread(Bx, sizeof(double), arraylength, input);
  key += fread(By, sizeof(double), arraylength, input);
  key += fread(Bz, sizeof(double), arraylength, input);

  key += fread(rho_star, sizeof(double), arraylength, input);
  key += fread(tau, sizeof(double), arraylength, input);
  key += fread(S_x, sizeof(double), arraylength, input);
  key += fread(S_y, sizeof(double), arraylength, input);
  key += fread(S_z, sizeof(double), arraylength, input);

  if( key != (10+9+5)*arraylength)
    grhayl_error("An error has occured with reading initial data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(input);

  // The output for this test is provided by IllinoisGRMHD via the ET
  // for validation with legacy code
  FILE* output = fopen("ET_Legacy_conservs_output.bin", "rb");
  check_file_was_successfully_open(output, "ET_Legacy_conservs_output.bin");

  key  = fread(rho_b_trusted, sizeof(double), arraylength, output);
  key += fread(press_trusted, sizeof(double), arraylength, output);
  key += fread(vx_trusted, sizeof(double), arraylength, output);
  key += fread(vy_trusted, sizeof(double), arraylength, output);
  key += fread(vz_trusted, sizeof(double), arraylength, output);

  key += fread(rho_star_trusted, sizeof(double), arraylength, output);
  key += fread(tau_trusted, sizeof(double), arraylength, output);
  key += fread(S_x_trusted, sizeof(double), arraylength, output);
  key += fread(S_y_trusted, sizeof(double), arraylength, output);
  key += fread(S_z_trusted, sizeof(double), arraylength, output);

  key += fread(Ttt_trusted,   sizeof(double), arraylength, output);
  key += fread(Ttx_trusted,   sizeof(double), arraylength, output);
  key += fread(Tty_trusted,   sizeof(double), arraylength, output);
  key += fread(Ttz_trusted,   sizeof(double), arraylength, output);
  key += fread(Txx_trusted,   sizeof(double), arraylength, output);
  key += fread(Txy_trusted,   sizeof(double), arraylength, output);
  key += fread(Txz_trusted,   sizeof(double), arraylength, output);
  key += fread(Tyy_trusted,   sizeof(double), arraylength, output);
  key += fread(Tyz_trusted,   sizeof(double), arraylength, output);
  key += fread(Tzz_trusted,   sizeof(double), arraylength, output);

  if( key != 20*arraylength)
    grhayl_error("An error has occured with reading trusted data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(output);

  output = fopen("ET_Legacy_conservs_output_pert.bin", "rb");
  check_file_was_successfully_open(output, "ET_Legacy_conservs_output_pert.bin");

  key  = fread(rho_b_pert, sizeof(double), arraylength, output);
  key += fread(press_pert, sizeof(double), arraylength, output);
  key += fread(vx_pert, sizeof(double), arraylength, output);
  key += fread(vy_pert, sizeof(double), arraylength, output);
  key += fread(vz_pert, sizeof(double), arraylength, output);

  key += fread(rho_star_pert, sizeof(double), arraylength, output);
  key += fread(tau_pert, sizeof(double), arraylength, output);
  key += fread(S_x_pert, sizeof(double), arraylength, output);
  key += fread(S_y_pert, sizeof(double), arraylength, output);
  key += fread(S_z_pert, sizeof(double), arraylength, output);

  key += fread(Ttt_pert,   sizeof(double), arraylength, output);
  key += fread(Ttx_pert,   sizeof(double), arraylength, output);
  key += fread(Tty_pert,   sizeof(double), arraylength, output);
  key += fread(Ttz_pert,   sizeof(double), arraylength, output);
  key += fread(Txx_pert,   sizeof(double), arraylength, output);
  key += fread(Txy_pert,   sizeof(double), arraylength, output);
  key += fread(Txz_pert,   sizeof(double), arraylength, output);
  key += fread(Tyy_pert,   sizeof(double), arraylength, output);
  key += fread(Tyz_pert,   sizeof(double), arraylength, output);
  key += fread(Tzz_pert,   sizeof(double), arraylength, output);

  if( key != 20*arraylength)
    grhayl_error("An error has occured with reading perturbed data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(output);

  //OMP
  for(int index=0; index<arraylength; index++) {
    // Define the various GRHayL structs for the unit tests
    con2prim_diagnostics diagnostics;
    initialize_diagnostics(&diagnostics);
    metric_quantities metric;
    primitive_quantities prims;
    conservative_quantities cons;
    stress_energy Tmunu;

    initialize_metric(lapse[index],
                      gxx[index], gxy[index], gxz[index],
                      gyy[index], gyz[index], gzz[index],
                      betax[index], betay[index], betaz[index],
                      &metric);

    initialize_primitives(rho_b[index], press[index], eps[index],
                          vx[index], vy[index], vz[index],
                          Bx[index], By[index], Bz[index],
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims);

    initialize_conservatives(rho_star[index], tau[index],
                             S_x[index], S_y[index], S_z[index],
                             poison, poison, &cons);

    // Enforce limits on primitive variables and recompute conservatives.
    enforce_primitive_limits_and_compute_u0(&params, &eos, &metric, &prims, &diagnostics);
    compute_conservs_and_Tmunu(&params, &eos, &metric, &prims, &cons, &Tmunu);

    // Now, we load the trusted/perturbed data for this index and validate the computed results.
    primitive_quantities prims_trusted, prims_pert;
    conservative_quantities cons_trusted, cons_pert;
    stress_energy Tmunu_trusted, Tmunu_pert;

    initialize_primitives(rho_b_trusted[index], press_trusted[index], prims.eps, // Old code has no eps variable
                          vx_trusted[index], vy_trusted[index], vz_trusted[index],
                          poison, poison, poison,
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims_trusted);

    initialize_primitives(rho_b_pert[index], press_pert[index], prims.eps, // Old code has no eps variable
                          vx_pert[index], vy_pert[index], vz_pert[index],
                          poison, poison, poison,
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims_pert);

    initialize_conservatives(rho_star_trusted[index], tau_trusted[index],
                             S_x_trusted[index], S_y_trusted[index], S_z_trusted[index],
                             poison, poison, &cons_trusted);

    initialize_conservatives(rho_star_pert[index], tau_pert[index],
                             S_x_pert[index], S_y_pert[index], S_z_pert[index],
                             poison, poison, &cons_pert);

    Tmunu_trusted.Ttt = Ttt_trusted[index];
    Tmunu_trusted.Ttx = Ttx_trusted[index];
    Tmunu_trusted.Tty = Tty_trusted[index];
    Tmunu_trusted.Ttz = Ttz_trusted[index];
    Tmunu_trusted.Txx = Txx_trusted[index];
    Tmunu_trusted.Txy = Txy_trusted[index];
    Tmunu_trusted.Txz = Txz_trusted[index];
    Tmunu_trusted.Tyy = Tyy_trusted[index];
    Tmunu_trusted.Tyz = Tyz_trusted[index];
    Tmunu_trusted.Tzz = Tzz_trusted[index];

    Tmunu_pert.Ttt = Ttt_pert[index];
    Tmunu_pert.Ttx = Ttx_pert[index];
    Tmunu_pert.Tty = Tty_pert[index];
    Tmunu_pert.Ttz = Ttz_pert[index];
    Tmunu_pert.Txx = Txx_pert[index];
    Tmunu_pert.Txy = Txy_pert[index];
    Tmunu_pert.Txz = Txz_pert[index];
    Tmunu_pert.Tyy = Tyy_pert[index];
    Tmunu_pert.Tyz = Tyz_pert[index];

    validate_primitives(params.evolve_entropy, &eos, &prims_trusted, &prims, &prims_pert);
    validate_conservatives(params.evolve_entropy, &cons_trusted, &cons, &cons_pert);
    validate_stress_energy(&Tmunu_trusted, &Tmunu, &Tmunu_pert);
  }
  return 0;
}
