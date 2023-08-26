// Thorn      : GRHayL
// File       : unit_test_data_con2prim.c
// Author(s)  : Leo Werneck & Samuel Cupp
// Description: In this file we provide an extensive unit test of
//              the Con2Prim gem.
#include "unit_tests.h"


int main(int argc, char **argv) {

  FILE* input = fopen_with_check("ET_Legacy_conservs_input.bin", "rb");

  int npoints;
  int key = fread(&npoints, sizeof(int), 1, input);
  if( key != 1 || npoints < 1 )
    ghl_error("An error has occured with reading the grid size. "
                 "Please check that Noble2D_initial_data.bin"
                 "is up-to-date with current test version.\n");
  const int arraylength = npoints*npoints;

  const double poison = 1e300;
  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {None,None,None};
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100; //Taken from magnetizedTOV.par

  const int neos = 1;
  const double W_max = 10.0; //IGM default
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300; //IGM default
  const double Gamma_th = 2.0; //Taken from magnetizedTOV.par
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.
  ghl_parameters params;
  ghl_initialize_params(Noble2D, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess,
                    Psi6threshold, 0 /*Cupp Fix*/, 0 /*Lorenz damping factor*/, &params);

  ghl_eos_parameters eos;
  ghl_initialize_hybrid_eos_functions_and_params(W_max,
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
    ghl_error("An error has occured with reading initial data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(input);

  // The output for this test is provided by IllinoisGRMHD via the ET
  // for validation with legacy code
  FILE* output = fopen_with_check("ET_Legacy_conservs_output.bin", "rb");

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
    ghl_error("An error has occured with reading trusted data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(output);

  output = fopen_with_check("ET_Legacy_conservs_output_pert.bin", "rb");

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
    ghl_error("An error has occured with reading perturbed data. "
                 "Please check that comparison data "
                 "is up-to-date with current test version.\n");

  fclose(output);

  //Parallelizing this also needs parallel sum of abs/rel error arrays
  for(int index=0; index<arraylength; index++) {
    // Define the various GRHayL structs for the unit tests
    ghl_con2prim_diagnostics diagnostics;
    ghl_initialize_diagnostics(&diagnostics);
    ghl_metric_quantities ADM_metric;
    ghl_primitive_quantities prims;
    ghl_conservative_quantities cons;
    ghl_stress_energy Tmunu;

    ghl_initialize_metric(lapse[index],
                      betax[index], betay[index], betaz[index],
                      gxx[index], gxy[index], gxz[index],
                      gyy[index], gyz[index], gzz[index],
                      &ADM_metric);

    ghl_ADM_aux_quantities metric_aux;
    ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

    ghl_initialize_primitives(rho_b[index], press[index], eps[index],
                          vx[index], vy[index], vz[index],
                          Bx[index], By[index], Bz[index],
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims);

    ghl_initialize_conservatives(rho_star[index], tau[index],
                             S_x[index], S_y[index], S_z[index],
                             poison, poison, &cons);

    // Enforce limits on primitive variables and recompute conservatives.
    diagnostics.speed_limited = ghl_enforce_primitive_limits_and_compute_u0(&params, &eos, &ADM_metric, &prims);
    ghl_compute_conservs_and_Tmunu(&ADM_metric, &metric_aux, &prims, &cons, &Tmunu);

    // Here, we call the ghl_return_* functions and then repack the struct from that data. These functions are too
    // small to need an individual test, and they are primarily used by Con2Prim anyways.
    double rho_tmp, press_tmp, eps_tmp, vx_tmp, vy_tmp, vz_tmp, Bx_tmp, By_tmp, Bz_tmp, ent_tmp, Ye_tmp, temp_tmp;
    double rhos_tmp, tau_tmp, Sx_tmp, Sy_tmp, Sz_tmp, ent_s_tmp, Ye_s_tmp;
    double Ttt_tmp, Ttx_tmp, Tty_tmp, Ttz_tmp, Txx_tmp, Txy_tmp, Txz_tmp, Tyy_tmp, Tyz_tmp, Tzz_tmp;

    ghl_return_primitives(&prims,
                      &rho_tmp, &press_tmp, &eps_tmp,
                      &vx_tmp, &vy_tmp, &vz_tmp,
                      &Bx_tmp, &By_tmp, &Bz_tmp,
                      &ent_tmp, &Ye_tmp, &temp_tmp);

    ghl_return_conservatives(&cons,
                         &rhos_tmp, &tau_tmp,
                         &Sx_tmp, &Sy_tmp, &Sz_tmp,
                         &ent_s_tmp, &Ye_s_tmp);

    ghl_return_stress_energy(&Tmunu,
                         &Ttt_tmp, &Ttx_tmp,
                         &Tty_tmp, &Ttz_tmp,
                         &Txx_tmp, &Txy_tmp,
                         &Txz_tmp, &Tyy_tmp,
                         &Tyz_tmp, &Tzz_tmp);

    ghl_initialize_primitives(rho_tmp, press_tmp, eps_tmp,
                          vx_tmp, vy_tmp, vz_tmp,
                          Bx_tmp, By_tmp, Bz_tmp,
                          ent_tmp, Ye_tmp, temp_tmp,
                          &prims);

    ghl_initialize_conservatives(rhos_tmp, tau_tmp,
                             Sx_tmp, Sy_tmp, Sz_tmp,
                             ent_s_tmp, Ye_s_tmp, &cons);

    ghl_initialize_stress_energy(Ttt_tmp, Ttx_tmp,
                             Tty_tmp, Ttz_tmp,
                             Txx_tmp, Txy_tmp,
                             Txz_tmp, Tyy_tmp,
                             Tyz_tmp, Tzz_tmp,
                             &Tmunu);

    // Now, we load the trusted/perturbed data for this index and ghl_pert_test_fail the computed results.
    ghl_primitive_quantities prims_trusted, prims_pert;
    ghl_conservative_quantities cons_trusted, cons_pert;
    ghl_stress_energy Tmunu_trusted, Tmunu_pert;

    ghl_initialize_primitives(rho_b_trusted[index], press_trusted[index], prims.eps, // Old code has no eps variable
                          vx_trusted[index], vy_trusted[index], vz_trusted[index],
                          poison, poison, poison,
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims_trusted);

    ghl_initialize_primitives(rho_b_pert[index], press_pert[index], prims.eps, // Old code has no eps variable
                          vx_pert[index], vy_pert[index], vz_pert[index],
                          poison, poison, poison,
                          poison, poison, poison, // entropy, Y_e, temp
                          &prims_pert);

    ghl_initialize_conservatives(rho_star_trusted[index], tau_trusted[index],
                             S_x_trusted[index], S_y_trusted[index], S_z_trusted[index],
                             poison, poison, &cons_trusted);

    ghl_initialize_conservatives(rho_star_pert[index], tau_pert[index],
                             S_x_pert[index], S_y_pert[index], S_z_pert[index],
                             poison, poison, &cons_pert);

    ghl_initialize_stress_energy(Ttt_trusted[index], Ttx_trusted[index],
                             Tty_trusted[index], Ttz_trusted[index],
                             Txx_trusted[index], Txy_trusted[index],
                             Txz_trusted[index], Tyy_trusted[index],
                             Tyz_trusted[index], Tzz_trusted[index],
                             &Tmunu_trusted);

    ghl_initialize_stress_energy(Ttt_pert[index], Ttx_pert[index],
                             Tty_pert[index], Ttz_pert[index],
                             Txx_pert[index], Txy_pert[index],
                             Txz_pert[index], Tyy_pert[index],
                             Tyz_pert[index], Tzz_pert[index],
                             &Tmunu_pert);

    ghl_pert_test_fail_primitives(params.evolve_entropy, &eos, &prims_trusted, &prims, &prims_pert);
    ghl_pert_test_fail_conservatives(params.evolve_entropy, &cons_trusted, &cons, &cons_pert);
    ghl_pert_test_fail_stress_energy(&Tmunu_trusted, &Tmunu, &Tmunu_pert);
  }
  ghl_info("ET_Legacy primitives-to-conservatives test has passed!\n");
  free(lapse);
  free(betax); free(betay); free(betaz);
  free(gxx); free(gxy); free(gxz);
  free(gyy); free(gyz); free(gzz);
  free(rho_b); free(press); free(eps);
  free(vx); free(vy); free(vz);
  free(Bx); free(By); free(Bz);
  free(rho_star); free(tau);
  free(S_x); free(S_y); free(S_z);
  free(rho_b_trusted); free(press_trusted);
  free(vx_trusted); free(vy_trusted); free(vz_trusted);
  free(rho_star_trusted); free(tau_trusted);
  free(S_x_trusted); free(S_y_trusted); free(S_z_trusted);
  free(Ttt_trusted); free(Ttx_trusted); free(Tty_trusted);
  free(Ttz_trusted); free(Txx_trusted); free(Txy_trusted);
  free(Txz_trusted); free(Tyy_trusted); free(Tyz_trusted); free(Tzz_trusted);
  free(rho_b_pert); free(press_pert);
  free(vx_pert); free(vy_pert); free(vz_pert);
  free(rho_star_pert); free(tau_pert);
  free(S_x_pert); free(S_y_pert); free(S_z_pert);
  free(Ttt_pert); free(Ttx_pert); free(Tty_pert);
  free(Ttz_pert); free(Txx_pert); free(Txy_pert);
  free(Txz_pert); free(Tyy_pert); free(Tyz_pert); free(Tzz_pert);
}
