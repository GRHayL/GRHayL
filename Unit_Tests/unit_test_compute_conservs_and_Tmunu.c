#include "unit_tests.h"

int main(int argc, char **argv) {

  FILE* initial_data = fopen("Noble2D_initial_data.bin","rb");
  check_file_was_successfully_open(initial_data, "Noble2D_initial_data.bin");

  //Number of sampling points in density and pressure/temperature
  int npoints;
  const int key = fread(&npoints, sizeof(int), 1, initial_data);
  if( key != 1)
    grhayl_error("An error has occured with reading the grid size. "
                 "Please check that Noble2D_initial_data.bin"
                 "is up-to-date with current test version.\n");
  const int arraylength = npoints*npoints;

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {None,None,None};
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100; //Taken from magnetizedTOV.par
  const int update_Tmunu = 1; //IGM default

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
  GRHayL_parameters params;
  initialize_GRHayL(None, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess,
                    Psi6threshold, update_Tmunu, 1 /*Cupp Fix*/, 0 /*Lorenz damping factor*/, &params);

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

  // Initialize the needed data files
  FILE* infile = fopen("compute_conservs_and_Tmunu.bin","rb");
  check_file_was_successfully_open(infile, "compute_conservs_and_Tmunu.bin");

  FILE* inpert = fopen("compute_conservs_and_Tmunu_pert.bin","rb");
  check_file_was_successfully_open(inpert, "compute_conservs_and_Tmunu_pert.bin");

  // Allocate memory for the metric data
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

  // Allocate memory for the initial primitive data
  double *rho_b = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vy = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);
  double *eps = (double*) malloc(sizeof(double)*arraylength);
  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  // These arrays may not be used, but it's simpler to just declare them
  // either way.
  double *entropy = (double*) malloc(sizeof(double)*arraylength);
  double *Y_e = (double*) malloc(sizeof(double)*arraylength);
  double *temperature = (double*) malloc(sizeof(double)*arraylength);

  // This function uses u^0, so we need to store it as well
  double *u0 = (double*) malloc(sizeof(double)*arraylength);

  // Allocate memory for trusted and perturbed output for the conservatives
  double *rho_star_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *tau_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_trusted = (double*) malloc(sizeof(double)*arraylength);

  double *rho_star_pert = (double*) malloc(sizeof(double)*arraylength);
  double *tau_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_x_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_y_pert = (double*) malloc(sizeof(double)*arraylength);
  double *S_z_pert = (double*) malloc(sizeof(double)*arraylength);

  // These arrays may not be used, but it's simpler to just declare them
  // either way.
  double *ent_cons_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Y_e_cons_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *ent_cons_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Y_e_cons_pert = (double*) malloc(sizeof(double)*arraylength);

  // Storage for the trusted and perturbed output for the stress-energy tensor
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

  // Can't parallelize because it could change the behavior of reading from the files
  for(int j=0;j<npoints;j++)
    for(int i=0;i<npoints;i++) {
      const int index = i + j*npoints;
      read_metric_binary(&lapse[index], &gxx[index], &gxy[index], &gxz[index],
                         &gyy[index], &gyz[index], &gzz[index], &betax[index],
                         &betay[index], &betaz[index], initial_data);

      read_primitive_binary(eos.eos_type, params.evolve_entropy, &rho_b[index], &press[index],
                            &vx[index], &vy[index], &vz[index], &eps[index],
                            &Bx[index], &By[index], &Bz[index],
                            &entropy[index], &Y_e[index], &temperature[index],
                            infile);

      int key = fread(&u0[index], sizeof(double), 1, infile);
      if( key != 1)
        grhayl_error("An error has occured with reading in trusted primitive data."
                     "Please check that comparison data "
                     "is up-to-date with current test version.\n");

      read_conservative_binary(params.evolve_entropy, &rho_star_trusted[index], &tau_trusted[index],
                               &S_x_trusted[index], &S_y_trusted[index], &S_z_trusted[index], &ent_cons_trusted[index],
                               infile);

      read_stress_energy_binary(&Ttt_trusted[index], &Ttx_trusted[index], &Tty_trusted[index], &Ttz_trusted[index],
                                &Txx_trusted[index], &Txy_trusted[index], &Txz_trusted[index], &Tyy_trusted[index],
                                &Tyz_trusted[index], &Tzz_trusted[index], infile);

      read_conservative_binary(params.evolve_entropy, &rho_star_pert[index], &tau_pert[index],
                               &S_x_pert[index], &S_y_pert[index], &S_z_pert[index], &ent_cons_pert[index],
                               inpert);

      read_stress_energy_binary(&Ttt_pert[index], &Ttx_pert[index], &Tty_pert[index], &Ttz_pert[index],
                                &Txx_pert[index], &Txy_pert[index], &Txz_pert[index], &Tyy_pert[index],
                                &Tyz_pert[index], &Tzz_pert[index], inpert);
  }
  fclose(initial_data);
  fclose(infile);
  fclose(inpert);


  for(int i=0;i<arraylength;i++) {

    // Define the various GRHayL structs for the unit tests
    con2prim_diagnostics diagnostics;
    initialize_diagnostics(&diagnostics);
    metric_quantities metric;
    primitive_quantities prims;
    conservative_quantities cons;
    stress_energy Tmunu;

    // Read initial data accompanying trusted output
    initialize_metric(lapse[i], gxx[i], gxy[i], gxz[i],
                      gyy[i], gyz[i], gzz[i], betax[i],
                      betay[i], betaz[i], &metric);

    initialize_primitives(
                      rho_b[i], press[i], eps[i],
                      vx[i], vy[i], vz[i],
                      Bx[i], By[i], Bz[i],
                      entropy[i], Y_e[i], temperature[i],
                      &prims);
    prims.u0 = u0[i];

    //This computes the conservatives and stress-energy tensor from the new primitives
    compute_conservs_and_Tmunu(&params, &metric, &prims, &cons, &Tmunu);

    conservative_quantities cons_trusted, cons_pert;
    stress_energy Tmunu_trusted, Tmunu_pert;

    initialize_conservatives(rho_star_trusted[i], tau_trusted[i],
                             S_x_trusted[i], S_y_trusted[i], S_z_trusted[i],
                             Y_e_cons_trusted[i], ent_cons_trusted[i], &cons_trusted);

    initialize_conservatives(rho_star_pert[i], tau_pert[i],
                             S_x_pert[i], S_y_pert[i], S_z_pert[i],
                             Y_e_cons_pert[i], ent_cons_pert[i], &cons_pert);

    initialize_stress_energy(Ttt_trusted[i], Ttx_trusted[i],
                             Tty_trusted[i], Ttz_trusted[i],
                             Txx_trusted[i], Txy_trusted[i],
                             Txz_trusted[i], Tyy_trusted[i],
                             Tyz_trusted[i], Tzz_trusted[i],
                             &Tmunu_trusted);

    initialize_stress_energy(Ttt_pert[i], Ttx_pert[i],
                             Tty_pert[i], Ttz_pert[i],
                             Txx_pert[i], Txy_pert[i],
                             Txz_pert[i], Tyy_pert[i],
                             Tyz_pert[i], Tzz_pert[i],
                             &Tmunu_pert);

    validate_conservatives(params.evolve_entropy, &cons_trusted, &cons, &cons_pert);
    validate_stress_energy(&Tmunu_trusted, &Tmunu, &Tmunu_pert);
  }
  grhayl_info("compute_conservs_and_Tmunu function test has passed!\n");
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

  free(Y_e);
  free(entropy);
  free(temperature);
  free(Y_e_cons_trusted); free(ent_cons_trusted);
  free(Y_e_cons_pert); free(ent_cons_pert);
}
