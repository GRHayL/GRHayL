#include "unit_tests.h"

int main(int argc, char **argv) {

  //Number of sampling points in density and pressure/temperature
  const int npoints = 80;
  const int arraylength = npoints*npoints;

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {None,None,None};
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100; //Taken from magnetizedTOV.par
  const int update_Tmunu = 1; //IGM default

  const int eos_type = 0; // Hybrid=0, Tabulated=1;
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
  initialize_GRHayL(None, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess, Psi6threshold, update_Tmunu, 1 /*Cupp Fix*/, &params);

  eos_parameters eos;
  initialize_general_eos(eos_type, W_max,
             rho_b_min, rho_b_min, rho_b_max,
             &eos);

  initialize_hybrid_functions(&eos);

  initialize_hybrid_eos(neos, rho_ppoly,
             Gamma_ppoly, k_ppoly0, Gamma_th,
             &eos);

  char filename[100];

  // Initialize the needed data files
  sprintf(filename,"C2P_Noble2D_initial_data.bin");
  FILE* initial_data = fopen(filename,"rb");
  check_file_was_successfully_open(initial_data, filename);

  sprintf(filename,"C2P_compute_conservs_and_Tmunu.bin");
  FILE* infile = fopen(filename,"rb");
  check_file_was_successfully_open(infile, filename);

  sprintf(filename,"C2P_compute_conservs_and_Tmunu_pert.bin");
  FILE* inpert = fopen(filename,"rb");
  check_file_was_successfully_open(inpert, filename);

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

  
#pragma omp parallel for
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

    //This computes the conservatives and stress-energy tensor from the new primitives
    compute_conservs_and_Tmunu(&params, &eos, &metric, &prims, u0[i], &cons, &Tmunu);

    conservative_quantities cons_trusted, cons_pert;
    stress_energy Tmunu_trusted, Tmunu_pert;

    initialize_conservatives(rho_star_trusted[i], tau_trusted[i],
                             S_x_trusted[i], S_y_trusted[i], S_z_trusted[i],
                             Y_e_cons_trusted[i], ent_cons_trusted[i], &cons_trusted);

    initialize_conservatives(rho_star_pert[i], tau_pert[i],
                             S_x_pert[i], S_y_pert[i], S_z_pert[i],
                             Y_e_cons_pert[i], ent_cons_pert[i], &cons_pert);

    Tmunu_trusted.Ttt = Ttt_trusted[i];
    Tmunu_trusted.Ttx = Ttx_trusted[i];
    Tmunu_trusted.Tty = Tty_trusted[i];
    Tmunu_trusted.Ttz = Ttz_trusted[i];
    Tmunu_trusted.Txx = Txx_trusted[i];
    Tmunu_trusted.Txy = Txy_trusted[i];
    Tmunu_trusted.Txz = Txz_trusted[i];
    Tmunu_trusted.Tyy = Tyy_trusted[i];
    Tmunu_trusted.Tyz = Tyz_trusted[i];
    Tmunu_trusted.Tzz = Tzz_trusted[i];

    Tmunu_pert.Ttt = Ttt_pert[i];
    Tmunu_pert.Ttx = Ttx_pert[i];
    Tmunu_pert.Tty = Tty_pert[i];
    Tmunu_pert.Ttz = Ttz_pert[i];
    Tmunu_pert.Txx = Txx_pert[i];
    Tmunu_pert.Txy = Txy_pert[i];
    Tmunu_pert.Txz = Txz_pert[i];
    Tmunu_pert.Tyy = Tyy_pert[i];
    Tmunu_pert.Tyz = Tyz_pert[i];
    Tmunu_pert.Tzz = Tzz_pert[i];

    validate_conservatives(params.evolve_entropy, &cons, &cons_trusted, &cons_pert);
    validate_stress_energy(&Tmunu, &Tmunu_trusted, &Tmunu_pert);
  }
  printf("Completed test for routine compute_conservs_and_Tmunu\n");
  return 0;
}
