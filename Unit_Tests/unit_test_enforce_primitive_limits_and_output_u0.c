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
  const double gamma_th = 2.0; //Taken from magnetizedTOV.par
  const double rho_ppoly[1] = {0.0};
  const double gamma_ppoly[1] = {2.0};
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
             gamma_ppoly, k_ppoly0, gamma_th,
             &eos);

  con2prim_diagnostics diagnostics;
  initialize_diagnostics(&diagnostics);

  char filename[100];

  // Initialize the needed data files
  sprintf(filename,"C2P_Noble2D_initial_data.bin");
  FILE* initial_data = fopen(filename,"rb");
  check_file_was_successfully_open(initial_data, filename);

  sprintf(filename,"C2P_enforce_primitive_limits_and_output_u0.bin");
  FILE* infile = fopen(filename,"rb");
  check_file_was_successfully_open(infile, filename);

  sprintf(filename,"C2P_enforce_primitive_limits_and_output_u0_pert.bin");
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

  // Allocate memory for the initial primitive data, trusted output,
  // and perturbed output
  double *rho_b = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vy = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);
  double *eps = (double*) malloc(sizeof(double)*arraylength);
  double *Bx = (double*) malloc(sizeof(double)*arraylength);
  double *By = (double*) malloc(sizeof(double)*arraylength);
  double *Bz = (double*) malloc(sizeof(double)*arraylength);

  double *rho_b_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *press_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vx_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vy_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *vz_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *eps_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *By_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_trusted = (double*) malloc(sizeof(double)*arraylength);

  double *rho_b_pert = (double*) malloc(sizeof(double)*arraylength);
  double *press_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vy_pert = (double*) malloc(sizeof(double)*arraylength);
  double *vz_pert = (double*) malloc(sizeof(double)*arraylength);
  double *eps_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Bx_pert = (double*) malloc(sizeof(double)*arraylength);
  double *By_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Bz_pert = (double*) malloc(sizeof(double)*arraylength);

  // These arrays may not be used, but it's simpler to just declare them
  // either way.
  double *entropy = (double*) malloc(sizeof(double)*arraylength);
  double *entropy_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *entropy_pert = (double*) malloc(sizeof(double)*arraylength);
  double *Y_e = (double*) malloc(sizeof(double)*arraylength);
  double *Y_e_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *Y_e_pert = (double*) malloc(sizeof(double)*arraylength);
  double *temperature = (double*) malloc(sizeof(double)*arraylength);
  double *temperature_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *temperature_pert = (double*) malloc(sizeof(double)*arraylength);

  // This function also output u^0, so we need to store trusted output
  // and perturbed output for this as well
  double *u0_trusted = (double*) malloc(sizeof(double)*arraylength);
  double *u0_pert = (double*) malloc(sizeof(double)*arraylength);

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

      read_primitive_binary(eos.eos_type, params.evolve_entropy, &rho_b_trusted[index], &press_trusted[index],
                            &vx_trusted[index], &vy_trusted[index], &vz_trusted[index], &eps_trusted[index],
                            &Bx_trusted[index], &By_trusted[index], &Bz_trusted[index],
                            &entropy_trusted[index], &Y_e_trusted[index], &temperature_trusted[index],
                            infile);

      int key = fread(&u0_trusted[index], sizeof(double), 1, infile);
      if( key != 1)
        grhayl_error("An error has occured with reading in trusted primitive data."
                     "Please check that comparison data "
                     "is up-to-date with current test version.\n");

      read_primitive_binary(eos.eos_type, params.evolve_entropy, &rho_b_pert[index], &press_pert[index],
                            &vx_pert[index], &vy_pert[index], &vz_pert[index], &eps_pert[index],
                            &Bx_pert[index], &By_pert[index], &Bz_pert[index],
                            &entropy_pert[index], &Y_e_pert[index], &temperature_pert[index],
                            inpert);

      key = fread(&u0_pert[index], sizeof(double), 1, inpert);
      if( key != 1)
        grhayl_error("An error has occured with reading in trusted primitive data."
                     "Please check that comparison data "
                     "is up-to-date with current test version.\n");
  }
  fclose(initial_data);
  fclose(infile);
  fclose(inpert);

  
#pragma omp parallel for
  for(int i=0;i<arraylength;i++) {

    // Define the various GRHayL structs for the unit tests
    metric_quantities metric;
    primitive_quantities prims;

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

    //This applies the inequality (or "Faber") fixes on the conservatives
    double u0;
    enforce_primitive_limits_and_output_u0(&params, &eos, &metric, &prims, &u0, &diagnostics);

    primitive_quantities prims_trusted, prims_pert;
    initialize_primitives(
                      rho_b_trusted[i], press_trusted[i], eps_trusted[i],
                      vx_trusted[i], vy_trusted[i], vz_trusted[i],
                      Bx_trusted[i], By_trusted[i], Bz_trusted[i],
                      entropy_trusted[i], Y_e_trusted[i], temperature_trusted[i],
                      &prims_trusted);

    initialize_primitives(
                      rho_b_pert[i], press_pert[i], eps_pert[i],
                      vx_pert[i], vy_pert[i], vz_pert[i],
                      Bx_pert[i], By_pert[i], Bz_pert[i],
                      entropy_pert[i], Y_e_pert[i], temperature_pert[i],
                      &prims_pert);

    validate_primitives(eos.eos_type, params.evolve_entropy, &prims, &prims_trusted, &prims_pert);
    // The value of u0 varies by more than roundoff because roundoff alone can trigger the velocity
    // limiting if it is on the edge, leading to a higher difference in some edge cases. Thus,
    // the default roundoff tolerance of validate() is not a valid floor for the possible
    // relative difference.
    if( validate(u0_trusted[i], u0, u0_pert[i]) )
      grhayl_error("Test has failed! The computed u0 does not fall within tolerance.\n"
                   "   u0_trusted %.15e\n"
                   "   u0_compute %.15e\n"
                   "   u0_perturb %.15e\n"
                   "   rel_diff %.15e %.15e\n", u0_trusted[i], u0, u0_pert[i], relative_error(u0_trusted[i], u0), relative_error(u0_trusted[i], u0_pert[i]));
  }
  printf("Completed test for routine enforce_primitive_limits_and_output_u0\n");
  return 0;
}
