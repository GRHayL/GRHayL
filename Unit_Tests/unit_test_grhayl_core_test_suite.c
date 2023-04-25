#include "unit_tests.h"

int main(int argc, char **argv) {

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {None,None,None};
  const bool calc_prims_guess = true;
  const double Psi6threshold = 1e100;
  const int update_Tmunu = 1;

  const int neos = 1;
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Here, we initialize the structs that are (usually) static during a simulation.
  GRHayL_parameters params;
  initialize_GRHayL(None, backup_routine, false /*evolve entropy*/, false /*evolve temperature*/, calc_prims_guess,
                    Psi6threshold, update_Tmunu, 1 /*Cupp Fix*/, 0.0 /*Lorenz damping factor*/, &params);

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

  // First test: reset_prims_to_atmosphere() function
  primitive_quantities prims;

  // poison initial primitive data so reset is obvious
  const double poison = 0.0/0.0;
  prims.rho = poison;
  prims.press = poison;
  prims.vx = poison;
  prims.vy = poison;
  prims.vz = poison;
  prims.eps = poison;
  prims.entropy = poison;
  prims.Y_e = poison;
  prims.temperature = poison;

  reset_prims_to_atmosphere(&eos, &prims);

  // Function just sets data to eos data
  if(prims.rho         != eos.rho_atm
  || prims.press       != eos.press_atm
  || prims.eps         != eos.eps_atm
  || prims.entropy     != eos.entropy_atm
  || prims.Y_e         != eos.Ye_atm
  || prims.temperature != eos.T_atm
  || prims.vx          != 0.0
  || prims.vy          != 0.0
  || prims.vz          != 0.0)
    grhayl_error("grhayl_core_test_suite has failed for reset_prims_to_atmosphere().\n"
                 "  Struct output: %e %e %e %e %e %e %e %e %e\n"
                 "  EOS atm data:  %e %e %e %e %e %e %e %e %e\n",
                 prims.rho, prims.press, prims.eps, prims.entropy, prims.Y_e,
                 prims.temperature, prims.vx, prims.vy, prims.vz,
                 eos.rho_atm, eos.press_atm, eos.eps_atm, eos.entropy_atm, eos.Ye_atm, eos.T_atm, 0.0, 0.0, 0.0);

  // Second test: GRHayL_enforce_detgtij_and_initialize_metric() function
  // Valid metrics should return near-identical (round-off level) metrics
  metric_quantities new_metric;

  // read in metric data
  FILE* infile = fopen("grhayl_core_test_suite_input.bin","rb");
  check_file_was_successfully_open(infile, "grhayl_core_test_suite_input.bin");

  int arraylength;
  int key = fread(&arraylength, sizeof(int), 1, infile);

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
    grhayl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  const double rel_tol = 1e-14;
  for(int i=0; i<arraylength; i++) {
    GRHayL_enforce_detgtij_and_initialize_metric(
        lapse[i],
        gxx[i], gxy[i], gxz[i],
        gyy[i], gyz[i], gzz[i],
        betax[i], betay[i], betaz[i],
        &new_metric);

  if(relative_error(lapse[i], new_metric.lapse) > rel_tol
  || relative_error(betax[i], new_metric.betax) > rel_tol
  || relative_error(betay[i], new_metric.betay) > rel_tol
  || relative_error(betaz[i], new_metric.betaz) > rel_tol
  || relative_error(gxx[i],   new_metric.adm_gxx) > rel_tol
  || relative_error(gxy[i],   new_metric.adm_gxy) > rel_tol
  || relative_error(gxz[i],   new_metric.adm_gxz) > rel_tol
  || relative_error(gyy[i],   new_metric.adm_gyy) > rel_tol
  || relative_error(gyz[i],   new_metric.adm_gyz) > rel_tol
  || relative_error(gzz[i],   new_metric.adm_gzz) > rel_tol)
    grhayl_error("grhayl_core_test_suite has failed for GRHayL_enforce_detgtij_and_initialize_metric().\n"
                 "  input metric:  %e %e %e %e %e %e %e %e %e %e\n"
                 "  output metric: %e %e %e %e %e %e %e %e %e %e\n",
                 lapse[i], betax[i], betay[i], betaz[i],
                 gxx[i], gxy[i], gxz[i], gyy[i], gyz[i], gzz[i],
                 new_metric.lapse, new_metric.betax, new_metric.betay, new_metric.betaz,
                 new_metric.adm_gxx, new_metric.adm_gxy, new_metric.adm_gxz,
                 new_metric.adm_gyy, new_metric.adm_gyz, new_metric.adm_gzz);

  }

  // For a final test, we trigger the -detg warning
  gxy[0] += 10.0;

    GRHayL_enforce_detgtij_and_initialize_metric(
        lapse[0],
        gxx[0], gxy[0], gxz[0],
        gyy[0], gyz[0], gzz[0],
        betax[0], betay[0], betaz[0],
        &new_metric);
 
  if(relative_error(lapse[0], new_metric.lapse) > rel_tol
  || relative_error(betax[0], new_metric.betax) > rel_tol
  || relative_error(betay[0], new_metric.betay) > rel_tol
  || relative_error(betaz[0], new_metric.betaz) > rel_tol
  || relative_error(gxx[0],   new_metric.adm_gxx) > rel_tol
  || relative_error(gxy[0],   new_metric.adm_gxy) > rel_tol
  || relative_error(gxz[0],   new_metric.adm_gxz) > rel_tol
  || relative_error(gyy[0],   new_metric.adm_gyy) > rel_tol
  || relative_error(gyz[0],   new_metric.adm_gyz) > rel_tol
  || relative_error(gzz[0],   new_metric.adm_gzz) > rel_tol)
    grhayl_error("grhayl_core_test_suite has failed for GRHayL_enforce_detgtij_and_initialize_metric().\n"
                 "  input metric:  %e %e %e %e %e %e %e %e %e %e\n"
                 "  output metric: %e %e %e %e %e %e %e %e %e %e\n",
                 lapse[0], betax[0], betay[0], betaz[0], 
                 gxx[0], gxy[0], gxz[0], gyy[0], gyz[0], gzz[0],
                 new_metric.lapse, new_metric.betax, new_metric.betay, new_metric.betaz,
                 new_metric.adm_gxx, new_metric.adm_gxy, new_metric.adm_gxz,
                 new_metric.adm_gyy, new_metric.adm_gyz, new_metric.adm_gzz);

  grhayl_info("grhayl_core_test_suite has passed!\n");
}
