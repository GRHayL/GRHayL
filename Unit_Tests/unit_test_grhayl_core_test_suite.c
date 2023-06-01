#include "unit_tests.h"

int main(int argc, char **argv) {

  const int neos = 1;
  const double W_max = 10.0;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  eos_parameters hybrid_eos;
  grhayl_initialize_hybrid_eos_functions_and_params(
      W_max,
      rho_b_min, rho_b_min, rho_b_max,
      neos, rho_ppoly, Gamma_ppoly,
      k_ppoly0, Gamma_th,
      &hybrid_eos);

  // TODO: we should add a simple table for tests like this
  //const double Ye_min = 1e-12;
  //const double T_min = 1e-12;
  //const double Ye_max = 1e300;
  //const double T_max = 1e300;
  //eos_parameters tabulated_eos;
  //grhayl_initialize_tabulated_eos_functions_and_params(
  //    W_max,
  //    rho_b_min, rho_b_min, rho_b_max,
  //    Ye_min, Ye_min, Ye_max,
  //    T_min, T_min, T_max,
  //    &tabulated_eos);

  // First test: grhayl_set_prims_to_constant_atm() function
  // Function just sets data to eos data, so no need
  // to store pre-computed comparison data
  primitive_quantities prims;

  for(int eos_it=0; eos_it<2; eos_it++) {
    // poison initial primitive data so reset is obvious
    const double poison = 0.0/0.0;
    prims.rho = poison;
    prims.press = poison;
    prims.vU[0] = poison;
    prims.vU[1] = poison;
    prims.vU[2] = poison;
    prims.eps = poison;
    prims.entropy = poison;
    prims.Y_e = poison;
    prims.temperature = poison;

    eos_parameters eos;
    if(eos_it==0) {
      eos = hybrid_eos;

      grhayl_set_prims_to_constant_atm(&eos, &prims);

      if(prims.rho         != eos.rho_atm
      || prims.press       != eos.press_atm
      || prims.vU[0]          != 0.0
      || prims.vU[1]          != 0.0
      || prims.vU[2]          != 0.0
      || prims.eps         != eos.eps_atm
      || prims.entropy     != eos.entropy_atm)
        grhayl_error("grhayl_core_test_suite has failed for grhayl_set_prims_to_constant_atm() with hybrid EOS.\n"
                     "  rho_b, pressure, vx, vy, vz, epsilon, entropy\n"
                     "  Struct output: %e %e %e %e %e %e %e %e %e\n"
                     "  EOS atm data:  %e %e %e %e %e %e %e %e %e\n",
                     prims.rho, prims.press, prims.vU[0], prims.vU[1], prims.vU[2],
                     prims.eps, prims.entropy, prims.Y_e, prims.temperature,
                     eos.rho_atm, eos.press_atm, 0.0, 0.0, 0.0, eos.eps_atm,
                     eos.entropy_atm, eos.Y_e_atm, eos.T_atm);
  //  } else if(eos_it==0) {
  //    eos = tabulated_eos;

  //    grhayl_set_prims_to_constant_atm(&eos, &prims);

  //    if(prims.rho         != eos.rho_atm
  //    || prims.press       != eos.press_atm
  //    || prims.eps         != eos.eps_atm
  //    || prims.entropy     != eos.entropy_atm
  //    || prims.Y_e         != eos.Ye_atm
  //    || prims.temperature != eos.T_atm
  //    || prims.vU[0]          != 0.0
  //    || prims.vU[1]          != 0.0
  //    || prims.vU[2]          != 0.0)
  //      grhayl_error("grhayl_core_test_suite has failed for grhayl_set_prims_to_constant_atm() with tabulated EOS.\n"
  //                   "  rho_b, pressure, vx, vy, vz, epsilon, entropy, Y_e, temperature\n"
  //                   "  Struct output: %e %e %e %e %e %e %e %e %e\n"
  //                   "  EOS atm data:  %e %e %e %e %e %e %e %e %e\n",
  //                   prims.rho, prims.press, prims.vU[0], prims.vU[1], prims.vU[2],
  //                   prims.eps, prims.entropy, prims.Y_e, prims.temperature,
  //                   eos.rho_atm, eos.press_atm, 0.0, 0.0, 0.0, eos.eps_atm,
  //                   eos.entropy_atm, eos.Ye_atm, eos.T_atm);
    }
  }

  // Second test: grhayl_enforce_detgtij_and_initialize_ADM_metric() function
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
    grhayl_enforce_detgtij_and_initialize_ADM_metric(
        lapse[i],
        betax[i], betay[i], betaz[i],
        gxx[i], gxy[i], gxz[i],
        gyy[i], gyz[i], gzz[i],
        &new_metric);

  if(relative_error(lapse[i], new_metric.lapse) > rel_tol
  || relative_error(betax[i], new_metric.betaU[0]) > rel_tol
  || relative_error(betay[i], new_metric.betaU[1]) > rel_tol
  || relative_error(betaz[i], new_metric.betaU[2]) > rel_tol
  || relative_error(gxx[i],   new_metric.gammaDD[0][0]) > rel_tol
  || relative_error(gxy[i],   new_metric.gammaDD[0][1]) > rel_tol
  || relative_error(gxz[i],   new_metric.gammaDD[0][2]) > rel_tol
  || relative_error(gyy[i],   new_metric.gammaDD[1][1]) > rel_tol
  || relative_error(gyz[i],   new_metric.gammaDD[1][2]) > rel_tol
  || relative_error(gzz[i],   new_metric.gammaDD[2][2]) > rel_tol)
    grhayl_error("unit_test_grhayl_core_test_suite has failed for grhayl_enforce_detgtij_and_initialize_ADM_metric().\n"
                 "  input metric:  %e %e %e %e %e %e %e %e %e %e\n"
                 "  output metric: %e %e %e %e %e %e %e %e %e %e\n",
                 lapse[i], betax[i], betay[i], betaz[i],
                 gxx[i], gxy[i], gxz[i], gyy[i], gyz[i], gzz[i],
                 new_metric.lapse, new_metric.betaU[0], new_metric.betaU[1], new_metric.betaU[2],
                 new_metric.gammaDD[0][0], new_metric.gammaDD[0][1], new_metric.gammaDD[0][2],
                 new_metric.gammaDD[1][1], new_metric.gammaDD[1][2], new_metric.gammaDD[2][2]);

  }

  // For a final test, we trigger the -detg warning
  gxy[0] += 10.0;

    grhayl_enforce_detgtij_and_initialize_ADM_metric(
        lapse[0],
        betax[0], betay[0], betaz[0],
        gxx[0], gxy[0], gxz[0],
        gyy[0], gyz[0], gzz[0],
        &new_metric);

  if(relative_error(lapse[0], new_metric.lapse) > rel_tol
  || relative_error(betax[0], new_metric.betaU[0]) > rel_tol
  || relative_error(betay[0], new_metric.betaU[1]) > rel_tol
  || relative_error(betaz[0], new_metric.betaU[2]) > rel_tol
  || relative_error(gxx[0],   new_metric.gammaDD[0][0]) > rel_tol
  || relative_error(gxy[0],   new_metric.gammaDD[0][1]) > rel_tol
  || relative_error(gxz[0],   new_metric.gammaDD[0][2]) > rel_tol
  || relative_error(gyy[0],   new_metric.gammaDD[1][1]) > rel_tol
  || relative_error(gyz[0],   new_metric.gammaDD[1][2]) > rel_tol
  || relative_error(gzz[0],   new_metric.gammaDD[2][2]) > rel_tol)
    grhayl_error("unit_test_grhayl_core_test_suite has failed for grhayl_enforce_detgtij_and_initialize_ADM_metric().\n"
                 "  input metric:  %e %e %e %e %e %e %e %e %e %e\n"
                 "  output metric: %e %e %e %e %e %e %e %e %e %e\n",
                 lapse[0], betax[0], betay[0], betaz[0],
                 gxx[0], gxy[0], gxz[0], gyy[0], gyz[0], gzz[0],
                 new_metric.lapse, new_metric.betaU[0], new_metric.betaU[1], new_metric.betaU[2],
                 new_metric.gammaDD[0][0], new_metric.gammaDD[0][1], new_metric.gammaDD[0][2],
                 new_metric.gammaDD[1][1], new_metric.gammaDD[1][2], new_metric.gammaDD[2][2]);

  grhayl_info("grhayl_core_test_suite has passed!\n");
}
