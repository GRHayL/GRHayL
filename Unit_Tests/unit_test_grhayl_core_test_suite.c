#include "ghl_unit_tests.h"

int main(int argc, char **argv) {

  const double rel_tol = 1e-14;

  const int neos = 1;
  const double rho_b_min = 1e-12;
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;

  // Checking that warnings trigger
  const double press_min = 1e-20;
  ghl_eos_parameters simple_eos;
  ghl_initialize_simple_eos_functions_and_params(
      rho_b_min, -1, -1,
      press_min, -1, -1,
      Gamma_th, &simple_eos);
  if(simple_eos.rho_min > 1e-50)
    ghl_error("Simple EOS failed to set default rho_min");
  if(simple_eos.rho_max < 9e299)
    ghl_error("Simple EOS failed to set default rho_max");
  if(simple_eos.press_min > 1e-50)
    ghl_error("Simple EOS failed to set default press_min");
  if(simple_eos.press_max < 9e299)
    ghl_error("Simple EOS failed to set default press_max");

  ghl_eos_parameters hybrid_eos;
  ghl_initialize_hybrid_eos_functions_and_params(
      rho_b_min, rho_b_min, rho_b_max,
      neos, rho_ppoly, Gamma_ppoly,
      k_ppoly0, Gamma_th,
      &hybrid_eos);

  //const double Ye_atm = 1e-12;
  //const double T_atm = 1e300;
  //ghl_eos_parameters tab_eos;
  //ghl_initialize_tabulated_eos_functions_and_params(
  //    table_filepath,
  //    rho_b_min, -1, -1,
  //    Ye_atm, -1, -1,
  //    T_atm, -1, -1,
  //    &tab_eos)
  ////TODO: add tests checking defaults based on table

  // ghl_set_prims_to_constant_atm() function
  // Function just sets data to eos data, so no need
  // to store pre-computed comparison data
  ghl_primitive_quantities prims;

  // TODO: when tabulated is added, add 1 to upper limit of loop
  for(int eos_type = 0; eos_type < 2; eos_type++) {
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

    ghl_eos_parameters eos;
    char *eos_name;
    switch(eos_type) {
      case 0:
        eos = simple_eos;
        eos_name = "simple";
        break;
      case 1:
        eos = hybrid_eos;
        eos_name = "hybrid";
        break;
      case 2:
        //eos = tab_eos;
        eos_name = "tabulated";
        break;
      default:
        exit(1);
    }

    ghl_set_prims_to_constant_atm(&eos, &prims);

    bool check1 = (fabs(prims.rho     - eos.rho_atm) > 1e-50
                || fabs(prims.press   - eos.press_atm) > 1e-50
                || fabs(prims.vU[0]) > 1e-50
                || fabs(prims.vU[1]) > 1e-50
                || fabs(prims.vU[2]) > 1e-50
                || fabs(prims.eps     - eos.eps_atm) > 1e-50
                || fabs(prims.entropy - eos.entropy_atm) > 1e-50);

    bool check2 = false;
    if(eos_type ==2)
      check2 = (fabs(prims.Y_e         - eos.Y_e_atm) > 1e-50
             || fabs(prims.temperature - eos.T_atm) > 1e-50);

    if(check1 || check2)
      ghl_error("grhayl_core_test_suite has failed for ghl_set_prims_to_constant_atm() with %s EOS.\n"
                   "  rho_b, pressure, vx, vy, vz, epsilon, entropy\n"
                   "  Struct output: %e %e %e %e %e %e %e %e %e\n"
                   "  EOS atm data:  %e %e %e %e %e %e %e %e %e\n",
                   eos_name, prims.rho, prims.press, prims.vU[0], prims.vU[1], prims.vU[2],
                   prims.eps, prims.entropy, prims.Y_e, prims.temperature,
                   eos.rho_atm, eos.press_atm, 0.0, 0.0, 0.0, eos.eps_atm,
                   eos.entropy_atm, eos.Y_e_atm, eos.T_atm);
  }

  // ghl_compute_vec2_from_vec4D() function
  // Use test quantities with known result
  double g4[4][4];
  g4[0][0]            =  1.0;
  g4[0][1] = g4[1][0] =  2.0;
  g4[0][2] = g4[2][0] =  3.0;
  g4[0][3] = g4[3][0] =  4.0;
  g4[1][1]            =  5.0;
  g4[1][2] = g4[2][1] =  6.0;
  g4[1][3] = g4[3][1] =  7.0;
  g4[2][2]            =  8.0;
  g4[2][3] = g4[3][2] =  9.0;
  g4[3][3]            = 10.0;

  const double vec[4] = {-1.0, 2.0, -3.0, 4.0};

  const double v2 = ghl_compute_vec2_from_vec4D(g4, vec);
  if(relative_error(v2, 55.0) > rel_tol) 
    ghl_error("unit_test_grhayl_core_test_suite has failed for ghl_compute_vec2_from_vec4D().\n"
              "  Expected result is 55, but computed value is %e.\n",
              v2);

  double vec_inv[4];
  ghl_raise_lower_vector_4D(g4, vec, vec_inv);
  const bool check_inverse = (relative_error(vec_inv[0], 10.0) > rel_tol
                           || relative_error(vec_inv[1], 18.0) > rel_tol
                           || relative_error(vec_inv[2], 21.0) > rel_tol
                           || relative_error(vec_inv[3], 23.0) > rel_tol);
  if(check_inverse) 
    ghl_error("unit_test_grhayl_core_test_suite has failed for ghl_raise_lower_vector_4D().\n"
              "  Expected values are 10, 18, 21, 23. Computed values are %e %e %e %e.\n",
              vec_inv[0], vec_inv[1], vec_inv[2], vec_inv[3]);

  // ghl_enforce_detgtij_and_initialize_ADM_metric() function
  // Valid metrics should return near-identical (round-off level) metrics
  ghl_metric_quantities new_metric;

  // read in metric data
  FILE* infile = fopen_with_check("grhayl_core_test_suite_input.bin","rb");

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
    ghl_error("An error has occured with reading in initial data. Please check that data\n"
                 "is up-to-date with current test version.\n");

  for(int i=0; i<arraylength; i++) {
    ghl_enforce_detgtij_and_initialize_ADM_metric(
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
    ghl_error("unit_test_grhayl_core_test_suite has failed for ghl_enforce_detgtij_and_initialize_ADM_metric().\n"
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

    ghl_enforce_detgtij_and_initialize_ADM_metric(
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
    ghl_error("unit_test_grhayl_core_test_suite has failed for ghl_enforce_detgtij_and_initialize_ADM_metric().\n"
              "  input metric:  %e %e %e %e %e %e %e %e %e %e\n"
              "  output metric: %e %e %e %e %e %e %e %e %e %e\n",
              lapse[0], betax[0], betay[0], betaz[0],
              gxx[0], gxy[0], gxz[0], gyy[0], gyz[0], gzz[0],
              new_metric.lapse, new_metric.betaU[0], new_metric.betaU[1], new_metric.betaU[2],
              new_metric.gammaDD[0][0], new_metric.gammaDD[0][1], new_metric.gammaDD[0][2],
              new_metric.gammaDD[1][1], new_metric.gammaDD[1][2], new_metric.gammaDD[2][2]);

  char *valid_char[12];
  valid_char[0] = "None";
  valid_char[1] = "Noble2D";
  valid_char[2] = "Noble1D";
  valid_char[3] = "Noble1D_entropy";
  valid_char[4] = "Noble1D_entropy2";
  valid_char[5] = "Font1D";
  valid_char[6] = "CerdaDuran2D";
  valid_char[7] = "CerdaDuran3D";
  valid_char[8] = "Palenzuela1D";
  valid_char[9] = "Palenzuela1D_entropy";
  valid_char[10] = "Newman1D";
  valid_char[11] = "Newman1D_entropy";
  for(int i = None; i <= Newman1D_entropy; i++) {
    const char *test_char = ghl_get_con2prim_routine_name(i);
    if(strcmp(valid_char[i+1], test_char)) {
      ghl_error("ghl_get_con2prim_routine_name failed:\n    expected %s, got %s\n", valid_char[i+1], test_char);
    }
  }

  ghl_info("grhayl_core_test_suite has passed!\n");
}
