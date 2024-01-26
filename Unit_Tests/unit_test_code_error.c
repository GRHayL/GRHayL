// clang-format off
#include "unit_tests.h"

void read_table_error_test(int test_key);

int main(int argc, char **argv) {
  const int test_key    = atoi(argv[1]);


  int backup_routine[3] = {None,None,None};
  bool evolve_entropy = false;
  bool evolve_temperature = false;
  bool calc_prims_guess = true;
  double Psi6threshold = 1e100;

  int neos = 1;
  double W_max = 10.0;
  double rho_b_atm = 1e-12;
  double rho_b_min = 1e-12;
  double rho_b_max = 1e300;
  double Gamma_th = 2.0;
  double rho_ppoly[1] = {0.0};
  double Gamma_ppoly[1] = {2.0};
  double k_ppoly0 = 1.0;
  double Lorenz_damping_factor = 0.0;

  /*
     Hybrid EOS setup:
       0: rho_b_atm too small
       1: rho_b_min > rho_b_max
  */
  switch (test_key) {
    case 0:
      rho_b_atm = -1;
      break;
    case 1:
      rho_b_min = 1e301;
      break;
  }

  ghl_parameters params;
  ghl_initialize_params(
        None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
        Psi6threshold, W_max, Lorenz_damping_factor, &params);

  ghl_eos_parameters hybrid_eos;
  ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_atm, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &hybrid_eos);

  evolve_temperature = true;
  ghl_parameters tab_params;
  ghl_initialize_params(
        None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
        Psi6threshold, W_max, Lorenz_damping_factor, &tab_params);

  /*
     NRPyLeakage_Fermi_Dirac_integrals:
       2: invalid choice of k (z<1e-3)
       3: invalid choice of k (z>1e-3)
  */
  double Fermi_Dirac_integral = 0.0;
  switch (test_key) {
    case 2:
      Fermi_Dirac_integral = NRPyLeakage_Fermi_Dirac_integrals(-1, 1e-4);
      break;
    case 3:
      Fermi_Dirac_integral = NRPyLeakage_Fermi_Dirac_integrals(-1, 1e-2);
      break;
  }

  ghl_metric_quantities ADM_metric;
  ghl_primitive_quantities prims;
  ghl_conservative_quantities cons;
  ghl_con2prim_diagnostics diagnostics;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        &ADM_metric);

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

  /*
     ghl_limit_v_and_compute_u0:
       4: u^0 is nan
     ghl_con2prim_select_method:
       5: invalid C2P key
  */
   int speed_limited = 0;
  switch (test_key) {
    case 4:
      prims.vU[0] = 0.0/0.0;
      prims.vU[1] = 0.0/0.0;
      prims.vU[2] = 0.0/0.0;
      speed_limited = ghl_limit_v_and_compute_u0(&params, &ADM_metric, &prims);
      break;
    case 5:
      if(ghl_con2prim_hybrid_select_method(-5, &params, &hybrid_eos, &ADM_metric, &metric_aux, &cons, &prims, &diagnostics) < 0)
        ghl_Error(100, "Unsupported c2p key (%d) with hybrid EOS.\n", -5);
      break;
  }
/*
rho: 2.718281828459045e+00, 1.096633158428459e+03
 T : 2.718281828459045e+00, 1.484131591025766e+02
Y_e: 1.000000000000000e+00, 3.000000000000000e+00
*/
  const char tablepath[] = "SLy4_3335_rho391_temp163_ye66.h5";
  double Y_e_atm   = 0.5;
  double Y_e_min   = 0.05;
  double Y_e_max   = Y_e_atm;
  double T_atm     = 1e-2;
  double T_min     = T_atm;
  double T_max     = 1e2;

  /*
     Tabulated EOS setup:
       6: rho_b_atm too small
       7: Y_e_atm too small
       8: T_atm too small
       9: rho_b_min > rho_b_max
      10: Y_e_min > Y_e_max
      11: T_min > T_max
   */
  switch (test_key) {
    case 6:
      rho_b_atm = -1;
      break;
    case 7:
      Y_e_atm = -1;
      break;
    case 8:
    T_atm = -1;
      break;
    case 9:
      rho_b_min = 1e301;
      break;
    case 10:
      Y_e_min = 1e1;
      break;
    case 11:
      T_min = 1e3;
      break;
  }

  ghl_eos_parameters tab_eos;
  /*
     NRPyEOS interpolators:
       NRPyEOS_from_rho_Ye_T_interpolate_n_quantities:
         12: TODO: still need to manually call w/ too many keys
         13: rho is too small with NRPyEOS_P_and_eps_from_rho_Ye_T
         14: rho is too large with NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T
         15: Y_e is too small with NRPyEOS_P_eps_and_S_from_rho_Ye_T
         16: Y_e is too large with NRPyEOS_P_eps_and_cs2_from_rho_Ye_T
         17: T is too small   with NRPyEOS_eps_from_rho_Ye_T
         18: T is too large   with NRPyEOS_P_eps_and_depsdT_from_rho_Ye_T
       Following tests check remaining interpolators with rho too small
         19: NRPyEOS_P_eps_muhat_mue_mup_and_mun_from_rho_Ye_T
         20: NRPyEOS_P_from_rho_Ye_T
         21: NRPyEOS_eps_from_rho_Ye_T
         22: NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T
   */
  double rho, Y_e, T, eps, S;
  rho = 1e-2;
  Y_e = eps = T = S = 0.1;
  double P, cs2, depsdT, muhat, mu_e, mu_p, mu_n, X_n, X_p;
  const int nvars = 30;
  const int keys[1];
  double outvars[1];
  NRPyEOS_error_report report;
  int error;

  if(test_key < 31) {
    ghl_initialize_tabulated_eos_functions_and_params(
          tablepath,
          rho_b_atm, rho_b_min, rho_b_max,
          Y_e_atm, Y_e_min, Y_e_max,
          T_atm, T_min, T_max, &tab_eos);
  }

  switch (test_key) {
    case 12:
      error = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(&tab_eos, nvars, rho, Y_e, T, keys, outvars, &report);
      if( error )
        ghl_error(report.message, error);
      break;
    case 13:
      rho = rho_b_min-1.0;
      NRPyEOS_P_and_eps_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps);
      break;
    case 14:
      rho = rho_b_max+1e2;
      NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &S, &cs2);
      break;
    case 15:
      Y_e = Y_e_min-1.0;
      NRPyEOS_P_eps_and_S_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &S);
      break;
    case 16:
      Y_e = Y_e_max+1e2;
      NRPyEOS_P_eps_and_cs2_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &cs2);
      break;
    case 17:
      T = T_min-1.0;
      NRPyEOS_eps_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &eps);
      break;
    case 18:
      T = T_max+1e5;
      NRPyEOS_P_eps_and_depsdT_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &depsdT);
      break;
    case 19:
      rho = rho_b_min-1.0;
      NRPyEOS_P_eps_muhat_mue_mup_and_mun_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &muhat, &mu_e, &mu_p, &mu_n);
      break;
    case 20:
      rho = rho_b_min-1.0;
      NRPyEOS_P_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P);
      break;
    case 21:
      rho = rho_b_min-1.0;
      NRPyEOS_eps_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &eps);
      break;
    case 22:
      rho = rho_b_min-1.0;
      NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &muhat, &mu_e, &mu_p, &mu_n, &X_n, &X_p);
      break;
  }

  /*
     NRPyEOS interpolators:
       NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities:
         23: TODO: still need to manually call w/ too many keys
         24: rho is too small with NRPyEOS_P_S_depsdT_and_T_from_rho_Ye_eps
         25: rho is too large with NRPyEOS_P_and_T_from_rho_Ye_S
         26: Y_e is too small with NRPyEOS_P_and_T_from_rho_Ye_eps
         27: Y_e is too large with NRPyEOS_P_cs2_and_T_from_rho_Ye_eps
       Following tests check remaining interpolators with rho too small
         28: NRPyEOS_P_eps_and_T_from_rho_Ye_S
         29: NRPyEOS_T_from_rho_Ye_eps
         30: NRPyEOS_eps_S_and_T_from_rho_Ye_P
         31: NRPyEOS_eps_cs2_and_T_from_rho_Ye_P
   */
  switch (test_key) {
    case 23:
      error = NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(&tab_eos, nvars, tab_eos.root_finding_precision,
                                                               rho, Y_e, eps, NRPyEOS_eps_key, keys, outvars, &T, &report);
      if( error )
        ghl_error(report.message, error);
      break;
    case 24:
      rho = rho_b_min-1.0;
      NRPyEOS_P_S_depsdT_and_T_from_rho_Ye_eps(&tab_eos, rho, Y_e, eps, &P, &S, &depsdT, &T);
      break;
    case 25:
      rho = rho_b_max+1e2;
      NRPyEOS_P_and_T_from_rho_Ye_S(&tab_eos, rho, Y_e, S, &P, &T);
      break;
    case 26:
      Y_e = Y_e_min-1.0;
      NRPyEOS_P_and_T_from_rho_Ye_eps(&tab_eos, rho, Y_e, eps, &P, &T);
      break;
    case 27:
      Y_e = Y_e_max+1e2;
      NRPyEOS_P_cs2_and_T_from_rho_Ye_eps(&tab_eos, rho, Y_e, eps, &P, &cs2, &T);
      break;
    case 28:
      rho = rho_b_min-1.0;
      NRPyEOS_P_eps_and_T_from_rho_Ye_S(&tab_eos, rho, Y_e, S, &P, &eps, &T);
      break;
    case 29:
      rho = rho_b_min-1.0;
      NRPyEOS_T_from_rho_Ye_eps(&tab_eos, rho, Y_e, eps, &T);
      break;
    case 30:
      rho = rho_b_min-1.0;
      NRPyEOS_eps_S_and_T_from_rho_Ye_P(&tab_eos, rho, Y_e, P, &eps, &S, &T);
      break;
    case 31:
      rho = rho_b_min-1.0;
      NRPyEOS_eps_cs2_and_T_from_rho_Ye_P(&tab_eos, rho, Y_e, P, &eps, &cs2, &T);
      break;
  }

  /*
     ghl_con2prim_*_select_method:
         32: Invalid tabulated C2P key
     ghl_get_con2prim_routine_name:
         33: Invalid C2P key
   */
  
  switch (test_key) {
    case 32:
      if(ghl_con2prim_tabulated_select_method(-10, &params, &tab_eos, &ADM_metric, &metric_aux, &cons, &prims, &diagnostics) < 0)
        ghl_Error(100, "Unsupported c2p key (%d) with hybrid EOS.\n", -10);
      break;
    case 33:
      printf("%s\n", ghl_get_con2prim_routine_name(-5));
      break;
  }

  if(test_key > 33 && test_key < 61) {
    read_table_error_test(test_key);

  /*
     initialize_*_eos_functions_and_params::
         61-63: Invalid rho_atm
         64-66: rho_min > rho_max
         67: Invalid P_atm
         68: P_min > P_max
         69: Invalid Y_e_atm
         70: Y_e_min > Y_e_max
         71: Invalid T_atm
         72: T_min > T_max
   */

  const double press_min = 1e-20;
  const double press_atm = 1e-16;
  const double press_max = 1e4;
  switch (test_key) {
    case 61:
      ghl_initialize_simple_eos_functions_and_params(
            -1, rho_b_min, rho_b_max,
            press_atm, press_min, press_max,
            Gamma_th, &hybrid_eos);
      break;
    case 62:
      ghl_initialize_hybrid_eos_functions_and_params(
            -1, rho_b_min, rho_b_max,
            neos, rho_ppoly, Gamma_ppoly,
            k_ppoly0, Gamma_th, &hybrid_eos);
      break;
    case 63:
      ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            -1, rho_b_min, rho_b_max,
            Y_e_atm, Y_e_min, Y_e_max,
            T_atm, T_min, T_max, &tab_eos);
      break;
    case 64:
      ghl_initialize_simple_eos_functions_and_params(
            rho_b_atm, rho_b_max, rho_b_min,
            press_atm, press_min, press_max,
            Gamma_th, &hybrid_eos);
      break;
    case 65:
      ghl_initialize_hybrid_eos_functions_and_params(
            rho_b_atm, rho_b_max, rho_b_min,
            neos, rho_ppoly, Gamma_ppoly,
            k_ppoly0, Gamma_th, &hybrid_eos);
      break;
    case 66:
      ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_max, rho_b_min,
            Y_e_atm, Y_e_min, Y_e_max,
            T_atm, T_min, T_max, &tab_eos);
      break;
    case 67:
      ghl_initialize_simple_eos_functions_and_params(
            rho_b_atm, rho_b_min, rho_b_max,
            -1, press_min, press_max,
            Gamma_th, &hybrid_eos);
      break;
    case 68:
      ghl_initialize_simple_eos_functions_and_params(
            rho_b_atm, rho_b_min, rho_b_max,
            press_atm, press_max, press_min,
            Gamma_th, &hybrid_eos);
      break;
    case 69:
      ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_min, rho_b_max,
            -1, Y_e_min, Y_e_max,
            T_atm, T_min, T_max, &tab_eos);
      break;
    case 70:
      ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_min, rho_b_max,
            Y_e_atm, Y_e_max, Y_e_min,
            T_atm, T_min, T_max, &tab_eos);
      break;
    case 71:
      ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_min, rho_b_max,
            Y_e_atm, Y_e_min, Y_e_max,
            -1, T_min, T_max, &tab_eos);
      break;
    case 72:
      ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_min, rho_b_max,
            Y_e_atm, Y_e_min, Y_e_max,
            T_atm, T_max, T_min, &tab_eos);
      break;
  }

  printf("We shouldn't be here, so I'll get rid of some compilation warnings :)\n"
         "%e %d %e %e %e %e\n", Fermi_Dirac_integral, speed_limited, rho, Y_e, eps, T);
  return 0;
}
// clang-format on

void create_dataset(char *name, int dim, hid_t datatype_id, hid_t file_id) {
  hsize_t dims[dim];
  for(int i=0;i<dim;i++) {
    dims[i] = 1;
  }
  hid_t dataspace_id = H5Screate_simple(dim, dims, NULL);
  hid_t dataset_id = H5Dcreate2(file_id, name, datatype_id, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if( datatype_id == H5T_NATIVE_INT ) {
    int data = 1;
    H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
  }
  else {
    double data = 42;
    H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
  }
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}

void read_table_error_test(int test_key) {

  test_key -= 34;
  FILE *fp = fopen("test.h5", "r");
  if(fp) {
    fclose(fp);
    remove("test.h5");
  }
  if(test_key==0) {
    ghl_eos_parameters eos;
    NRPyEOS_read_table_set_EOS_params("test.h5", &eos);
  }
  test_key--;

  if (test_key < 0 || test_key > 25) {
    ghl_info("read_table_error_test: unknown test_key\n");
    return;
  }

  char *binnames[] = {"pointsrho", "pointstemp",  "pointsye", "logpress",
                      "logenergy", "entropy",     "munu",     "cs2",
                      "dedt",      "dpdrhoe",     "dpderho",  "muhat",
                      "mu_e",      "mu_p",        "mu_n",     "Xa",
                      "Xh",        "Xn",          "Xp",       "Abar",
                      "Zbar",      "gamma",       "logrho",   "logtemp",
                      "ye",        "energy_shift"};

  hid_t file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  int imax = test_key > 3 ? 3 : test_key;
  for (int i = 0; i < imax; i++) {
    create_dataset(binnames[i], 1, H5T_NATIVE_INT, file_id);
  }
  imax = test_key > 22 ? 22 : test_key;
  for (int i = 3; i < imax; i++) {
    create_dataset(binnames[i], 3, H5T_NATIVE_DOUBLE, file_id); 
  }
  for(int i=22;i<test_key;i++) {
    create_dataset(binnames[i], 1, H5T_NATIVE_DOUBLE, file_id);
  }
  H5Fclose(file_id);

  ghl_eos_parameters eos;
  NRPyEOS_read_table_set_EOS_params("test.h5", &eos);
}
