// clang-format off
#include "ghl_unit_tests.h"

#ifndef GHL_DISABLE_HDF5
void read_table_error_test(int test_key);
#endif

static ghl_error_codes_t expected_error_code(const int test_key);
static void expect_error_code(ghl_error_codes_t error, int test_key, const char *call);

#ifdef GHL_DISABLE_HDF5
static bool hdf5_only_test(int test_key) {
  return (test_key >= 6 && test_key <= 32)
      || (test_key >= 34 && test_key <= 60)
      || test_key == 63
      || test_key == 66
      || (test_key >= 69 && test_key <= 82);
}
#endif

static void pass_test(int test_key, const char *message) {
  printf("Test %d passed: %s\n", test_key, message);
  exit(1);
}

static void fail_test(int test_key, const char *message) {
  fprintf(stderr, "Test %d failed: %s\n", test_key, message);
  exit(0);
}

int main(int argc, char **argv) {

  const int test_key = atoi(argv[1]);
#ifdef GHL_DISABLE_HDF5
  if(hdf5_only_test(test_key)) {
    pass_test(test_key, "requires HDF5; skipped in no-HDF5 build");
  }
#endif
#ifndef GHL_DISABLE_HDF5
  if((test_key > 33 && test_key < 61) || test_key == 73 || test_key == 74) {
    read_table_error_test(test_key);
    fail_test(test_key, "read_table_error_test returned unexpectedly");
  }
#endif

  ghl_error_codes_t error = ghl_success;

  const ghl_con2prim_id_t None = ghl_con2prim_id_None;
  const ghl_con2prim_id_t backup_routine[3] = {None, None, None};
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
  bool speed_limited = false;

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
  error = ghl_initialize_hybrid_eos_functions_and_params(
        rho_b_atm, rho_b_min, rho_b_max,
        neos, rho_ppoly, Gamma_ppoly,
        k_ppoly0, Gamma_th, &hybrid_eos);
  if(test_key == 0 || test_key == 1) {
    expect_error_code(error, test_key, "ghl_initialize_hybrid_eos_functions_and_params");
  }
  else if(error != ghl_success) {
    fprintf(stderr, "Unexpected setup error %d before test %d\n", error, test_key);
    return 1;
  }

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
      error = NRPyLeakage_Fermi_Dirac_integrals(-1, 1e-4, &Fermi_Dirac_integral);
      expect_error_code(error, test_key, "NRPyLeakage_Fermi_Dirac_integrals");
      break;
    case 3:
      error = NRPyLeakage_Fermi_Dirac_integrals(-1, 1e-2, &Fermi_Dirac_integral);
      expect_error_code(error, test_key, "NRPyLeakage_Fermi_Dirac_integrals");
      break;
  }

  ghl_metric_quantities metric_adm;
  ghl_primitive_quantities prims;
  ghl_conservative_quantities cons;
  ghl_con2prim_diagnostics diagnostics;
  ghl_initialize_metric(
        1.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 0.0, 1.0,
        &metric_adm);

  ghl_ADM_aux_quantities metric_aux;
  ghl_compute_ADM_auxiliaries(&metric_adm, &metric_aux);

  /*
     ghl_limit_v_and_compute_u0:
       4: u^0 is nan
     ghl_con2prim_select_method:
       5: invalid C2P key
  */
  switch (test_key) {
    case 4:
      prims.vU[0] = 0.0/0.0;
      prims.vU[1] = 0.0/0.0;
      prims.vU[2] = 0.0/0.0;
      error = ghl_limit_v_and_compute_u0(&params, &metric_adm, &prims, &speed_limited);
      expect_error_code(error, test_key, "ghl_limit_v_and_compute_u0");
      break;
    case 5:
      error = ghl_con2prim_hybrid_select_method(-5, &params, &hybrid_eos, &metric_adm, &metric_aux, &cons, &prims, &diagnostics);
      expect_error_code(error, test_key, "ghl_con2prim_hybrid_select_method");
      break;
  }
#ifndef GHL_DISABLE_HDF5
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

  /*
     NRPyEOS_read_table_set_EOS_params:
         75: NULL EOS parameter struct
         76: non-tabulated EOS type
         77: invalid table type
     Beta-equilibrium interpolators:
         78: rho too small with ghl_tabulated_compute_Ye_from_rho
         79: rho too large with ghl_tabulated_compute_P_from_rho
         80: P too small with ghl_tabulated_compute_rho_from_P
         81: rho too small with ghl_tabulated_compute_dP_drho_from_rho
         82: rho too large with ghl_tabulated_compute_deps_dP_from_rho
   */
  ghl_eos_parameters tab_eos = { 0 };
  tab_eos.clean_sound_speed = true;
  switch (test_key) {
    case 75:
      error = NRPyEOS_read_table_set_EOS_params(tablepath, NULL);
      expect_error_code(error, test_key, "NRPyEOS_read_table_set_EOS_params");
      break;
    case 76:
      tab_eos.eos_type = ghl_eos_hybrid;
      error = NRPyEOS_read_table_set_EOS_params(tablepath, &tab_eos);
      expect_error_code(error, test_key, "NRPyEOS_read_table_set_EOS_params");
      break;
    case 77:
      tab_eos.eos_type = ghl_eos_tabulated;
      tab_eos.table_type = ghl_eos_table_types;
      error = NRPyEOS_read_table_set_EOS_params(tablepath, &tab_eos);
      expect_error_code(error, test_key, "NRPyEOS_read_table_set_EOS_params");
      break;
  }

  tab_eos.table_type = ghl_eos_table_stellarcollapse;
  error = ghl_initialize_tabulated_eos_functions_and_params(
        tablepath,
        rho_b_atm, rho_b_min, rho_b_max,
        Y_e_atm, Y_e_min, Y_e_max,
        T_atm, T_min, T_max, &tab_eos);
  if(test_key >= 6 && test_key <= 11) {
    expect_error_code(error, test_key, "ghl_initialize_tabulated_eos_functions_and_params");
  }
  else if(error != ghl_success) {
    fprintf(stderr, "Unexpected tabulated EOS setup error %d before test %d\n", error, test_key);
    return 0;
  }

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

  switch (test_key) {
    case 12:
      error = NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(&tab_eos, nvars, rho, Y_e, T, keys, outvars);
      expect_error_code(error, test_key, "NRPyEOS_from_rho_Ye_T_interpolate_n_quantities");
      break;
    case 13:
      rho = rho_b_min-1.0;
      error = NRPyEOS_P_and_eps_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps);
      expect_error_code(error, test_key, "NRPyEOS_P_and_eps_from_rho_Ye_T");
      break;
    case 14:
      rho = rho_b_max+1e2;
      error = NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &S, &cs2);
      expect_error_code(error, test_key, "NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T");
      break;
    case 15:
      Y_e = Y_e_min-1.0;
      error = NRPyEOS_P_eps_and_S_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &S);
      expect_error_code(error, test_key, "NRPyEOS_P_eps_and_S_from_rho_Ye_T");
      break;
    case 16:
      Y_e = Y_e_max+1e2;
      error = NRPyEOS_P_eps_and_cs2_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &cs2);
      expect_error_code(error, test_key, "NRPyEOS_P_eps_and_cs2_from_rho_Ye_T");
      break;
    case 17:
      T = T_min-1.0;
      error = NRPyEOS_eps_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &eps);
      expect_error_code(error, test_key, "NRPyEOS_eps_from_rho_Ye_T");
      break;
    case 18:
      T = T_max+1e5;
      error = NRPyEOS_P_eps_and_depsdT_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &depsdT);
      expect_error_code(error, test_key, "NRPyEOS_P_eps_and_depsdT_from_rho_Ye_T");
      break;
    case 19:
      rho = rho_b_min-1.0;
      error = NRPyEOS_P_eps_muhat_mue_mup_and_mun_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P, &eps, &muhat, &mu_e, &mu_p, &mu_n);
      expect_error_code(error, test_key, "NRPyEOS_P_eps_muhat_mue_mup_and_mun_from_rho_Ye_T");
      break;
    case 20:
      rho = rho_b_min-1.0;
      error = NRPyEOS_P_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &P);
      expect_error_code(error, test_key, "NRPyEOS_P_from_rho_Ye_T");
      break;
    case 21:
      rho = rho_b_min-1.0;
      error = NRPyEOS_eps_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &eps);
      expect_error_code(error, test_key, "NRPyEOS_eps_from_rho_Ye_T");
      break;
    case 22:
      rho = rho_b_min-1.0;
      error = NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T(&tab_eos, rho, Y_e, T, &muhat, &mu_e, &mu_p, &mu_n, &X_n, &X_p);
      expect_error_code(error, test_key, "NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T");
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
                                                               rho, Y_e, eps, NRPyEOS_eps_key, keys, outvars, &T);
      expect_error_code(error, test_key, "NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities");
      break;
    case 24:
      rho = rho_b_min-1.0;
      error = NRPyEOS_P_S_depsdT_and_T_from_rho_Ye_eps(&tab_eos, rho, Y_e, eps, &P, &S, &depsdT, &T);
      expect_error_code(error, test_key, "NRPyEOS_P_S_depsdT_and_T_from_rho_Ye_eps");
      break;
    case 25:
      rho = rho_b_max+1e2;
      error = NRPyEOS_P_and_T_from_rho_Ye_S(&tab_eos, rho, Y_e, S, &P, &T);
      expect_error_code(error, test_key, "NRPyEOS_P_and_T_from_rho_Ye_S");
      break;
    case 26:
      Y_e = Y_e_min-1.0;
      error = NRPyEOS_P_and_T_from_rho_Ye_eps(&tab_eos, rho, Y_e, eps, &P, &T);
      expect_error_code(error, test_key, "NRPyEOS_P_and_T_from_rho_Ye_eps");
      break;
    case 27:
      Y_e = Y_e_max+1e2;
      error = NRPyEOS_P_cs2_and_T_from_rho_Ye_eps(&tab_eos, rho, Y_e, eps, &P, &cs2, &T);
      expect_error_code(error, test_key, "NRPyEOS_P_cs2_and_T_from_rho_Ye_eps");
      break;
    case 28:
      rho = rho_b_min-1.0;
      error = NRPyEOS_P_eps_and_T_from_rho_Ye_S(&tab_eos, rho, Y_e, S, &P, &eps, &T);
      expect_error_code(error, test_key, "NRPyEOS_P_eps_and_T_from_rho_Ye_S");
      break;
    case 29:
      rho = rho_b_min-1.0;
      error = NRPyEOS_T_from_rho_Ye_eps(&tab_eos, rho, Y_e, eps, &T);
      expect_error_code(error, test_key, "NRPyEOS_T_from_rho_Ye_eps");
      break;
    case 30:
      rho = rho_b_min-1.0;
      error = NRPyEOS_eps_S_and_T_from_rho_Ye_P(&tab_eos, rho, Y_e, P, &eps, &S, &T);
      expect_error_code(error, test_key, "NRPyEOS_eps_S_and_T_from_rho_Ye_P");
      break;
    case 31:
      rho = rho_b_min-1.0;
      error = NRPyEOS_eps_cs2_and_T_from_rho_Ye_P(&tab_eos, rho, Y_e, P, &eps, &cs2, &T);
      expect_error_code(error, test_key, "NRPyEOS_eps_cs2_and_T_from_rho_Ye_P");
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
      error = ghl_con2prim_tabulated_select_method(-10, &params, &tab_eos, &metric_adm, &metric_aux, &cons, &prims, &diagnostics);
      expect_error_code(error, test_key, "ghl_con2prim_tabulated_select_method");
      break;
  }
#endif

  switch (test_key) {
    case 33:
    {
      const char *routine_name = ghl_get_con2prim_routine_name(-5);
      if(routine_name == NULL) {
        pass_test(test_key, "ghl_get_con2prim_routine_name returned NULL for invalid key");
      }
      fprintf(stderr,
              "Test %d failed: ghl_get_con2prim_routine_name returned \"%s\" for invalid key\n",
              test_key, routine_name);
      return 1;
    }
  }

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
      error = ghl_initialize_simple_eos_functions_and_params(
            -1, rho_b_min, rho_b_max,
            press_atm, press_min, press_max,
            Gamma_th, &hybrid_eos);
      expect_error_code(error, test_key, "ghl_initialize_simple_eos_functions_and_params");
      break;
    case 62:
      error = ghl_initialize_hybrid_eos_functions_and_params(
            -1, rho_b_min, rho_b_max,
            neos, rho_ppoly, Gamma_ppoly,
            k_ppoly0, Gamma_th, &hybrid_eos);
      expect_error_code(error, test_key, "ghl_initialize_hybrid_eos_functions_and_params");
      break;
#ifndef GHL_DISABLE_HDF5
    case 63:
      error = ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            -1, rho_b_min, rho_b_max,
            Y_e_atm, Y_e_min, Y_e_max,
            T_atm, T_min, T_max, &tab_eos);
      expect_error_code(error, test_key, "ghl_initialize_tabulated_eos_functions_and_params");
      break;
#endif
    case 64:
      error = ghl_initialize_simple_eos_functions_and_params(
            rho_b_atm, rho_b_max, rho_b_min,
            press_atm, press_min, press_max,
            Gamma_th, &hybrid_eos);
      expect_error_code(error, test_key, "ghl_initialize_simple_eos_functions_and_params");
      break;
    case 65:
      error = ghl_initialize_hybrid_eos_functions_and_params(
            rho_b_atm, rho_b_max, rho_b_min,
            neos, rho_ppoly, Gamma_ppoly,
            k_ppoly0, Gamma_th, &hybrid_eos);
      expect_error_code(error, test_key, "ghl_initialize_hybrid_eos_functions_and_params");
      break;
#ifndef GHL_DISABLE_HDF5
    case 66:
      error = ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_max, rho_b_min,
            Y_e_atm, Y_e_min, Y_e_max,
            T_atm, T_min, T_max, &tab_eos);
      expect_error_code(error, test_key, "ghl_initialize_tabulated_eos_functions_and_params");
      break;
#endif
    case 67:
      error = ghl_initialize_simple_eos_functions_and_params(
            rho_b_atm, rho_b_min, rho_b_max,
            -1, press_min, press_max,
            Gamma_th, &hybrid_eos);
      expect_error_code(error, test_key, "ghl_initialize_simple_eos_functions_and_params");
      break;
    case 68:
      error = ghl_initialize_simple_eos_functions_and_params(
            rho_b_atm, rho_b_min, rho_b_max,
            press_atm, press_max, press_min,
            Gamma_th, &hybrid_eos);
      expect_error_code(error, test_key, "ghl_initialize_simple_eos_functions_and_params");
      break;
#ifndef GHL_DISABLE_HDF5
    case 69:
      error = ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_min, rho_b_max,
            -1, Y_e_min, Y_e_max,
            T_atm, T_min, T_max, &tab_eos);
      expect_error_code(error, test_key, "ghl_initialize_tabulated_eos_functions_and_params");
      break;
    case 70:
      error = ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_min, rho_b_max,
            Y_e_atm, Y_e_max, Y_e_min,
            T_atm, T_min, T_max, &tab_eos);
      expect_error_code(error, test_key, "ghl_initialize_tabulated_eos_functions_and_params");
      break;
    case 71:
      error = ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_min, rho_b_max,
            Y_e_atm, Y_e_min, Y_e_max,
            -1, T_min, T_max, &tab_eos);
      expect_error_code(error, test_key, "ghl_initialize_tabulated_eos_functions_and_params");
      break;
    case 72:
      error = ghl_initialize_tabulated_eos_functions_and_params(
            tablepath,
            rho_b_atm, rho_b_min, rho_b_max,
            Y_e_atm, Y_e_min, Y_e_max,
            T_atm, T_max, T_min, &tab_eos);
      expect_error_code(error, test_key, "ghl_initialize_tabulated_eos_functions_and_params");
      break;
#endif
#ifndef GHL_DISABLE_HDF5
    case 78:
    case 79:
    case 80:
    case 81:
    case 82:
      error = ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(exp(4.0), &tab_eos);
      if(error != ghl_success) {
        char message[128];
        snprintf(message, sizeof(message), "Unexpected tabulated EOS setup error %d", error);
        fail_test(test_key, message);
      }
      break;
#endif
  }

#ifndef GHL_DISABLE_HDF5
  switch (test_key) {
    case 78:
      error = ghl_tabulated_compute_Ye_from_rho(&tab_eos, 0.5 * tab_eos.table_rho_min, &Y_e);
      expect_error_code(error, test_key, "ghl_tabulated_compute_Ye_from_rho");
      break;
    case 79:
      error = ghl_tabulated_compute_P_from_rho(&tab_eos, 2.0 * tab_eos.table_rho_max, &P);
      expect_error_code(error, test_key, "ghl_tabulated_compute_P_from_rho");
      break;
    case 80:
      error = ghl_tabulated_compute_P_from_rho(&tab_eos, tab_eos.table_rho_min, &P);
      if(error != ghl_success) {
        fprintf(stderr, "Failed to compute minimum pressure before test %d: %d\n", test_key, error);
        return 1;
      }
      error = ghl_tabulated_compute_rho_from_P(&tab_eos, 0.5 * P, &rho);
      expect_error_code(error, test_key, "ghl_tabulated_compute_rho_from_P");
      break;
    case 81:
      error = ghl_tabulated_compute_dP_drho_from_rho(&tab_eos, 0.5 * tab_eos.table_rho_min, &depsdT);
      expect_error_code(error, test_key, "ghl_tabulated_compute_dP_drho_from_rho");
      break;
    case 82:
      error = ghl_tabulated_compute_deps_dP_from_rho(&tab_eos, 2.0 * tab_eos.table_rho_max, &depsdT);
      expect_error_code(error, test_key, "ghl_tabulated_compute_deps_dP_from_rho");
      break;
  }
#endif

  // Silence warnings
  (void)Fermi_Dirac_integral;
  (void)speed_limited;
#ifndef GHL_DISABLE_HDF5
  (void)rho;
  (void)Y_e;
  (void)eps;
  (void)T;
#endif

  fail_test(test_key, "did not trigger its expected failure path");
}
// clang-format on

#ifndef GHL_DISABLE_HDF5
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

void create_scalar_dataset(char *name, hid_t datatype_id, hid_t file_id) {
  hid_t dataspace_id = H5Screate(H5S_SCALAR);
  hid_t dataset_id = H5Dcreate2(file_id, name, datatype_id, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  int data = 1;
  H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}

void create_opaque_dataset(char *name, hid_t file_id) {
  hsize_t dims[1] = {1};
  hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
  hid_t datatype_id = H5Tcreate(H5T_OPAQUE, 1);
  hid_t dataset_id = H5Dcreate2(file_id, name, datatype_id, dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  char data = 1;
  H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
  H5Dclose(dataset_id);
  H5Tclose(datatype_id);
  H5Sclose(dataspace_id);
}

#endif

static ghl_error_codes_t expected_error_code(const int test_key) {
  switch(test_key) {
    case  0:
    case  6:
    case 61:
    case 62:
    case 63: return ghl_error_invalid_rho_atm;
    case  1:
    case  9:
    case 64:
    case 65:
    case 66: return ghl_error_rho_min_gt_rho_max;
    case  2:
    case  3: return ghl_error_invalid_fermi_dirac_integral_key;
    case  4: return ghl_error_u0_singular;
    case  5:
    case 32: return ghl_error_invalid_c2p_key;
    case  7:
    case 69: return ghl_error_invalid_Y_e_atm;
    case  8:
    case 71: return ghl_error_invalid_T_atm;
    case 10:
    case 70: return ghl_error_Y_e_min_gt_Y_e_max;
    case 11:
    case 72: return ghl_error_T_min_gt_T_max;
    case 12:
    case 23: return ghl_error_exceed_table_vars;
    case 13:
    case 19:
    case 20:
    case 21:
    case 22:
    case 24:
    case 28:
    case 29:
    case 30:
    case 31: return ghl_error_table_min_rho;
    case 14:
    case 25: return ghl_error_table_max_rho;
    case 15:
    case 26: return ghl_error_table_min_ye;
    case 16:
    case 27: return ghl_error_table_max_ye;
    case 17: return ghl_error_table_min_T;
    case 18: return ghl_error_table_max_T;
    case 34: return ghl_error_could_not_open_file;
    case 35:
    case 36:
    case 37:
    case 38:
    case 39:
    case 40:
    case 41:
    case 42:
    case 43:
    case 44:
    case 45:
    case 46:
    case 47:
    case 48:
    case 49:
    case 50:
    case 51:
    case 52:
    case 53:
    case 54:
    case 55:
    case 56:
    case 57:
    case 58:
    case 59:
    case 60: return ghl_error_hdf5_dataset_could_not_open;
    case 67: return ghl_error_invalid_press_atm;
    case 68: return ghl_error_press_min_gt_press_max;
    case 73: return ghl_error_hdf5_dataset_could_not_open;
    case 74: return ghl_error_hdf5_dataset_could_not_read;
    case 75: return ghl_error_eos_struct_is_null;
    case 76: return ghl_error_invalid_eos_type;
    case 77: return ghl_error_invalid_eos_table_type;
    case 78:
    case 79:
    case 80:
    case 81:
    case 82: return ghl_error_root_not_bracketed;
  }

  fprintf(stderr, "No expected error code configured for test %d\n", test_key);
  exit(0);
}

static void expect_error_code(ghl_error_codes_t error, int test_key, const char *call) {
  const ghl_error_codes_t expected = expected_error_code(test_key);
  if(error == ghl_success) {
    fprintf(stderr,
            "Test %d failed: %s returned ghl_success, expected error code %d\n",
            test_key, call, expected);
    exit(0);
  }

  if(error != expected) {
    fprintf(stderr,
            "Test %d failed: %s returned error code %d, expected %d\n",
            test_key, call, error, expected);
    exit(0);
  }

  printf("Test %d passed: %s returned expected error code %d\n",
         test_key, call, error);
  ghl_abort_if_error(error);
}

#ifndef GHL_DISABLE_HDF5

void read_table_error_test(int test_key) {

  printf("read_table_error_test: %d\n", test_key);

  int original_test_key = test_key;
  test_key -= 34;
  FILE *fp = fopen("test.h5", "r");
  if(fp) {
    fclose(fp);
    remove("test.h5");
  }
  if(test_key==0) {
    ghl_eos_parameters eos = { 0 };
    eos.eos_type = ghl_eos_tabulated;
    eos.table_type = ghl_eos_table_stellarcollapse;
    eos.clean_sound_speed = true;
    ghl_error_codes_t error = NRPyEOS_read_table_set_EOS_params("test.h5", &eos);
    expect_error_code(error, original_test_key, "NRPyEOS_read_table_set_EOS_params");
  }
  test_key--;
  if(original_test_key == 73 || original_test_key == 74) {
    test_key = original_test_key;
  }

  if (test_key < 0 || (test_key > 25 && test_key != 73 && test_key != 74)) {
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
  if(test_key < 73) {
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
  }
  else if(test_key == 73) {
    // ndims failure test, line 21 of NRPyEOS_hdf5_helpers.c
    create_scalar_dataset(binnames[0], H5T_NATIVE_INT, file_id);
  }
  else if(test_key == 74) {
    // status < 0 failure test, line 39 of current version of NRPyEOS_hdf5_helpers.c
    create_opaque_dataset(binnames[0], file_id);
  }
  H5Fclose(file_id);

  ghl_eos_parameters eos = { 0 };
  eos.eos_type = ghl_eos_tabulated;
  eos.table_type = ghl_eos_table_stellarcollapse;
  eos.clean_sound_speed = true;
  ghl_error_codes_t error = NRPyEOS_read_table_set_EOS_params("test.h5", &eos);
  expect_error_code(error, original_test_key, "NRPyEOS_read_table_set_EOS_params");
}

#endif
