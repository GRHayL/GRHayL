#include "unit_tests.h"

int main(int argc, char **argv) {
  const int test_key    = atoi(argv[1]);

  int backup_routine[3] = {None,None,None};
  bool evolve_entropy = false;
  bool evolve_temperature = false;
  bool calc_prims_guess = true;
  double Psi6threshold = 1e100;
  bool Cupp_fix = true;

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
  ghl_initialize_params(None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
                    Psi6threshold, Cupp_fix, Lorenz_damping_factor, &params);

  eos_parameters hybrid_eos;
  ghl_initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_atm, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &hybrid_eos);

  evolve_temperature = true;
  ghl_parameters tab_params;
  ghl_initialize_params(None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
                    Psi6threshold, Cupp_fix, Lorenz_damping_factor, &tab_params);

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

  metric_quantities ADM_metric;
  primitive_quantities prims;
  conservative_quantities cons;
  con2prim_diagnostics diagnostics; 
  ghl_initialize_metric(1.0, 0.0, 0.0, 0.0,
                    1.0, 0.0, 0.0,
                    1.0, 0.0, 1.0,
                    &ADM_metric);

  ADM_aux_quantities metric_aux;
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
      ghl_limit_v_and_compute_u0(&hybrid_eos, &ADM_metric, &prims, &speed_limited);
      break;
    case 5:
      ghl_con2prim_select_method(-5, &params, &hybrid_eos, &ADM_metric, &metric_aux, &cons, &prims, &diagnostics);
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

  eos_parameters tab_eos;
  ghl_initialize_tabulated_eos_functions_and_params(tablepath, W_max,
                                                rho_b_atm, rho_b_min, rho_b_max,
                                                Y_e_atm, Y_e_min, Y_e_max,
                                                T_atm, T_min, T_max, &tab_eos);

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

  printf("We shouldn't be here, so I'll get rid of some compilation warnings :)\n"
         "%e %d %e %e %e %e\n", Fermi_Dirac_integral, speed_limited, rho, Y_e, eps, T);
  return 0;
}
