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

  if(test_key == 0) {
    rho_b_atm = -1;
  } else if(test_key == 1) {
    rho_b_min = 1e301;
  }

  GRHayL_parameters hybrid_params;
  grhayl_initialize_params(None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
                    Psi6threshold, Cupp_fix, Lorenz_damping_factor, &hybrid_params);

  eos_parameters hybrid_eos;
  initialize_hybrid_eos_functions_and_params(W_max,
                                             rho_b_atm, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &hybrid_eos);

  evolve_temperature = true;
  GRHayL_parameters tab_params;
  grhayl_initialize_params(None, backup_routine, evolve_entropy, evolve_temperature, calc_prims_guess,
                    Psi6threshold, Cupp_fix, Lorenz_damping_factor, &tab_params);

  double Fermi_Dirac_integral = 0.0;
  if(test_key == 8) {
    int k = -1;
    double z = 1e-4;
    Fermi_Dirac_integral = NRPyLeakage_Fermi_Dirac_integrals(k, z);
  } else if(test_key == 9) {
    int k = -1;
    double z = 1e-2;
    Fermi_Dirac_integral = NRPyLeakage_Fermi_Dirac_integrals(k, z);
  }

  metric_quantities metric;
  initialize_metric(1.0,
                    1.0, 0.0, 0.0,
                    1.0, 0.0, 1.0,
                    0.0, 0.0, 0.0,
                    &metric);

  int speed_limited = 0;
  if(test_key == 10) {
    primitive_quantities prims;
    prims.vx = 0.0/0.0;
    prims.vy = 0.0/0.0;
    prims.vz = 0.0/0.0;
    limit_v_and_compute_u0(&hybrid_eos, &metric, &prims, &speed_limited);
  }

  const char tablepath[] = "SLy4_3335_rho391_temp163_ye66.h5";
  double Y_e_atm   = 0.5;
  double Y_e_min   = 0.05;
  double Y_e_max   = Y_e_atm;
  double T_atm     = 1e-2;
  double T_min     = T_atm;
  double T_max     = 1e2;

  if(test_key == 2) {
    rho_b_atm = -1;
  } else if(test_key == 3) {
    Y_e_atm = -1;
  } else if(test_key == 4) {
    T_atm = -1;
  } else if(test_key == 5) {
    rho_b_min = 1e301;
  } else if(test_key == 6) {
    Y_e_min = 1e1;
  } else if(test_key == 7) {
    T_min = 1e3;
  }


  eos_parameters tab_eos;
  //if( rho_min > rho_max ) grhayl_error("rho_min cannot be greater than rho_max\n");
  //if( Y_e_min > Y_e_max ) grhayl_error("Y_e_min cannot be greater than Y_e_max\n");
  //if(   T_min >   T_max ) grhayl_error("T_min cannot be greater than T_max\n");
  initialize_tabulated_eos_functions_and_params(tablepath, W_max,
                                                rho_b_atm, rho_b_min, rho_b_max,
                                                Y_e_atm, Y_e_min, Y_e_max,
                                                T_atm, T_min, T_max, &tab_eos);

  printf("We shouldn't be here, so I'll get rid of some compilation warnings :)\n"
         "%e %d\n", Fermi_Dirac_integral, speed_limited);
  return 0;
}
