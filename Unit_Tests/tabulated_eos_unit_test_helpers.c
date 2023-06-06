#include "nrpyeos_tabulated.h"

double get_table_quantity(
      const int which_var,
      const double rho,
      const double Y_e,
      const double T ) {
  switch (which_var) {
    case NRPyEOS_press_key:
      return rho + Y_e + T;
    case NRPyEOS_eps_key:
      return rho + Y_e - T;
    case NRPyEOS_entropy_key:
      return rho - Y_e + T;
    case NRPyEOS_munu_key:
      return rho - Y_e - T;
    case NRPyEOS_cs2_key:
      return -rho + Y_e + T;
    case NRPyEOS_depsdT_key:
      return -rho + Y_e - T;
    case NRPyEOS_dPdrho_key:
      return -rho - Y_e + T;
    case NRPyEOS_dPdeps_key:
      return -rho - Y_e - T;
    case NRPyEOS_muhat_key:
      return 2*rho + Y_e + T;
    case NRPyEOS_mu_e_key:
      return 2*rho + Y_e - T;
    case NRPyEOS_mu_p_key:
      return 2*rho - Y_e + T;
    case NRPyEOS_mu_n_key:
      return 2*rho - Y_e - T;
    case NRPyEOS_X_a_key:
      return -2*rho + Y_e + T;
    case NRPyEOS_X_h_key:
      return -2*rho + Y_e - T;
    case NRPyEOS_X_n_key:
      return -2*rho - Y_e + T;
    case NRPyEOS_X_p_key:
      return -2*rho - Y_e - T;
    case NRPyEOS_Abar_key:
      return rho + 3*Y_e - T;
    case NRPyEOS_Zbar_key:
      return rho - 3*Y_e + T;
    case NRPyEOS_Gamma_key:
      return rho - 3*Y_e - T;
    default:
      ghl_error("Invalid variable %d\n", which_var);
      return -1;
  }
}
