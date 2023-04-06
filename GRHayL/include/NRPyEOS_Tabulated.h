#ifndef NRPYEOS_TABULATED_H_
#define NRPYEOS_TABULATED_H_
#include "GRHayL.h"
#ifndef GRHAYL_USE_HDF5
#define HDF5_ERROR_IF_USED \
  grhayl_error("HDF5 is disabled, so this function cannot be used\n")
#else

#include <string.h>
#include <hdf5.h>
#define H5_USE_16_API 1

// Table keys
enum eos_keys {
  NRPyEOS_press_key, NRPyEOS_eps_key, NRPyEOS_entropy_key,
  NRPyEOS_munu_key, NRPyEOS_cs2_key, NRPyEOS_depsdT_key,
  NRPyEOS_dPdrho_key, NRPyEOS_dPdeps_key, NRPyEOS_muhat_key,
  NRPyEOS_mu_e_key, NRPyEOS_mu_p_key, NRPyEOS_mu_n_key,
  NRPyEOS_X_a_key, NRPyEOS_X_h_key, NRPyEOS_X_n_key,
  NRPyEOS_X_p_key, NRPyEOS_Abar_key, NRPyEOS_Zbar_key,
  NRPyEOS_Gamma_key, NRPyEOS_ntablekeys
};

// Unit conversion
#define CODE_TO_CGS_DENSITY  6.17714470405638e+17 // [Density]    = M_sun / L^3
#define CODE_TO_CGS_ENERGY   8.98755178736818e+20 // [Energy]     = c^2
#define CODE_TO_CGS_PRESSURE 5.55174079257738e+38 // [Pressure]   = M_sun / (L T^2)
#define CGS_TO_CODE_LENGTH   6.77269222552442e-06 // 1/[Length]   = c^2 / (G M_sun)
#define CGS_TO_CODE_TIME     2.03040204956746e+05 // 1/[Time]     = c / L
#define CGS_TO_CODE_DENSITY  1.61887093132742e-18 // 1/[Density]  = L^3 / M_sun
#define CGS_TO_CODE_PRESSURE 1.80123683248503e-39 // 1/[Pressure] = (L T^2) / M_sun
#define CGS_TO_CODE_ENERGY   1.11265005605362e-21 // 1/[Energy]   = 1 / c^2

// Name of the variables. This is only used to print
// information about the keys during startup
static const char table_var_names[NRPyEOS_ntablekeys][10] = {
  "logpress","logenergy","entropy","munu","cs2","dedt",
  "dpdrhoe", "dpderho", "muhat", "mu_e", "mu_p", "mu_n",
  "Xa","Xh","Xn","Xp","Abar","Zbar","Gamma"
};
//********************************************
#endif // GRHAYL_USE_HDF5

// Error handling struct
typedef struct NRPyEOS_error_report {
  bool error;
  int error_key;
  char message[1024];
} NRPyEOS_error_report;

#ifdef __cplusplus
extern "C" {
#endif

// Function prototypes
void NRPyEOS_read_table_set_EOS_params(
      const char *nuceos_table_name,
      eos_parameters *restrict eos_params);

void NRPyEOS_free_memory(eos_parameters *restrict eos_params);

void NRPyEOS_P_from_rho_Ye_T(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P);

void NRPyEOS_eps_from_rho_Ye_T(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict eps);

void NRPyEOS_P_and_eps_from_rho_Ye_T(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps);

void NRPyEOS_P_eps_and_S_from_rho_Ye_T(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict S);

void NRPyEOS_P_eps_and_cs2_from_rho_Ye_T(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict cs2);

void NRPyEOS_P_eps_S_and_cs2_from_rho_Ye_T(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict S,
      double *restrict cs2);

void NRPyEOS_P_eps_and_depsdT_from_rho_Ye_T(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict depsdT);

void NRPyEOS_P_eps_muhat_mue_mup_and_mun_from_rho_Ye_T(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict P,
      double *restrict eps,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n);

void NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p);

void NRPyEOS_P_and_T_from_rho_Ye_eps(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict P,
      double *restrict T);

void NRPyEOS_P_and_T_from_rho_Ye_S(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double S,
      double *restrict P,
      double *restrict T);

void NRPyEOS_P_S_depsdT_and_T_from_rho_Ye_eps(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double eps,
      double *restrict P,
      double *restrict S,
      double *restrict depsdT,
      double *restrict T);

void NRPyEOS_eps_S_and_T_from_rho_Ye_P(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double P,
      double *restrict eps,
      double *restrict S,
      double *restrict T);

void NRPyEOS_P_eps_and_T_from_rho_Ye_S(
      const eos_parameters *restrict eos_params,
      const double rho,
      const double Y_e,
      const double S,
      double *restrict P,
      double *restrict eps,
      double *restrict T);

void NRPyEOS_from_rho_Ye_T_interpolate_n_quantities(
      const eos_parameters *restrict eos_params,
      const int n,
      const double rho,
      const double Y_e,
      const double T,
      const int *restrict tablevars_keys,
      double *restrict tablevars,
      NRPyEOS_error_report *restrict report);

void NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities(
      const eos_parameters *restrict eos_params,
      const int n,
      const double prec,
      const double rho,
      const double Y_e,
      const double tablevar_in,
      const int tablevar_in_key,
      const int *restrict tablevars_keys,
      double *restrict tablevars,
      double *restrict T,
      NRPyEOS_error_report *restrict report);

void NRPyEOS_initialize_tabulated_functions(eos_parameters *restrict eos);

void NRPyEOS_tabulated_compute_enthalpy_and_cs2(
      eos_parameters const *restrict eos,
      primitive_quantities const *restrict prims,
      double *restrict enthalpy_ptr,
      double *restrict cs2_ptr );

#ifdef __cplusplus
}
#endif
#endif // NRPYEOS_TABULATED_H_
