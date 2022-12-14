#ifndef TABULATED_EOS_H_
#define TABULATED_EOS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "../../GRHayL_Core/GRHayL.h"
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

// Error handling struct
typedef struct NRPyEOS_error_report {
  bool error;
  int error_key;
  char message[512];
} NRPyEOS_error_report;
//********************************************

#include "NRPy_function_prototypes.h"

#endif // TABULATED_EOS_H_
