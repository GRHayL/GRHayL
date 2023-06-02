#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLib.h"

grhayl_parameters *grhayl_params;
eos_parameters *grhayl_eos;

int parse_C2P_routine_keyword(const char *restrict routine_name);

typedef struct {
  double rho_atm, rho_min, rho_max;
  double Y_e_atm, Y_e_min, Y_e_max;
  double T_atm  , T_min  , T_max;
} grhayl_params_checked;

void paramcheck() {

  DECLARE_CCTK_PARAMETERS;

  if( rho_b_atm == -1 )
    CCTK_ERROR("Parameter rho_b_atm must be set in the parameter file");
  else if( rho_b_atm < 0 )
    CCTK_ERROR("Parameter rho_b_atm must be non-negative");

  if( rho_b_min < 0 && rho_b_min != -1 )
    CCTK_ERROR("Parameter rho_b_min must be non-negative");

  if( rho_b_max < 0 && rho_b_max != -1 )
    CCTK_ERROR("Parameter rho_b_max must be non-negative");

  if( CCTK_EQUALS(EOS_type, "tabulated") ) {
    if( Y_e_atm == -1 )
      CCTK_ERROR("Parameter Y_e_atm must be set in the parameter file");
    else if( Y_e_atm < 0 )
      CCTK_ERROR("Parameter Y_e_atm must be non-negative");

    if( Y_e_min < 0 && Y_e_min != -1 )
      CCTK_ERROR("Parameter Y_e_min must be non-negative");

    if( Y_e_max < 0 && Y_e_max != -1 )
      CCTK_ERROR("Parameter Y_e_max must be non-negative");

    if( T_atm == -1 )
      CCTK_ERROR("Parameter T_atm must be set in the parameter file");
    else if( T_atm < 0 )
      CCTK_ERROR("Parameter T_atm must be non-negative");

    if( T_min < 0 && T_min != -1 )
      CCTK_ERROR("Parameter T_min must be non-negative");

    if( T_max < 0 && T_max != -1 )
      CCTK_ERROR("Parameter T_max must be non-negative");
  }
}

void GRHayLib_initialize(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;

  grhayl_params = (grhayl_parameters *)malloc(sizeof(grhayl_parameters));
  grhayl_eos = (eos_parameters *)malloc(sizeof(eos_parameters));

  const int main = parse_C2P_routine_keyword(con2prim_routine);
  int backups[3];
  for(int i=0; i<3; i++)
    backups[i] = parse_C2P_routine_keyword(con2prim_backup_routines[i]);

  paramcheck();

  CCTK_VINFO("%.15e %.15e %.15e", rho_b_atm, rho_b_min, rho_b_max);
  CCTK_VINFO("%.15e %.15e %.15e", Y_e_atm, Y_e_min, Y_e_max);
  CCTK_VINFO("%.15e %.15e %.15e", T_atm, T_min, T_max);

  grhayl_initialize_params(
      main, backups,
      evolve_entropy, evolve_temperature,
      calc_primitive_guess, Psi6threshold,
      Cupp_Fix, Lorenz_damping_factor,
      grhayl_params);

  if (CCTK_EQUALS(EOS_type, "hybrid")) {

    grhayl_initialize_hybrid_eos_functions_and_params(
          W_max,
          rho_b_atm, rho_b_min, rho_b_max,
          neos, rho_ppoly_in,
          Gamma_ppoly_in, k_ppoly0,
          Gamma_th, grhayl_eos);
  } else if (CCTK_EQUALS(EOS_type, "tabulated")) {
    if( CCTK_EQUALS(EOS_tablepath, "") )
      CCTK_ERROR("Parameter EOS_tablepath uninitialized.");

    grhayl_initialize_tabulated_eos_functions_and_params(
          EOS_tablepath, W_max,
          rho_b_atm, rho_b_min, rho_b_max,
          Y_e_atm, Y_e_min, Y_e_max,
          T_atm, T_min, T_max,
          grhayl_eos);
  } else if (CCTK_EQUALS(EOS_type, "")) {
    CCTK_VERROR("GRHayLib parameter EOS_type is unset. Please set an EOS type.");
  } else {
    CCTK_VERROR("GRHayLib parameter EOS_type has an unsupported type. Please check"
                " the list of parameter options in the param.ccl.");
  }
}

void GRHayLib_terminate(CCTK_ARGUMENTS) {
  if(grhayl_eos->eos_type == grhayl_eos_tabulated) {
    free(grhayl_eos->table_all);
    free(grhayl_eos->table_logrho);
    free(grhayl_eos->table_logT);
    free(grhayl_eos->table_Y_e);
    free(grhayl_eos->table_eps);
  }
  free(grhayl_eos);
  free(grhayl_params);
}

int parse_C2P_routine_keyword(const char *restrict routine_name) {
  if (CCTK_EQUALS(routine_name, "None")) {
    return None;
  } else if (CCTK_EQUALS(routine_name, "Noble2D")) {
    return Noble2D;
  } else if (CCTK_EQUALS(routine_name, "Noble1D")) {
    return Noble1D;
  } else if (CCTK_EQUALS(routine_name, "Noble1D_entropy")) {
    return Noble1D_entropy;
  } else if (CCTK_EQUALS(routine_name, "Noble1D_entropy2")) {
    return Noble1D_entropy2;
  } else if (CCTK_EQUALS(routine_name, "FontFix")) {
    return FontFix;
  } else if (CCTK_EQUALS(routine_name, "CerdaDuran2D")) {
    return CerdaDuran2D;
  } else if (CCTK_EQUALS(routine_name, "CerdaDuran3D")) {
    return CerdaDuran3D;
  } else if (CCTK_EQUALS(routine_name, "Palenzuela1D")) {
    return Palenzuela1D;
  } else if (CCTK_EQUALS(routine_name, "Palenzuela1D_entropy")) {
    return Palenzuela1D_entropy;
  } else if (CCTK_EQUALS(routine_name, "Newman1D")) {
    return Newman1D;
  }
  return -100;
}
