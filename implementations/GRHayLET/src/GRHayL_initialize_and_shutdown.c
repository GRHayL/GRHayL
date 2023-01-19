#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLET.h"

GRHayL_parameters *grhayl_params;
eos_parameters *grhayl_eos;

int parse_C2P_routine_keyword(const char *restrict routine_name);

void GRHayLET_initialize(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;

  grhayl_params = (GRHayL_parameters *)malloc(sizeof(GRHayL_parameters));
  grhayl_eos = (eos_parameters *)malloc(sizeof(eos_parameters));

  const int main = parse_C2P_routine_keyword(con2prim_routine);
  int backups[3];
  for(int i=0; i<3; i++)
    backups[i] = parse_C2P_routine_keyword(con2prim_backup_routines[i]);

  initialize_GRHayL(main, backups,
                    evolve_entropy, evolve_temperature,
                    calc_primitive_guess, Psi6threshold,
                    update_Tmunu, Cupp_Fix, grhayl_params);

  int grhayl_eos_type;
  if (CCTK_EQUALS(EOS_type, "hybrid")) {
    grhayl_eos_type = 0;
  } else if (CCTK_EQUALS(EOS_type, "tabulated")) {
    grhayl_eos_type = 1;
  } else if (CCTK_EQUALS(EOS_type, "")) {
    CCTK_VERROR("GRHayLET parameter EOS_type is unset. Please set an EOS type.");
  } else {
    CCTK_VERROR("GRHayLET parameter EOS_type has an unsupported type. Please check"
                " the list of parameter options in the param.ccl.");
  }

  initialize_general_eos(grhayl_eos_type, gamma_speed_limit,
             rho_b_atm, rho_b_atm, rho_b_max,
             grhayl_eos);

  if(grhayl_eos->eos_type == 0) {
    initialize_hybrid_eos(
               neos, rho_ppoly_in,
               gamma_ppoly_in, k_ppoly0,
               gamma_th, grhayl_eos);
  } else if(grhayl_eos->eos_type == 1) {
double Ye_max = 100; //Temporary until the parameter is set up
    initialize_tabulated_functions(grhayl_eos);
    grhayl_eos->tabulated_read_table_set_EOS_params(EOS_tablepath, grhayl_eos);
    initialize_tabulated_eos(
               EOS_root_finding_precision,
               depsdT_threshold,
               Ye_atm, Ye_atm, Ye_max,
               T_atm, T_atm, T_max,
               grhayl_eos);
  }
}

void GRHayLET_terminate(CCTK_ARGUMENTS) {
  if(grhayl_eos->eos_type == 1) {
    free(grhayl_eos->table_all);
    free(grhayl_eos->table_logrho);
    free(grhayl_eos->table_logT);
    free(grhayl_eos->table_Ye);
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

