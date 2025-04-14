#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "GRHayLib.h"

ghl_parameters *ghl_params;
ghl_eos_parameters *ghl_eos;

int parse_C2P_routine_keyword(const char *restrict routine_name);

void GRHayLib_paramcheck() {
  DECLARE_CCTK_PARAMETERS;

  if(rho_b_atm < 0)
    CCTK_ERROR("Parameter rho_b_atm must be set in the parameter file and be non-negative.");

  if(rho_b_min < 0 && rho_b_min != -1)
    CCTK_ERROR("Parameter rho_b_min must be non-negative.");

  if(rho_b_max < 0 && rho_b_max != -1)
    CCTK_ERROR("Parameter rho_b_max must be non-negative.");

  if( CCTK_EQUALS(EOS_type, "Simple") ) {
    if(Gamma < 0)
      CCTK_ERROR("Parameter Gamma must be set in the parameter file and be non-negative.");

    if(P_atm < 0)
      CCTK_ERROR("Parameter P_atm must be set in the parameter file and be non-negative.");
  
    if(P_min < 0 && P_min != -1)
      CCTK_ERROR("Parameter P_min must be non-negative.");
  
    if(P_max < 0 && P_max != -1)
      CCTK_ERROR("Parameter P_max must be non-negative.");

    if( CCTK_EQUALS(con2prim_routine, "Newman1D_energy")  ||
        CCTK_EQUALS(con2prim_routine, "Font1D") ||
        CCTK_EQUALS(con2prim_routine, "Newman1D_entropy") ) {
      CCTK_VERROR("Selected parameter option for con2prim_routine %s is incompatible with\n"
                  "simple EOS. Please change the routine.", con2prim_routine);
    }
    for(int i=0; i<3; i++) {
      if( CCTK_EQUALS(con2prim_backup_routines[i], "Newman1D_energy")  ||
          CCTK_EQUALS(con2prim_backup_routines[i], "Font1D") ||
          CCTK_EQUALS(con2prim_backup_routines[i], "Newman1D_entropy") ) {
         CCTK_VERROR("Selected parameter option for con2prim_backup_routines[%d] = %s is incompatible with\n"
                     "simple EOS. Please change the routine.", i, con2prim_backup_routines[0]);
      }
    }
  } else if( CCTK_EQUALS(EOS_type, "Hybrid") ) {
    if(Gamma_th < 0)
      CCTK_ERROR("Parameter Gamma_th must be set in the parameter file and be non-negative.");

    if(Gamma_ppoly_in[0] < 0)
      CCTK_ERROR("Parameter Gamma_ppoly_in[0] must be set in the parameter file and be non-negative.");

    if(k_ppoly0 < 0)
      CCTK_ERROR("Parameter k_ppoly0 must be set in the parameter file and be non-negative.");

    if(neos > 1) {
      if(rho_ppoly_in[0] < 0)
        CCTK_ERROR("Parameter rho_ppoly_in[0] must be set in the parameter file and be non-negative.");

      for(int i=1; i<neos-1; i++) {
        if(Gamma_ppoly_in[i] < 0)
          CCTK_VERROR("Parameter Gamma_ppoly_in[%d] must be set in the parameter file and be non-negative.", i);

        if(rho_ppoly_in[i] < 0)
          CCTK_VERROR("Parameter rho_ppoly_in[%d] must be set in the parameter file and be non-negative.", i);
      }
      if(Gamma_ppoly_in[neos-1] < 0)
        CCTK_VERROR("Parameter Gamma_ppoly_in[%d] must be set in the parameter file and be non-negative.", neos-1);
    }

    if( CCTK_EQUALS(con2prim_routine, "Newman1D_energy")  ||
        CCTK_EQUALS(con2prim_routine, "Newman1D_entropy") ) {
      CCTK_VERROR("Selected parameter option for con2prim_routine %s is incompatible with\n"
                  "hybrid EOS. Please change the routine.", con2prim_routine);
    }
    for(int i=0; i<3; i++) {
      if( CCTK_EQUALS(con2prim_backup_routines[i], "Newman1D_energy")  ||
          CCTK_EQUALS(con2prim_backup_routines[i], "Newman1D_entropy") ) {
         CCTK_VERROR("Selected parameter option for con2prim_backup_routines[%d] = %s is incompatible with\n"
                     "hybrid EOS. Please change the routine.", i, con2prim_backup_routines[0]);
      }
    }
  } else if( CCTK_EQUALS(EOS_type, "Tabulated") ) {
    if(Y_e_atm == -1)
      CCTK_ERROR("Parameter Y_e_atm must be set in the parameter file");
    else if(Y_e_atm < 0)
      CCTK_ERROR("Parameter Y_e_atm must be non-negative.");

    if(Y_e_min < 0 && Y_e_min != -1)
      CCTK_ERROR("Parameter Y_e_min must be non-negative.");

    if(Y_e_max < 0 && Y_e_max != -1)
      CCTK_ERROR("Parameter Y_e_max must be non-negative.");

    if(T_atm == -1)
      CCTK_ERROR("Parameter T_atm must be set in the parameter file");
    else if(T_atm < 0)
      CCTK_ERROR("Parameter T_atm must be non-negative.");

    if(T_min < 0 && T_min != -1)
      CCTK_ERROR("Parameter T_min must be non-negative.");

    if(T_max < 0 && T_max != -1)
      CCTK_ERROR("Parameter T_max must be non-negative.");

    if( CCTK_EQUALS(con2prim_routine, "Noble2D")  ||
        CCTK_EQUALS(con2prim_routine, "Noble1D") ||
        CCTK_EQUALS(con2prim_routine, "Noble1D_entropy") ||
        CCTK_EQUALS(con2prim_routine, "Font1D") ) {
      CCTK_VERROR("Selected parameter option for con2prim_routine %s is incompatible with\n"
                  "tabulated EOS. Please change the routine.", con2prim_routine);
    }
    for(int i=0; i<3; i++) {
      if( CCTK_EQUALS(con2prim_backup_routines[i], "Noble2D")  ||
          CCTK_EQUALS(con2prim_backup_routines[i], "Noble1D") ||
          CCTK_EQUALS(con2prim_backup_routines[i], "Noble1D_entropy") ||
          CCTK_EQUALS(con2prim_backup_routines[i], "Font1D") ) {
         CCTK_VERROR("Selected parameter option for con2prim_backup_routines[%d] = %s is incompatible with\n"
                     "tabulated EOS. Please change the routine.", i, con2prim_backup_routines[0]);
      }
    }
  }

  if(!evolve_entropy) {
    if( CCTK_EQUALS(con2prim_routine, "Noble1D_entropy")  ||
        CCTK_EQUALS(con2prim_routine, "Newman1D_entropy") ||
        CCTK_EQUALS(con2prim_routine, "Palenzuela1D_entropy") ) {
      CCTK_VERROR("Selected parameter option for con2prim_routine %s requires entropy.\n"
                  "Please set evolve_entropy=\"yes\" or change the routine.", con2prim_routine);
    }
    for(int i=0; i<3; i++) {
      if( CCTK_EQUALS(con2prim_backup_routines[i], "Noble1D_entropy")  ||
          CCTK_EQUALS(con2prim_backup_routines[i], "Newman1D_entropy") ||
          CCTK_EQUALS(con2prim_backup_routines[i], "Palenzuela1D_entropy") ) {
         CCTK_VERROR("Selected parameter option for con2prim_backup_routines[%d] = %s requires entropy.\n"
                     "Please set evolve_entropy=\"yes\" or change the routine.", i, con2prim_backup_routines[0]);
      }
    }
  }
}

void GRHayLib_initialize(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;

  ghl_params = (ghl_parameters *)malloc(sizeof(ghl_parameters));
  ghl_eos = (ghl_eos_parameters *)malloc(sizeof(ghl_eos_parameters));

  const int main = parse_C2P_routine_keyword(con2prim_routine);
  int backups[3];
  for(int i=0; i<3; i++)
    backups[i] = parse_C2P_routine_keyword(con2prim_backup_routines[i]);

  GRHayLib_paramcheck();

  ghl_initialize_params(
      main, backups,
      evolve_entropy, evolve_temperature,
      calc_primitive_guess, Psi6threshold,
      max_Lorentz_factor, Lorenz_damping_factor,
      ghl_params);

      ghl_params->ppm_flattening_epsilon = ppm_flattening_epsilon;
      ghl_params->ppm_flattening_omega1  = ppm_flattening_omega1;
      ghl_params->ppm_flattening_omega2  = ppm_flattening_omega2;

      ghl_params->ppm_shock_epsilon = ppm_shock_epsilon;
      ghl_params->ppm_shock_eta1    = ppm_shock_eta1;
      ghl_params->ppm_shock_eta2    = ppm_shock_eta2;
      ghl_params->ppm_shock_k0      = ppm_shock_k0;

  if (CCTK_EQUALS(EOS_type, "Simple")) {
    if(main == Font1D || backups[0] == Font1D || backups[1] == Font1D || backups[2] == Font1D)
      CCTK_VERROR("Error: Font1D routine is incompatible with ideal fluid EOS. Please choose a different Con2Prim routine.");
    ghl_con2prim_multi_method = ghl_con2prim_hybrid_multi_method;
    ghl_initialize_simple_eos_functions_and_params(
          rho_b_atm, rho_b_min, rho_b_max,
          P_atm, P_min, P_max,
          Gamma, ghl_eos);
  } else if (CCTK_EQUALS(EOS_type, "Hybrid")) {
    ghl_con2prim_multi_method = ghl_con2prim_hybrid_multi_method;
    ghl_initialize_hybrid_eos_functions_and_params(
          rho_b_atm, rho_b_min, rho_b_max,
          neos, rho_ppoly_in,
          Gamma_ppoly_in, k_ppoly0,
          Gamma_th, ghl_eos);
  } else if (CCTK_EQUALS(EOS_type, "Tabulated")) {
    if( CCTK_EQUALS(EOS_tablepath, "") )
      CCTK_ERROR("Parameter EOS_tablepath uninitialized.");

    ghl_con2prim_multi_method = ghl_con2prim_tabulated_multi_method;
    ghl_initialize_tabulated_eos_functions_and_params(
          EOS_tablepath,
          rho_b_atm, rho_b_min, rho_b_max,
          Y_e_atm, Y_e_min, Y_e_max,
          T_atm, T_min, T_max,
          ghl_eos);
  } else if (CCTK_EQUALS(EOS_type, "")) {
    CCTK_VERROR("GRHayLib parameter EOS_type is unset. Please set an EOS type.");
  } else {
    CCTK_VERROR("GRHayLib parameter EOS_type has an unsupported type. Please check"
                " the list of parameter options in the param.ccl.");
  }
}

void GRHayLib_terminate(CCTK_ARGUMENTS) {
  if(ghl_eos->eos_type == ghl_eos_tabulated) {
    free(ghl_eos->table_all);
    free(ghl_eos->table_logrho);
    free(ghl_eos->table_logT);
    free(ghl_eos->table_Y_e);
    free(ghl_eos->table_eps);
  }
  free(ghl_eos);
  free(ghl_params);
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
  } else if (CCTK_EQUALS(routine_name, "Font1D")) {
    return Font1D;
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
  } else if (CCTK_EQUALS(routine_name, "Newman1D_entropy")) {
    return Newman1D_entropy;
  }
  return -100;
}
