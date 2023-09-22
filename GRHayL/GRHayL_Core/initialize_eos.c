#include "con2prim.h"
#include "nrpyeos_hybrid.h"
#include "nrpyeos_tabulated.h"
#include "ghl_eos_functions_declaration.h"

#define init_common_eos_quantities         \
  eos->rho_atm           = rho_atm;        \
  eos->rho_min           = rho_min;        \
  eos->rho_max           = rho_max;

int (*ghl_con2prim_multi_method)(
      const ghl_parameters *restrict params,
      const ghl_eos_parameters *restrict eos,
      const ghl_metric_quantities *restrict ADM_metric,
      const ghl_ADM_aux_quantities *restrict metric_aux,
      const ghl_conservative_quantities *restrict cons,
      ghl_primitive_quantities *restrict prim,
      ghl_con2prim_diagnostics *restrict diagnostics);

/*
 * Function    : ghl_initialize_eos_functions()
 * Description : Initializes function pointers in EOS struct to NRPyEOS
 *
 * Input/Output: eos - ghl_eos_parameters struct with the function pointers
 *                     initialized
 */
void ghl_initialize_eos_functions(
      ghl_eos_t const eos_type,
      ghl_eos_parameters *restrict eos) {

  // Step 1: Hybrid EOS functions (always available)
  NRPyEOS_initialize_hybrid_functions(eos);

  // Step 2: Tabulated EOS functions (always available)
  NRPyEOS_initialize_tabulated_functions(eos);

  // Step 3: General functions (same interface for all EOSs)
  if( eos_type == ghl_eos_hybrid || eos_type == ghl_eos_simple) {
    ghl_con2prim_multi_method = ghl_con2prim_hybrid_multi_method;
    ghl_compute_h_and_cs2 = NRPyEOS_hybrid_compute_enthalpy_and_cs2;
  }
  else if( eos_type == ghl_eos_tabulated ) {
    ghl_con2prim_multi_method = ghl_con2prim_tabulated_multi_method;
    ghl_compute_h_and_cs2 = NRPyEOS_tabulated_compute_enthalpy_and_cs2;
  }
}

/*
 * Function    : ghl_initialize_simple_eos()
 * Description : Initializes EOS struct elements for a simple EOS
*/
void ghl_initialize_simple_eos(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double press_atm,
      double press_min,
      double press_max,
      const double Gamma,
      ghl_eos_parameters *restrict eos ) {

  // Step 0: Enforce default values
  if( rho_atm < 0 ) ghl_error("rho_atm must be specified\n");
  if( rho_min < 0 ) {
    ghl_warn("Minimum density not provided. Disabling density floor (rho_min = 0)\n");
    rho_min = 0.0;
  }
  if( rho_max < 0 ) {
    ghl_warn("Maximum density not provided. Disabling density ceiling (rho_max = 1e300)\n");
    rho_max = 1e300;
  }
  if( rho_min > rho_max ) ghl_error("rho_min cannot be greater than rho_max\n");

  if( press_atm < 0 ) ghl_error("press_atm must be specified\n");
  if( press_min < 0 ) {
    ghl_warn("Minimum pressure not provided. Disabling pressure floor (press_min = 0)\n");
    press_min = 0.0;
  }
  if( press_max < 0 ) {
    ghl_warn("Maximum pressure not provided. Disabling pressure ceiling (press_max = 1e300)\n");
    press_max = 1e300;
  }
  if( press_min > press_max ) ghl_error("press_min cannot be greater than press_max\n");

  // Step 1: Set EOS type to Ideal Fluid
  eos->eos_type = ghl_eos_simple;

  // Step 2: Initialize quantities which are common to all EOSs.
  init_common_eos_quantities;

  // Step 3: Set basic ideal fluid EOS parameters.
  eos->press_atm = press_atm;
  eos->press_min = press_min;
  eos->press_max = press_max;
  eos->Gamma_th = eos->Gamma_ppoly[0] = Gamma;
  eos->neos = 1;
  eos->K_ppoly[0] = 1;
  eos->rho_ppoly[0] = 0.0;
  eos->eps_integ_const[0] = 0.0;

  const double Gm1 = Gamma - 1.0;
  // -------------- Ceilings --------------
  eos->eps_max = eos->press_max/(eos->rho_max*Gm1);
  eos->entropy_max = ghl_hybrid_compute_entropy_function(eos, eos->rho_max, eos->press_max);

  // --------------- Floors ---------------
  eos->eps_min = eos->press_min/(eos->rho_min*Gm1);
  eos->entropy_min = ghl_hybrid_compute_entropy_function(eos, eos->rho_min, eos->press_min);

  // --------- Atmospheric values ---------
  eos->eps_atm = eos->press_atm/(eos->rho_atm*Gm1);
  eos->entropy_atm = ghl_hybrid_compute_entropy_function(eos, eos->rho_atm, eos->press_atm);

  // Compute atmospheric tau
  eos->tau_atm = eos->rho_atm * eos->eps_atm;
  // --------------------------------------
}

/*
 * Function    : ghl_initialize_hybrid_eos()
 * Description : Initializes EOS struct elements for a hybrid EOS
*/
void ghl_initialize_hybrid_eos(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      ghl_eos_parameters *restrict eos ) {

  // Step 0: Enforce default values
  if( rho_atm < 0 ) ghl_error("rho_atm must be specified\n");
  if( rho_min < 0 ) {
    ghl_warn("Minimum density not provided. Disabling density floor (rho_min = 0)\n");
    rho_min = 0.0;
  }
  if( rho_max < 0 ) {
    ghl_warn("Maximum density not provided. Disabling density ceiling (rho_max = 1e300)\n");
    rho_max = 1e300;
  }
  if( rho_min > rho_max ) ghl_error("rho_min cannot be greater than rho_max\n");

  // Step 1: Set EOS type to Hybrid
  eos->eos_type = ghl_eos_hybrid;

  // Step 2: Initialize quantities which are common to all EOSs.
  init_common_eos_quantities;

  // Step 3: Set basic Hybrid EOS parameters.
  eos->neos = neos;
  eos->Gamma_th = Gamma_th;
  eos->K_ppoly[0] = K_ppoly0;
  if( neos==1 ) {
    eos->rho_ppoly[0] = rho_ppoly[0];
    eos->eps_integ_const[0] = 0.0;
  } else {
    for(int j=0; j<=neos-2; j++) eos->rho_ppoly[j] = rho_ppoly[j];
  }
  for(int j=0; j<=neos-1; j++) eos->Gamma_ppoly[j] = Gamma_ppoly[j];

  // Step 4: Initialize {K_{j}}, j>=1, and {eps_integ_const_{j}}
  ghl_hybrid_set_K_ppoly_and_eps_integ_consts(eos);

  // -------------- Ceilings --------------
  // Compute maximum P and eps
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, eos->rho_max, &eos->press_max, &eos->eps_max);

  // Compute maximum entropy
  eos->entropy_max = ghl_hybrid_compute_entropy_function(eos, eos->rho_max, eos->press_max);
  // --------------------------------------

  // --------------- Floors ---------------
  // Compute maximum P and eps
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, eos->rho_min, &eos->press_min, &eos->eps_min);

  // Compute maximum entropy
  eos->entropy_min = ghl_hybrid_compute_entropy_function(eos, eos->rho_min, eos->press_min);
  // --------------------------------------

  // --------- Atmospheric values ---------
  // Compute atmospheric P and eps
  ghl_hybrid_compute_P_cold_and_eps_cold(eos, eos->rho_atm, &eos->press_atm, &eos->eps_atm);

  // Compute atmospheric entropy
  eos->entropy_atm = ghl_hybrid_compute_entropy_function(eos, eos->rho_atm, eos->press_atm);

  // Compute atmospheric tau
  eos->tau_atm = eos->rho_atm * eos->eps_atm;
  // --------------------------------------
}

/*
 * Function    : ghl_initialize_tabulated_eos()
 * Description : Initializes EOS struct elements for tabulated EOS
*/
void ghl_initialize_tabulated_eos(
      const char *table_filepath,
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double Y_e_atm,
      double Y_e_min,
      double Y_e_max,
      const double T_atm,
      double T_min,
      double T_max,
      ghl_eos_parameters *restrict eos ) {

  // Step 1: Set EOS type to Tabulated.
  eos->eos_type = ghl_eos_tabulated;

  // Step 2: Read the EOS table
  ghl_tabulated_read_table_set_EOS_params(table_filepath, eos);

  // Step 3: Enforce default values for (rho, Y_e, T) min, max, and atm
  // Step 3.a: Atmosphere values
  if( rho_atm < 0 ) ghl_error("rho_atm must be specified\n");
  if( Y_e_atm < 0 ) ghl_error("Y_e_atm must be specified\n");
  if(   T_atm < 0 ) ghl_error("T_atm must be specified\n");

  // Step 3.b: Minimum values
  if( rho_min < 0 ) {
    ghl_warn("Minimum density not provided; using table bounds (%.15e)\n", eos->table_rho_min);
    rho_min = eos->table_rho_min;
  }
  if( Y_e_min < 0 ) {
    ghl_warn("Minimum electron fraction not provided; using table bounds (%.15e)\n", eos->table_Y_e_min);
    Y_e_min = eos->table_Y_e_min;
  }
  if( T_min < 0 ) {
    ghl_warn("Minimum temperature not provided; using table bounds (%.15e)\n", eos->table_T_min);
    T_min = eos->table_T_min;
  }

  // Step 3.c: Maximum values
  if( rho_max < 0 ) {
    ghl_warn("Maximum density not provided; using table bounds (%.15e)\n", eos->table_rho_max);
    rho_max = eos->table_rho_max;
  }
  if( Y_e_max < 0 ) {
    ghl_warn("Maximum electron fraction not provided; using table bounds (%.15e)\n", eos->table_Y_e_max);
    Y_e_max = eos->table_Y_e_max;
  }
  if( T_max < 0 ) {
    ghl_warn("Maximum temperature not provided; using table bounds (%.15e)\n", eos->table_T_max);
    T_max = eos->table_T_max;
  }

  // Step 3.d: Sanity check mins and maxs
  if( rho_min > rho_max ) ghl_error("rho_min cannot be greater than rho_max\n");
  if( Y_e_min > Y_e_max ) ghl_error("Y_e_min cannot be greater than Y_e_max\n");
  if(   T_min >   T_max ) ghl_error("T_min cannot be greater than T_max\n");

  // Step 4: Initialize quantities which are common to all EOSs.
  init_common_eos_quantities;

  // Step 5: Set parameters specific to Tabulated EOS.
  eos->Y_e_atm = Y_e_atm;
  eos->Y_e_min = Y_e_min;
  eos->Y_e_max = Y_e_max;
  eos->T_atm   = T_atm;
  eos->T_min   = T_min;
  eos->T_max   = T_max;
  ghl_tabulated_compute_P_eps_S_from_T(eos,
                                        eos->rho_atm,
                                        Y_e_atm, T_atm,
                                        &eos->press_atm,
                                        &eos->eps_atm,
                                        &eos->entropy_atm);

  // Step 6: These parameters are manually set here, but
  //         can be overwritten later.
  eos->root_finding_precision = 1e-10;

  // Step 7: Set minimum values for eps, P, and S
  eos->press_min   = eos->table_P_min;
  eos->press_max   = eos->table_P_max;
  eos->eps_min     = eos->table_eps_min;
  eos->eps_max     = eos->table_eps_max;
  eos->entropy_min = eos->table_ent_min;
  eos->entropy_max = eos->table_ent_max;
  eos->tau_atm     = eos->rho_min * eos->eps_min;
}

/*
 * Function    : ghl_initialize_hybrid_eos_functions_and_params()
 * Description : Fully initializes EOS struct elements for a hybrid EOS
*/
void ghl_initialize_simple_eos_functions_and_params(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const double press_atm,
      double press_min,
      double press_max,
      const double Gamma,
      ghl_eos_parameters *restrict eos ) {

  // Step 1: Initialize Hybrid EOS functions
  ghl_initialize_eos_functions(ghl_eos_simple, eos);

  // Step 2: Initialize Hybrid EOS parameters
  ghl_initialize_simple_eos(
        rho_atm, rho_min, rho_max,
        press_atm, press_min, press_max,
        Gamma, eos);
}

/*
 * Function    : ghl_initialize_hybrid_eos_functions_and_params()
 * Description : Fully initializes EOS struct elements for a hybrid EOS
*/
void ghl_initialize_hybrid_eos_functions_and_params(
      const double rho_atm,
      double rho_min,
      double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      ghl_eos_parameters *restrict eos ) {

  // Step 1: Initialize Hybrid EOS functions
  ghl_initialize_eos_functions(ghl_eos_hybrid, eos);

  // Step 2: Initialize Hybrid EOS parameters
  ghl_initialize_hybrid_eos(
        rho_atm, rho_min, rho_max,
        neos, rho_ppoly, Gamma_ppoly,
        K_ppoly0, Gamma_th, eos);
}

/* Function    : ghl_initialize_tabulated_eos()
 * Description : Initializes EOS struct elements for tabulated EOS
*/
void ghl_initialize_tabulated_eos_functions_and_params(
      const char *table_filepath,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      const double T_atm,
      const double T_min,
      const double T_max,
      ghl_eos_parameters *restrict eos ) {

  eos->eos_type = ghl_eos_tabulated;

  // Step 1: Initialize Tabulated EOS functions
  ghl_initialize_eos_functions(ghl_eos_tabulated, eos);

  // Step 2: Initialize Tabulated EOS parameters
  ghl_initialize_tabulated_eos(table_filepath,
                           rho_atm, rho_min, rho_max,
                           Ye_atm, Ye_min, Ye_max,
                           T_atm, T_min, T_max,
                           eos);
}
