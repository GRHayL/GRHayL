#include "NRPyEOS_Hybrid.h"
#include "NRPyEOS_Tabulated.h"

#define init_common_eos_quantities         \
  eos->W_max             = W_max;          \
  eos->inv_W_max_squared = 1.0/SQR(W_max); \
  eos->rho_atm           = rho_atm;        \
  eos->rho_min           = rho_min;        \
  eos->rho_max           = rho_max;

/*
 * Function    : initialize_eos_functions()
 * Description : Initializes function pointers in EOS struct to NRPyEOS
 *
 * Input/Output: eos - eos_parameters struct with the function pointers
 *                     initialized
 */
void initialize_eos_functions(
      grhayl_eos_t const eos_type,
      eos_parameters *restrict eos) {

  // Step 1: Hybrid EOS functions (always available)
  NRPyEOS_initialize_hybrid_functions(eos);

  // Step 2: Tabulated EOS functions (always available)
  NRPyEOS_initialize_tabulated_functions(eos);

  // Step 3: General functions (same interface for all EOSs)
  if( eos_type == grhayl_eos_hybrid ) {
    eos->compute_h_and_cs2 = &NRPyEOS_hybrid_compute_enthalpy_and_cs2;
  }
  else if( eos_type == grhayl_eos_tabulated ) {
    eos->compute_h_and_cs2 = &NRPyEOS_tabulated_compute_enthalpy_and_cs2;
  }
}

/*
 * Function    : initialize_hybrid_eos()
 * Description : Initializes EOS struct elements for a hybrid EOS
 *
 * Inputs      : W_max          - Maximum allowed Lorentz factor
 *             : rho_atm        - atmospheric value for rho
 *             : rho_min        - minimum allowable value for rho
 *             : rho_max        - maximum allowable value for rho
 *             : neos           - number of pieces in the piecewise
 *                                polytrope
 *             : rho_ppoly      - pointer to the array containing the
 *                                minimum rho_b for each polytropic piece
 *             : Gamma_ppoly    - pointer to the array containing the
 *                                minimum rho_b for each polytropic piece
 *             : K_ppoly0       - TODO: comment
 *             : Gamma_th       - TODO: comment
 *
 * Outputs     : eos            - eos_parameters struct with the above inputs
 *                                initialized
 */
void initialize_hybrid_eos(
      const double W_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      eos_parameters *restrict eos ) {

  // Step 1: Set EOS type to Hybrid
  eos->eos_type = grhayl_eos_hybrid;

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
  eos->hybrid_set_K_ppoly_and_eps_integ_consts(eos);

  // -------------- Ceilings --------------
  // Compute maximum P and eps
  eos->hybrid_compute_P_cold_and_eps_cold(eos, eos->rho_max, &eos->press_max, &eos->eps_max);

  // Compute maximum entropy
  eos->hybrid_compute_entropy_function(eos, eos->rho_max, eos->press_max, &eos->entropy_max);
  // --------------------------------------

  // --------------- Floors ---------------
  // Compute maximum P and eps
  eos->hybrid_compute_P_cold_and_eps_cold(eos, eos->rho_min, &eos->press_min, &eos->eps_min);

  // Compute maximum entropy
  eos->hybrid_compute_entropy_function(eos, eos->rho_min, eos->press_min, &eos->entropy_min);
  // --------------------------------------

  // --------- Atmospheric values ---------
  // Compute atmospheric P and eps
  eos->hybrid_compute_P_cold_and_eps_cold(eos, eos->rho_atm, &eos->press_atm, &eos->eps_atm);

  // Compute atmospheric entropy
  eos->hybrid_compute_entropy_function(eos, eos->rho_atm, eos->press_atm, &eos->entropy_atm);

  // Compute atmospheric tau
  eos->tau_atm = eos->rho_atm * eos->eps_atm;
  // --------------------------------------
}

/*
 * Function    : initialize_tabulated_eos()
 * Description : Initializes EOS struct elements for tabulated EOS
 *
 * Inputs      : W_max          - maximum allowed Lorentz factor
 *             : rho_atm        - atmospheric value for rho
 *             : rho_min        - minimum allowable value for rho
 *             : rho_max        - maximum allowable value for rho
 *             : Y_e_atm        - atmospheric value for Y_e
 *             : Y_e_min        - minimum allowable value for Y_e
 *             : Y_e_max        - maximum allowable value for Y_e
 *             : T_atm          - atmospheric value for temperature
 *             : T_min          - minimum allowable value for temperature
 *             : T_max          - maximum allowable value for temperature
 *
 * Outputs     : eos            - eos_parameters struct with the above inputs
 *                                initialized
 */
void initialize_tabulated_eos(
      const char *table_filepath,
      const double W_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      const double T_atm,
      const double T_min,
      const double T_max,
      eos_parameters *restrict eos ) {

  // Step 1: Set EOS type to Tabulated.
  eos->eos_type = grhayl_eos_hybrid;

  // Step 2: Read the EOS table
  eos->tabulated_read_table_set_EOS_params(table_filepath, eos);

  // Step 3: Initialize quantities which are common to all EOSs.
  init_common_eos_quantities;

  // Step 4: Set parameters specific to Tabulated EOS.
  eos->Ye_atm                 = Ye_atm;
  eos->Ye_min                 = Ye_min;
  eos->Ye_max                 = Ye_max;
  eos->T_atm                  = T_atm;
  eos->T_min                  = T_min;
  eos->T_max                  = T_max;
  eos->tabulated_compute_P_eps_S_from_T(eos,
                                        eos->rho_atm,
                                        Ye_atm, T_atm,
                                        &eos->press_atm,
                                        &eos->eps_atm,
                                        &eos->entropy_atm);

  // Step 5: These parameters are manually set here, but
  //         can be overwritten later.
  eos->root_finding_precision = 1e-15;
  eos->depsdT_threshold       = 1e-6;
}

/*
 * Function    : initialize_hybrid_eos_functions_and_params()
 * Description : Fully initializes EOS struct elements for a hybrid EOS
 *
 * Inputs      : W_max          - Maximum allowed Lorentz factor
 *             : rho_atm        - atmospheric value for rho
 *             : rho_min        - minimum allowable value for rho
 *             : rho_max        - maximum allowable value for rho
 *             : neos           - number of pieces in the piecewise
 *                                polytrope
 *             : rho_ppoly      - pointer to the array containing the
 *                                minimum rho_b for each polytropic piece
 *             : Gamma_ppoly    - pointer to the array containing the
 *                                minimum rho_b for each polytropic piece
 *             : K_ppoly0       - TODO: comment
 *             : Gamma_th       - TODO: comment
 *
 * Outputs     : eos            - eos_parameters struct with the above inputs
 *                                initialized
 */
void initialize_hybrid_eos_functions_and_params(
      const double W_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      eos_parameters *restrict eos ) {

  // Step 1: Initialize Hybrid EOS functions
  initialize_eos_functions(grhayl_eos_hybrid, eos);

  // Step 2: Initialize Hybrid EOS parameters
  initialize_hybrid_eos(W_max, rho_atm, rho_min, rho_max,
                        neos, rho_ppoly, Gamma_ppoly,
                        K_ppoly0, Gamma_th, eos);
}

/* Function    : initialize_tabulated_eos()
 * Description : Initializes EOS struct elements for tabulated EOS
 *
 * Inputs      : W_max          - maximum allowed Lorentz factor
 *             : rho_atm        - atmospheric value for rho
 *             : rho_min        - minimum allowable value for rho
 *             : rho_max        - maximum allowable value for rho
 *             : Y_e_atm        - atmospheric value for Y_e
 *             : Y_e_min        - minimum allowable value for Y_e
 *             : Y_e_max        - maximum allowable value for Y_e
 *             : T_atm          - atmospheric value for temperature
 *             : T_min          - minimum allowable value for temperature
 *             : T_max          - maximum allowable value for temperature
 *
 * Outputs     : eos            - eos_parameters struct with the above inputs
 *                                initialized
 */
void initialize_tabulated_eos_functions_and_params(
      const char *table_filepath,
      const double W_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      const double T_atm,
      const double T_min,
      const double T_max,
      eos_parameters *restrict eos ) {

  // Step 1: Initialize Tabulated EOS functions
  initialize_eos_functions(grhayl_eos_tabulated, eos);

  // Step 2: Initialize Tabulated EOS parameters
  initialize_tabulated_eos(table_filepath, W_max,
                           rho_atm, rho_min, rho_max,
                           Ye_atm, Ye_min, Ye_max,
                           T_atm, T_min, T_max,
                           eos);
}
