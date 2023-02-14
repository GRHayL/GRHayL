#include "NRPyEOS_Hybrid.h"
#include "NRPyEOS_Tabulated.h"

// The initialize_general_eos() function sets the parameters in the eos_parameters struct
// which are independent of EOS. The functions initialize_hybrid_functions() and
// initialize_tabulated_functions() set the function pointers for eos_parameters. These
// must be run before initialize_hybrid_eos() or initialize_tabulated_eos(), respectively.
// To replace these functions with in-house variants, simply set the pointers to the new
// functions. The initialize_hybrid_eos() or initialize_tabulated_eos() funcitons set
// the parameters for hybrid and tabulated EOS, respectively. For more information on
// the arguments and other properties of eos_parameters, see the definition of the struct
// in GRHayl.h.

/* Function    : initialize_general_eos()
 * Description : Initializes EOS struct elements which are independent of the type of EOS
 *
 * Inputs      : type           - type of EOS (0=Hybrid, 1=Tabulated)
 *             : tau_atm        - atmospheric value for \tilde{tau} TODO: give definition of tau
 *             : W_max          - maximum allowable Lorenz factor W
 *             : rho_atm        - atmospheric value for rho_*
 *             : rho_min        - minimum allowable value for rho_*
 *             : rho_max        - maximum allowable value for rho_*
 *
 * Outputs     : eos            - eos_parameters struct with the above inputs
 *                                initialized
 *
 */

void initialize_general_eos(
      const int type,
      const double W_max,
      const double rho_atm,
      const double rho_min,
      const double rho_max,
      eos_parameters *restrict eos){
  eos->eos_type          = type;
  eos->W_max             = W_max;
  eos->inv_W_max_squared = 1.0/SQR(W_max);
  eos->rho_atm           = rho_atm;
  eos->rho_min           = rho_min;
  eos->rho_max           = rho_max;
}

/* Function    : initialize_hybrid_eos()
 * Description : Initializes EOS struct elements for a hybrid EOS
 *
 * Inputs      : neos           - number of pieces in the piecewise
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
 *
 */

void initialize_hybrid_eos(
      const int neos,
      const double *restrict rho_ppoly,
      const double *restrict Gamma_ppoly,
      const double K_ppoly0,
      const double Gamma_th,
      eos_parameters *restrict eos) {

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

  // Initialize {K_{j}}, j>=1, and {eps_integ_const_{j}}
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

//TODO: Eventually, improve this using initialize_Tabulated_EOS_parameters_from_input()
/* Function    : initialize_tabulated_eos()
 * Description : Initializes EOS struct elements for tabulated EOS
 *
 * Inputs      : root_finding_precision - TODO:
 *             : depsdT_threshold       - TODO:
 *             : Y_e_atm                - atmospheric value for Y_e
 *             : Y_e_min                - minimum allowable value for Y_e
 *             : Y_e_max                - maximum allowable value for Y_e
 *             : T_atm                  - atmospheric value for temperature
 *             : T_min                  - minimum allowable value for temperature
 *             : T_max                  - maximum allowable value for temperature
 *
 * Outputs     : eos                    - eos_parameters struct with the above inputs
 *                                        initialized
 *
 */

void initialize_tabulated_eos(
      const double root_finding_precision,
      const double depsdT_threshold,
      const double Ye_atm,
      const double Ye_min,
      const double Ye_max,
      const double T_atm,
      const double T_min,
      const double T_max,
      eos_parameters *restrict eos) {

  eos->root_finding_precision = root_finding_precision;
  eos->depsdT_threshold       = depsdT_threshold;
  eos->Ye_atm                 = Ye_atm;
  eos->Ye_min                 = Ye_min;
  eos->Ye_max                 = Ye_max;
  eos->T_atm                  = T_atm;
  eos->T_min                  = T_min;
  eos->T_max                  = T_max;
  eos->tabulated_compute_P_eps_S_from_T(eos, eos->rho_atm, Ye_atm, T_atm, &eos->press_atm, &eos->eps_atm, &eos->entropy_atm);
}

//TODO: Eventually, improve this using initialize_Tabulated_EOS_parameters_from_input()
/* Function    : initialize_hybrid_eos()
 * Description : Sets EOS hybrid function pointers to point to NRPyEOS
 *
 * Outputs     : eos            - eos_parameters struct hybrid function pointers
 *                                pointing to NRPyEOS
 *
 */

void initialize_hybrid_functions(eos_parameters *restrict eos) {
  NRPyEOS_initialize_hybrid_functions(eos);
}

//TODO: Eventually, improve this using initialize_Tabulated_EOS_parameters_from_input()
/* Function    : initialize_tabulated_functions()
 * Description : Sets EOS tabulated function pointers to point to NRPyEOS
 *
 * Outputs     : eos            - eos_parameters struct tabulated function pointers
 *                                pointing to NRPyEOS
 *
 */

void initialize_tabulated_functions(eos_parameters *restrict eos) {
  NRPyEOS_initialize_tabulated_functions(eos);
}
