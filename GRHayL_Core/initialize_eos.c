#include "GRHayL.h"

/* Function    : initialize_general_eos()
 * Authors     : Leo Werneck & Samuel Cupp
 * Description : This function initializes the quantities in the
 *               EOS struct which are independent of the type of EOS
 *
 * Inputs      : type           - the type of EOS (0=Hybrid, 1=Tabulated)
 *             : tau_atm        - the atmospheric value for the conservative
 *                                tau tilde TODO: give definition of tau
 *             : W_max          - the maximum allowable value for the Lorenz
 *                                factor W
 *             : entropy_atm    - the atmospheric value for the entropy TODO: entropy or densitized entropy?
 *             : entropy_min    - the minimum allowable value for the entropy
 *             : entropy_max    - the maximum allowable value for the entropy
 *             : rho_atm        - the atmospheric value for rho_*
 *             : rho_min        - the minimum allowable value for rho_*
 *             : rho_max        - the maximum allowable value for rho_*
 *
 * Outputs     : eos            - an eos_parameters struct with the above inputs
 *                                initialized
 The functions initialize_general_eos, initialize_hybrid_eos, and initialize_tabulated_eos
   fill the struct eos_parameters with data. Depending on which eos is being using, an eos-
   specific function must be called to fill in the eos parameters. The initialize_general_eos
   should be called regardless of eos, as it contains the common parameters that are generally
   needed. For more information on the arguments, see the definition of the struct in new_header.h. */
void initialize_general_eos(const int type,
                const double tau_atm, const double W_max,
                const double entropy_atm, const double entropy_min, const double entropy_max,
                const double rho_atm, const double rho_min, const double rho_max, //let explicit setting of min, or assume atm?
                eos_parameters *restrict eos) {
  eos->eos_type = type;
  eos->tau_atm = tau_atm;
  eos->W_max = W_max;
  eos->inv_W_max_squared = 1.0/SQR(W_max);
  eos->entropy_atm = entropy_atm;
  eos->entropy_min = entropy_min;
  eos->entropy_max = entropy_max;
  eos->rho_atm = rho_atm;
  eos->rho_min = eos->rho_atm;
//  eos->rho_min = rho_min;
  eos->rho_max = rho_max;
}

void initialize_hybrid_eos(const int neos,
                const double rho_ppoly[], const double Gamma_ppoly[],
                const double K_ppoly0, const double Gamma_th,
                eos_parameters *restrict eos) {

  // --------- Function Prototypes --------
  eos->hybrid_find_polytropic_index            = &GRHayL_find_polytropic_index;
  eos->hybrid_get_K_and_Gamma                  = &GRHayL_get_K_and_Gamma;
  eos->hybrid_set_K_ppoly_and_eps_integ_consts = &GRHayL_set_K_ppoly_and_eps_integ_consts;
  eos->hybrid_compute_P_cold                   = &GRHayL_compute_P_cold;
  eos->hybrid_compute_P_cold_and_eps_cold      = &GRHayL_compute_P_cold_and_eps_cold;
  eos->hybrid_compute_entropy_function         = &GRHayL_compute_entropy_function;
  // --------------------------------------

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
  (*eos->hybrid_set_K_ppoly_and_eps_integ_consts)(eos);

  // --------- Atmospheric values ---------
  // Compute atmospheric P and eps
  double press, eps;
  (*eos->hybrid_compute_P_cold_and_eps_cold)(eos, eos->rho_atm, &press, &eps);
  // Set atmospheric values
  eos->press_atm = press;
  eos->eps_atm = eps;
//Leo computes tau_atm like this, but this is a parameter as well. We should choose one.
//  eos->tau_atm = eos->rho_atm * eos->eps_atm;
  // --------------------------------------

  // -------------- Ceilings --------------
  // Compute maximum P and eps
  (*eos->hybrid_compute_P_cold_and_eps_cold)(eos, eos->rho_max, &press, &eps);
  // Set maximum values
  eos->press_max = press;
  eos->eps_max = eps;
  // --------------------------------------

  // --------------- Floors ---------------
  // We'll choose the minimal values to be atmosphere
  eos->press_min   = eos->press_atm;
  eos->eps_min = eos->eps_atm;
  // --------------------------------------
}


//Eventually, improve this using initialize_Tabulated_EOS_parameters_from_input()
void initialize_tabulated_eos(const double precision, const double threshold,
                const double T_atm, const double T_min, const double T_max,
                const double Ye_atm, const double Ye_min, const double Ye_max,
                eos_parameters *restrict eos) {
  eos->T_atm = T_atm;
  eos->T_min = T_min;
  eos->T_max = T_max;
  eos->Ye_atm = Ye_atm;
  eos->Ye_min = Ye_min;
  eos->Ye_max = Ye_max;
  eos->root_finding_precision = precision;
  eos->depsdT_threshold = threshold;
}
