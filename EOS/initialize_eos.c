#include "con2prim_header.h"
#include "Hybrid/EOS_hybrid_header.h"

/* The functions initialize_general_eos, initialize_hybrid_eos, and initialize_tabulated_eos
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
  eos->inv_W_max_squared = 1.0/W_max;
  eos->entropy_atm = entropy_atm;
  eos->entropy_min = entropy_min;
  eos->entropy_max = entropy_max;
  eos->rho_atm = rho_atm;
//  eos->rho_min = rho_min;
  eos->rho_max = rho_max;
}

void initialize_hybrid_eos(const int neos,
                const double rho_ppoly_tab[], const double Gamma_ppoly_tab[],
                const double K_ppoly_tab0, const double Gamma_th,
                eos_parameters *restrict eos) {

  eos->neos = neos;
  eos->Gamma_th = Gamma_th;
  eos->K_ppoly_tab[0] = K_ppoly_tab0;
  if( neos==1 ) {
    eos->rho_ppoly_tab[0] = rho_ppoly_tab[0];
    eos->eps_integ_const[0] = 0.0;
  } else {
    for(int j=0; j<=neos-2; j++) eos->rho_ppoly_tab[j] = rho_ppoly_tab[j];
  }
  for(int j=0; j<=neos-1; j++) eos->Gamma_ppoly_tab[j] = Gamma_ppoly_tab[j];

  // Initialize {K_{j}}, j>=1, and {eps_integ_const_{j}}
  setup_K_ppoly_tab_and_eps_integ_consts(eos);

  // --------- Atmospheric values ---------
  // Compute atmospheric P and eps
  double press, eps;
  compute_P_cold_and_eps_cold(eos, eos->rho_atm, &press, &eps);
  // Set atmospheric values
  eos->press_atm = press;
  eos->eps_atm = eps;
//Leo computes tau_atm like this, but this is a parameter as well. We should choose one.
//  eos->tau_atm = eos->rho_atm * eos->eps_atm;
  // --------------------------------------

  // -------------- Ceilings --------------
  // Compute maximum P and eps
  compute_P_cold_and_eps_cold(eos, eos->rho_max, &press, &eps);
  // Set maximum values
  eos->press_max = press;
  eos->eps_max = eps;
  // --------------------------------------

  // --------------- Floors ---------------
  // We'll choose the minimal values to be atmosphere
  eos->rho_min = eos->rho_atm;
  eos->press_min   = eos->press_atm;
  eos->eps_min = eos->eps_atm;
  // --------------------------------------
}


//Eventually, improve this using initialize_Tabulated_EOS_parameters_from_input()
void initialize_tabulated_eos(const double precision, const double threshold,
                const double temp_atm, const double temp_min, const double temp_max,
                const double Ye_atm, const double Ye_min, const double Ye_max,
                eos_parameters *restrict eos) {
  eos->temp_atm = temp_atm;
  eos->temp_min = temp_min;
  eos->temp_max = temp_max;
  eos->Ye_atm = Ye_atm;
  eos->Ye_min = Ye_min;
  eos->Ye_max = Ye_max;
  eos->root_finding_precision = precision;
  eos->depsdT_threshold = threshold;
}
