#include "cctk.h"
#include "con2prim_header.h"

/* The functions initialize_general_eos, initialize_hybrid_eos, and initialize_tabulated_eos
   fill the struct eos_parameters with data. Depending on which eos is being using, an eos-
   specific function must be called to fill in the eos parameters. The initialize_general_eos
   should be called regardless of eos, as it contains the common parameters that are generally
   needed. For more information on the arguments, see the definition of the struct in new_header.h. */
void initialize_general_eos(eos_parameters *restrict eos, const int type,
                const double tau_atm, const double W_max,
                const double eps_atm, const double eps_min, const double eps_max,
                const double press_atm, const double press_min, const double press_max,
                const double entropy_atm, const double entropy_min, const double entropy_max,
                const double rho_atm, const double rho_min, const double rho_max) {
  eos->eos_type = type;
  eos->tau_atm = tau_atm;
  eos->W_max = W_max;
  eos->inv_W_max_squared = 1.0/W_max;
  eos->eps_atm = eps_atm;
  eos->eps_min = eps_min;
  eos->eps_max = eps_max;
  eos->press_atm = press_atm;
  eos->press_min = press_min;
  eos->press_max = press_max;
  eos->entropy_atm = entropy_atm;
  eos->entropy_min = entropy_min;
  eos->entropy_max = entropy_max;
  eos->rho_atm = rho_atm;
  eos->rho_min = rho_min;
  eos->rho_max = rho_max;
}

void initialize_hybrid_eos(eos_parameters *restrict eos, const int neos,
                const double rho_ppoly_tab[], const double Gamma_ppoly_tab[],
                const double K_ppoly_tab[], const double eps_integ_const[],
                const double Gamma_th) {
  eos->neos = neos;
  for (int i=0; i<neos-1; i++) {
    eos->rho_ppoly_tab[i] = rho_ppoly_tab[i];
    eos->Gamma_ppoly_tab[i] = Gamma_ppoly_tab[i];
    eos->K_ppoly_tab[i] = K_ppoly_tab[i];
    eos->eps_integ_const[i] = eps_integ_const[i];
  }
  eos->Gamma_ppoly_tab[neos-1] = Gamma_ppoly_tab[neos-1];
  eos->K_ppoly_tab[neos-1] = K_ppoly_tab[neos-1];
  eos->eps_integ_const[neos-1] = eps_integ_const[neos-1];
  if(neos==1)
    eos->eps_integ_const[0] = 0.0;
  eos->Gamma_th = Gamma_th;
}

void initialize_tabulated_eos(eos_parameters *restrict eos, const double precision, const double threshold,
                const double temp_atm, const double temp_min, const double temp_max,
                const double Ye_atm, const double Ye_min, const double Ye_max) {
  eos->temp_atm = temp_atm;
  eos->temp_min = temp_min;
  eos->temp_max = temp_max;
  eos->Ye_atm = Ye_atm;
  eos->Ye_min = Ye_min;
  eos->Ye_max = Ye_max;
  eos->root_finding_precision = precision;
  eos->depsdT_threshold = threshold;
}
