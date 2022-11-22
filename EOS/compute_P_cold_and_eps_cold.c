// These were taken from Leo Werneck's IllinoisGRMHD code for con2prim to work,
// but they should actually be provided by the EOS code.

#include "EOS_hybrid_header.h"

void compute_P_cold_and_eps_cold(const eos_parameters *restrict eos, const double rho_in,
                              double *restrict P_cold_ptr, double *restrict eps_cold_ptr) {

  if(rho_in==0) {
    *P_cold_ptr   = 0.0;
    *eps_cold_ptr = 0.0;
    return;
  }

  int polytropic_index      = find_polytropic_K_and_Gamma_index(eos,rho_in);
  double K_ppoly_tab     = eos->K_ppoly_tab[polytropic_index];
  double Gamma_ppoly_tab = eos->Gamma_ppoly_tab[polytropic_index];
  double eps_integ_const = eos->eps_integ_const[polytropic_index];

  double P_cold = K_ppoly_tab*pow(rho_in,Gamma_ppoly_tab);

  *eps_cold_ptr = P_cold/(rho_in*(Gamma_ppoly_tab-1.0)) + eps_integ_const;

  *P_cold_ptr = P_cold;
}
