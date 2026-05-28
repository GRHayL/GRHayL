#include "ghl_nrpyeos_tabulated.h"
#include "NRPyEOS_stellarcollapse.h"

void NRPyEOS_stellarcollapse_to_ghl(NRPyEOS_stellarcollapse_t *restrict sc, ghl_eos_parameters *restrict eos) {

  eos->N_rho = sc->n_rho;
  eos->N_T = sc->n_temperature;
  eos->N_Ye = sc->n_ye;

  const int npoints = eos->N_rho * eos->N_T * eos->N_Ye;

  eos->table_logrho = (double *)malloc(sizeof(double) * eos->N_rho);
  eos->table_logT = (double *)malloc(sizeof(double) * eos->N_T);
  eos->table_Y_e = (double *)malloc(sizeof(double) * eos->N_Ye);
  eos->table_all = (double *)malloc(sizeof(double) * npoints * NRPyEOS_ntablekeys);
  eos->table_eps = (double *)malloc(sizeof(double) * npoints);
  eos->table_logh = (double *)malloc(sizeof(double) * npoints);

  if(eos->table_logrho == NULL || eos->table_logT == NULL || eos->table_Y_e == NULL
     || eos->table_all == NULL || eos->table_eps == NULL || eos->table_logh == NULL) {
    ghl_error("Cannot allocate memory for EOS table\n");
  }

  for(int ir = 0; ir < eos->N_rho; ir++) {
    eos->table_logrho[ir] = sc->log10_rho[ir] * log(10.0) + log(CGS_TO_CODE_DENSITY);
  }
  for(int it = 0; it < eos->N_T; it++) {
    eos->table_logT[it] = sc->log10_temperature[it] * log(10.0);
  }
  for(int iy = 0; iy < eos->N_Ye; iy++) {
    eos->table_Y_e[iy] = sc->ye[iy];
  }

  eos->energy_shift = sc->energy_shift * CGS_TO_CODE_ENERGY;

  const int map_sc_to_ghl[NRPyEOS_sc_n_quantities] = {
    [NRPyEOS_sc_Abar] = NRPyEOS_Abar_key,
    [NRPyEOS_sc_Xa] = NRPyEOS_X_a_key,
    [NRPyEOS_sc_Xh] = NRPyEOS_X_h_key,
    [NRPyEOS_sc_Xn] = NRPyEOS_X_n_key,
    [NRPyEOS_sc_Xp] = NRPyEOS_X_p_key,
    [NRPyEOS_sc_Zbar] = NRPyEOS_Zbar_key,
    [NRPyEOS_sc_cs2] = NRPyEOS_cs2_key,
    [NRPyEOS_sc_dedt] = NRPyEOS_depsdT_key,
    [NRPyEOS_sc_dpderho] = NRPyEOS_dPdeps_key,
    [NRPyEOS_sc_dpdrhoe] = NRPyEOS_dPdrho_key,
    [NRPyEOS_sc_entropy] = NRPyEOS_entropy_key,
    [NRPyEOS_sc_gamma] = NRPyEOS_Gamma_key,
    [NRPyEOS_sc_logenergy] = NRPyEOS_eps_key,
    [NRPyEOS_sc_logpress] = NRPyEOS_press_key,
    [NRPyEOS_sc_mu_e] = NRPyEOS_mu_e_key,
    [NRPyEOS_sc_mu_n] = NRPyEOS_mu_n_key,
    [NRPyEOS_sc_mu_p] = NRPyEOS_mu_p_key,
    [NRPyEOS_sc_muhat] = NRPyEOS_muhat_key,
    [NRPyEOS_sc_munu] = NRPyEOS_munu_key,
  };

  for(int i = 0; i < npoints * NRPyEOS_ntablekeys; i++) {
    eos->table_all[i] = 0.0;
  }

  for(int n = 0; n < NRPyEOS_sc_n_quantities; n++) {
    const int key = map_sc_to_ghl[n];
    for(int i = 0; i < npoints; i++) {
      eos->table_all[key + NRPyEOS_ntablekeys * i] = sc->data[n][i];
    }
  }
}
