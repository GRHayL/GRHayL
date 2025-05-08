#include "ghl.h"
#include "ghl_radiation.h"

void ghl_m1_set_equilibrium(
    const double rho,
    const double T,
    const double eta_nue,
    const double eta_nuae,
    const double eta_nux,
    const double w_lorentz,
    const double *u4D,
    const double *n4U,
    const ghl_radiation_metric_tensor *proj4,
    const ghl_ADM_aux_quantities *adm_aux,
    const ghl_metric_quantities *metric,
    double *rN,
    double *rE) {

  const double hc = 1.24e-4 //hc in units of MeVcm

  //Calculate neutrino density (0 for number, 1 for energy)
  //eta refers to the chemical potential, not the emissions.
  double nue_dens0 = 4*M_PI/pow(hc, 3)*pow(T, 3)*NRPyLeakage_Fermi_Dirac_integrals(2, eta_nue);
  double nuae_dens0 = 4*M_PI/pow(hc, 3)*pow(T, 3)*NRPyLeakage_Fermi_Dirac_integrals(2, eta_nuae);
  double nux_dens0 = 16*M_PI/pow(hc, 3)*pow(T, 3)*NRPyLeakage_Fermi_Dirac_integrals(2, eta_nux);
  double nue_dens1 = 4*M_PI/pow(hc, 3)*pow(T, 4)*NRPyLeakage_Fermi_Dirac_integrals(3, eta_nue);
  double nuae_dens1 = 4*M_PI/pow(hc, 3)*pow(T, 4)*NRPyLeakage_Fermi_Dirac_integrals(3, eta_nuae);
  double nux_dens1 = 16*M_PI/pow(hc, 3)*pow(T, 4)*NRPyLeakage_Fermi_Dirac_integrals(3, eta_nux);

  //Now set to equilibrium in the fluid frame
  const double J = (nue_dens1 + nuae_dens1 + nux_dens1)*metric->detgamma; //FIXME: Is this right? I need a single J for the TDD calc, but I have 3 different dens1's?
  const double chi = 1.0/3.0; //Assuming Eddington closure

  //Compute Tmunu in lab frame from equilibrium.
  ghl_stress_energy rT4DD = {0};
  for (int i = 0; i < 4; ++i){
  for(int j = 0; j < 4; ++j) {
    rT4DD.T4[i][j] = ((4.0 / 3.0) * J * u4D[i] * u4D[j])
                     + ((1.0 / 3.0) * J * adm_aux->g4DD[i][j]);
  }
  }

  *rN = max(w_lorentz*nu_dens*metric->detgamma, Nmin); //FIXME: Same question as with J's calc.
  *rE = calc_J_from_rT(&n4U, &proj4, &rT4DD); //This works as long as you pass the normal vectors, rather than the four-velocity. Compare eqs (1) and (2) in Radice 2021
  // TODO: why is this being done?
  // ghl_radiation_flux_vector F4 = {0};
  // ghl_radiation_pressure_tensor P4 = {0};
  // calc_H4D_from_rT(&n4U, &proj4, &rT4DD, &F4);
  // calc_K4DD_from_rT(&n4U, &proj4, &rT4DD, &P4);
}
