#include "ghl.h"
#include "ghl_radiation.h"

void calc_neutrino_densities(const bool heavy,
                             const double eta,
                             const double T,
                             double *num_dens,
                             double *e_dens){
  
  const double hc = 1.24e-4; //hc in units of MeVcm

  if(heavy == false){
    *num_dens = 4*M_PI/pow(hc, 3)*pow(T, 3)*NRPyLeakage_Fermi_Dirac_integrals(2, eta);
    *e_dens = 4*M_PI/pow(hc, 3)*pow(T, 4)*NRPyLeakage_Fermi_Dirac_integrals(3, eta);
  }
  else{
    *num_dens = 16*M_PI/pow(hc, 3)*pow(T, 3)*NRPyLeakage_Fermi_Dirac_integrals(2, eta);
    *e_dens = 16*M_PI/pow(hc, 3)*pow(T, 4)*NRPyLeakage_Fermi_Dirac_integrals(3, eta);
  }

}

void ghl_m1_set_equilibrium(
    const double rho,
    const double T,
    const double num_dens,
    const double e_dens,
    const double w_lorentz,
    const double *u4D,
    const double *n4U,
    const ghl_radiation_metric_tensor *proj4,
    const ghl_ADM_aux_quantities *adm_aux,
    const ghl_metric_quantities *metric,
    double *rN,
    double *rnnu,
    double *rE,
    double *rJ,
    ghl_radiation_flux_vector *rF4,
    ghl_radiation_flux_vector *rH4,
    ghl_radiation_pressure_tensor *rP4) {

  //Now set to equilibrium in the fluid frame
  *rJ = e_dens*metric->detgamma;
  rH4->D[0] = 0.0;
  rH4->D[1] = 0.0;
  rH4->D[2] = 0.0;
  rH4->D[3] = 0.0;

  //Compute Tmunu in lab frame from equilibrium.
  ghl_stress_energy rT4DD = {0};
  for (int i = 0; i < 4; ++i){
    for(int j = 0; j < 4; ++j) {
      rT4DD.T4[i][j] = ((4.0 / 3.0) * (*rJ) * u4D[i] * u4D[j])
                       + ((1.0 / 3.0) * (*rJ) * adm_aux->g4DD[i][j]);
    }
  }

  //Update the GF's for the equilibrium data
  *rnnu = num_dens*metric->detgamma;
  *rN = fmax(w_lorentz*(*rnnu), 0.0); //TODO(DRB): Add in Nfloor.
  *rE = calc_J_from_rT(n4U, proj4, &rT4DD); //This works as long as you pass the normal vectors, rather than the four-velocity. Compare eqs (1) and (2) in Radice 2021
  calc_H4D_from_rT(n4U, proj4, &rT4DD, rF4);
  calc_K4DD_from_rT(n4U, proj4, &rT4DD, rP4);
}
