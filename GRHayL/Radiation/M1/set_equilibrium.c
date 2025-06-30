#include "ghl.h"
#include "ghl_radiation.h"

void ghl_calc_neutrino_densities(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Y_e,
      const double T,
      double *restrict nudens_0,
      double *restrict nudens_1) {

  // Compute chemical potentials, needed to compute eta's
  double press_ = 0.0, eps_ = 0.0, muhat_ = 0.0; // unused stuff
  double mu_e = 0.0, mu_p = 0.0, mu_n = 0.0;
  ghl_tabulated_compute_P_eps_muhat_mue_mup_mun_from_T(
        eos, rho, Y_e, T, &press_, &eps_, &muhat_, &mu_e, &mu_p, &mu_n);

  // Compute eta's
  const double eta_nue = (mu_p + mu_e - mu_n) / T;
  const double eta_nua = -eta_nue;
  const double eta_nux = 0.0;

  // This follows the logic from NeutrinoDens_cgs.F90 from thcextra/WeakRates
  const double hc3 = 1.905895198207215e-30; // (hc)^3 in units of (MeVcm)^3
  const double prefactor = 4 * M_PI / hc3 * pow(T, 3);

  // These are the same as n_nue, n_nua, n_nux
  nudens_0[0] = prefactor * NRPyLeakage_Fermi_Dirac_integrals(2, eta_nue);
  nudens_0[1] = prefactor * NRPyLeakage_Fermi_Dirac_integrals(2, eta_nua);
  nudens_0[2] = prefactor * NRPyLeakage_Fermi_Dirac_integrals(2, eta_nux) * 4;

  // There are the same as en_nue, en_nua, en_nux
  nudens_1[0] = prefactor * T * NRPyLeakage_Fermi_Dirac_integrals(3, eta_nue);
  nudens_1[1] = prefactor * T * NRPyLeakage_Fermi_Dirac_integrals(3, eta_nua);
  nudens_1[2] = prefactor * T * NRPyLeakage_Fermi_Dirac_integrals(3, eta_nux) * 4;

  // nudens_0 has units of 1/cm^3. nudens_1 has units of MeV/cm3. Here we change
  // them to geometrized units.
  for(int i = 0; i < 3; i++) {
    nudens_0[i] *= NRPyLeakage_units_cgs_to_geom_nudens_0;
    nudens_1[i] *= NRPyLeakage_units_cgs_to_geom_nudens_1;
  }
}

/*
 * rho       : baryonic density (HydroBase::rho)
 *  T        : temperature (HydroBase::temperature)
 * num_dens  : number density
 * e_dens    :
 * w_lorentz : Lorentz factor (HydroBase::w_lorentz)
 * u4D       : fluid 4-velocity (use GRHayL struct and u^i = u^0 v^i)
 * n4U       : normal vector
 * proj4     : P^mu_nu = u^mu u_nu
 */
void ghl_m1_set_equilibrium(
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims,
      const ghl_metric_quantities *restrict adm_metric,
      const ghl_ADM_aux_quantities *restrict aux_metric,
      const double nudens_0,
      const double nudens_1,
      double *restrict rN,
      double *restrict rnnu,
      double *restrict rE,
      double *restrict rJ,
      ghl_radiation_flux_vector *restrict rF4,
      ghl_radiation_flux_vector *restrict rH4,
      ghl_radiation_pressure_tensor *restrict rP4) {

  // Compute fluid 4-velocity u^a
  const double u4U[4] = {
    prims->u0,
    prims->u0 * prims->vU[0],
    prims->u0 * prims->vU[1],
    prims->u0 * prims->vU[2],
  };

  // Compute fluid 4-velocity u_a = g_ab u^b
  double u4D[4] = { 0.0, 0.0, 0.0, 0.0 };
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      u4D[a] += aux_metric->g4DD[a][b] * u4U[b];
    }
  }

  // Now compute "projection" operator P^a_b = u^a u_b
  ghl_radiation_metric_tensor proj4 = { 0 };
  for(int a = 0; a < 4; a++) {
    for(int b = 0; b < 4; b++) {
      proj4.UD[a][b] = (int)(a == b) + u4U[a] * u4D[b];
    }
  }

  // Set unit normal vector
  const double n4U[4] = {
    +adm_metric->lapseinv,
    -adm_metric->lapseinv * adm_metric->betaU[0],
    -adm_metric->lapseinv * adm_metric->betaU[1],
    -adm_metric->lapseinv * adm_metric->betaU[2],
  };

  // Now set to equilibrium in the fluid frame
  *rJ = nudens_1 * adm_metric->detgamma;
  rH4->D[0] = 0.0;
  rH4->D[1] = 0.0;
  rH4->D[2] = 0.0;
  rH4->D[3] = 0.0;

  // W = - n_a u^a = alpha u^0, since n_a = (-alpha, 0, 0, 0).
  const double W_lorentz = adm_metric->lapse * prims->u0;

  // Compute Tmunu in lab frame from equilibrium.
  ghl_stress_energy rT4DD = { 0 };
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < 4; ++j) {
      rT4DD.T4[i][j] = ((4.0 / 3.0) * (*rJ) * u4D[i] * u4D[j])
                       + ((1.0 / 3.0) * (*rJ) * aux_metric->g4DD[i][j]);
    }
  }

  // Update the GF's for the equilibrium data
  *rnnu = nudens_0 * adm_metric->detgamma;
  *rN = fmax(W_lorentz * (*rnnu), 0.0); // TODO(DRB): Add in Nfloor.
  *rE = calc_J_from_rT(
        n4U, &proj4,
        &rT4DD); // This works as long as you pass the normal vectors, rather than the
                 // four-velocity. Compare eqs (1) and (2) in Radice 2021
  calc_H4D_from_rT(n4U, &proj4, &rT4DD, rF4);
  calc_K4DD_from_rT(n4U, &proj4, &rT4DD, rP4);
}
