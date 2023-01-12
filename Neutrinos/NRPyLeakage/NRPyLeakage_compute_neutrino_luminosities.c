#include "Neutrinos.h"

static double EnsureFinite(const double x) {
  if(robust_isfinite(x))
    return x;
  else
    return 1e-15;
}

/*
 * (c) Leo Werneck
 * Compute neutrino luminosities following Siegel & Metzger (2018)
 * https://arxiv.org/pdf/1711.00868.pdf
 * Neutrino rates: https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
 */
void NRPyLeakage_compute_neutrino_luminosities(
      const eos_parameters *restrict eos,
      const double alpha,
      const double gammaxx,
      const double gammaxy,
      const double gammaxz,
      const double gammayy,
      const double gammayz,
      const double gammazz,
      const double rho,
      const double Y_e,
      const double T,
      const double W,
      const neutrino_optical_depths *restrict tau,
      neutrino_luminosities *restrict lum ) {

  // Step 1: Get chemical potentials and mass
  //         fractions using the EOS
  double muhat, mu_e, mu_p, mu_n, X_n, X_p;
  eos->tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(eos, rho, Y_e, T, &muhat, &mu_e, &mu_p, &mu_n, &X_n, &X_p);

  // Step 2: Compute rho in cgs units
  const double rho_cgs = rho * NRPyLeakage_units_geom_to_cgs_D;

  // Step 3: Compute Y_{pn} and Y_{np}
  const double Y_p = Y_e;
  const double Y_n = 1-Y_e;
  const double exp_metahat = exp(-muhat/T);
  // Step 3.a: Compute Y_{np}
  double Y_np = (Y_e < 0.5) ? (2.0*Y_e-1.0)/(exp_metahat-1.0) : Y_n;

  // Step 3.b: Compute Y_{pn}
  double Y_pn = (Y_e > 0.5) ? exp_metahat*(2.0*Y_e-1.0)/(exp_metahat-1.0) : Y_p;

  // Step 3.c: Make sure both Y_np and Y_pn are non-negative
  if( Y_np < 0.0 ) Y_np = Y_n;
  if( Y_pn < 0.0 ) Y_pn = Y_p;

  // Step 4: Compute the source terms
  //         Note: The code below is generated by NRPy+
  const double tmp_0 = (1.0/(T));
  const double tmp_1 = mu_e*tmp_0;
  const double tmp_2 = NRPyLeakage_Fermi_Dirac_integrals(4, tmp_1);
  const double tmp_3 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_1)/tmp_2;
  const double tmp_4 = exp(-tau->nue[0]);
  const double tmp_6 = NRPyLeakage_eta_nue_0*tmp_4 + (1 - tmp_4)*(-muhat*tmp_0 + tmp_1);
  const double tmp_7 = ((NRPyLeakage_alpha)*(NRPyLeakage_alpha));
  const double tmp_8 = M_PI/NRPyLeakage_hc3;
  const double tmp_10 = 8*NRPyLeakage_N_A*NRPyLeakage_beta*((T)*(T)*(T)*(T)*(T))*rho_cgs*tmp_8*((3.0/8.0)*tmp_7 + 1.0/8.0);
  const double tmp_11 = EnsureFinite(NRPyLeakage_Brems_C2*NRPyLeakage_enable_brems_nui_anui*T*EnsureFinite(NRPyLeakage_Brems_C1*NRPyLeakage_Brems_zeta*pow(T, 4.5)*rho_cgs*(((X_n)*(X_n)) + (28.0/3.0)*X_n*X_p + ((X_p)*(X_p))))/NRPyLeakage_Brems_C1);
  const double tmp_12 = exp(-tau->anue[0]);
  const double tmp_14 = NRPyLeakage_eta_anue_0*tmp_12 + (1 - tmp_12)*(muhat*tmp_0 - tmp_1);
  const double tmp_15 = ((M_PI)*(M_PI));
  const double tmp_17 = (1.0/3.0)*tmp_15 + ((mu_e)*(mu_e))/((T)*(T));
  const double tmp_18 = NRPyLeakage_gamma_0*sqrt(tmp_17);
  const double tmp_20 = ((NRPyLeakage_gamma_0)*(NRPyLeakage_gamma_0))*tmp_17/(tmp_18 + 1);
  const double tmp_21 = -1.0/2.0*tmp_20 - 1;
  const double tmp_23 = (1.0/((NRPyLeakage_hc3)*(NRPyLeakage_hc3)));
  const double tmp_24 = pow(T, 8);
  const double tmp_26 = (1.0/3.0)*((M_PI)*(M_PI)*(M_PI))*NRPyLeakage_beta*pow(NRPyLeakage_gamma_0, 6)*((tmp_17)*(tmp_17)*(tmp_17))*tmp_23*tmp_24*(tmp_18 + 1)*exp(-tmp_18)/NRPyLeakage_alpha_fs;
  const double tmp_27 = (1.0/2.0)*T*(tmp_20 + 2);
  const double tmp_28 = NRPyLeakage_Fermi_Dirac_integrals(3, tmp_1);
  const double tmp_29 = (1.0/(tmp_28));
  const double tmp_30 = NRPyLeakage_Fermi_Dirac_integrals(4, -tmp_1);
  const double tmp_31 = NRPyLeakage_Fermi_Dirac_integrals(3, -tmp_1);
  const double tmp_32 = (1.0/(tmp_31));
  const double tmp_33 = -1.0/2.0*tmp_2*tmp_29 - 1.0/2.0*tmp_30*tmp_32;
  const double tmp_34 = tmp_15*tmp_23*tmp_28;
  const double tmp_35 = (16.0/9.0)*NRPyLeakage_beta*tmp_24*tmp_31*tmp_34;
  const double tmp_36 = 32*pow(T, 9);
  const double tmp_37 = (1.0/64.0)*((NRPyLeakage_hc3)*(NRPyLeakage_hc3))*tmp_29*tmp_32*(tmp_15*tmp_2*tmp_23*tmp_31*tmp_36 + tmp_30*tmp_34*tmp_36)/(tmp_15*tmp_24);
  const double tmp_38 = tmp_11 + EnsureFinite(NRPyLeakage_enable_pair_nue_anue*tmp_37*EnsureFinite(NRPyLeakage_C1pC2_nue_anue*tmp_35/((exp(tmp_14 + tmp_33) + 1)*(exp(tmp_33 + tmp_6) + 1)))) + EnsureFinite(NRPyLeakage_enable_plasmon_nue_anue*tmp_27*EnsureFinite(((NRPyLeakage_C_V)*(NRPyLeakage_C_V))*tmp_26/((exp(tmp_14 + tmp_21) + 1)*(exp(tmp_21 + tmp_6) + 1))));
  const double tmp_39 = tmp_38 + EnsureFinite(NRPyLeakage_enable_beta_nue*T*tmp_3*EnsureFinite(Y_pn*tmp_10*tmp_2/(exp(-tmp_3 + tmp_6) + 1)));
  const double tmp_40 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_6);
  const double tmp_41 = NRPyLeakage_Fermi_Dirac_integrals(3, tmp_6);
  const double tmp_42 = NRPyLeakage_N_A*NRPyLeakage_sigma_0*((T)*(T))*rho_cgs/((NRPyLeakage_m_e_c2)*(NRPyLeakage_m_e_c2));
  const double tmp_43 = tmp_40*tmp_42/tmp_41;
  const double tmp_45 = (1 - Y_e)*((5.0/24.0)*tmp_7 + 1.0/24.0)/((2.0/3.0)*MAX(mu_n*tmp_0, 0) + 1);
  const double tmp_46 = ((NRPyLeakage_C_V - 1)*(NRPyLeakage_C_V - 1));
  const double tmp_47 = Y_e*((1.0/6.0)*tmp_46 + (5.0/24.0)*tmp_7)/((2.0/3.0)*MAX(mu_p*tmp_0, 0) + 1);
  const double tmp_48 = (3.0/4.0)*tmp_7 + 1.0/4.0;
  const double tmp_49 = 4*((T)*(T)*(T)*(T))*tmp_8;
  const double tmp_50 = 6/NRPyLeakage_units_geom_to_cgs_L;
  const double tmp_52 = NRPyLeakage_units_cgs_to_geom_Q*W*((alpha)*(alpha))*sqrt(gammaxx*gammayy*gammazz - gammaxx*((gammayz)*(gammayz)) - ((gammaxy)*(gammaxy))*gammazz + 2*gammaxy*gammaxz*gammayz - ((gammaxz)*(gammaxz))*gammayy);
  const double tmp_53 = NRPyLeakage_Fermi_Dirac_integrals(5, -tmp_1)/tmp_30;
  const double tmp_54 = tmp_38 + EnsureFinite(NRPyLeakage_enable_beta_anue*T*tmp_53*EnsureFinite(Y_np*tmp_10*tmp_30/(exp(tmp_14 - tmp_53) + 1)));
  const double tmp_55 = NRPyLeakage_Fermi_Dirac_integrals(5, tmp_14);
  const double tmp_56 = NRPyLeakage_Fermi_Dirac_integrals(3, tmp_14);
  const double tmp_57 = tmp_55/tmp_56;
  const double tmp_60 = NRPyLeakage_Fermi_Dirac_integrals(3, 0);
  const double tmp_61 = NRPyLeakage_Fermi_Dirac_integrals(5, 0)/tmp_60;
  lum->nue = tmp_39*tmp_52/(((tau->nue[1])*(tau->nue[1]))*tmp_39*tmp_50/((EnsureFinite(tmp_43*tmp_45) + EnsureFinite(tmp_43*tmp_47) + EnsureFinite(Y_np*tmp_43*tmp_48/(exp(tmp_1 - tmp_40/NRPyLeakage_Fermi_Dirac_integrals(4, tmp_6)) + 1)))*MAX(tmp_41*tmp_49, 1.0000000000000001e-15)) + 1);
  lum->anue = tmp_52*tmp_54/(((tau->anue[1])*(tau->anue[1]))*tmp_50*tmp_54/((EnsureFinite(tmp_42*tmp_45*tmp_57) + EnsureFinite(tmp_42*tmp_47*tmp_57) + EnsureFinite(Y_pn*tmp_42*tmp_48*tmp_57/(exp(-tmp_1 - tmp_55/NRPyLeakage_Fermi_Dirac_integrals(4, tmp_14)) + 1)))*MAX(tmp_49*tmp_56, 1.0000000000000001e-15)) + 1);
  lum->nux = tmp_52*(tmp_11 + EnsureFinite(NRPyLeakage_enable_pair_nux_anux*tmp_37*EnsureFinite(NRPyLeakage_C1pC2_nux_anux*tmp_35/((exp(tmp_33) + 1)*(exp(tmp_33) + 1)))) + EnsureFinite(NRPyLeakage_enable_plasmon_nux_anux*tmp_27*EnsureFinite(tmp_26*tmp_46/((exp(tmp_21) + 1)*(exp(tmp_21) + 1)))))/(((tau->nux[1])*(tau->nux[1]))*tmp_39*tmp_50/((EnsureFinite(tmp_42*tmp_45*tmp_61) + EnsureFinite(tmp_42*tmp_47*tmp_61))*MAX(tmp_49*tmp_60, 1.0000000000000001e-15)) + 1);
}
