#include "ghl_m1_nrpyleakage_kernel.h"
#include "ghl_nrpyeos_tabulated.h"

static double ghl_m1_nrpyleakage_finite_rate(const double x,
                                      const bool strict_rates,
                                      bool *const rate_failure) {
  if(isfinite(x))
    return x;
  if(strict_rates)
    *rate_failure = true;
  return 1.0e-15;
}

static bool raw_value_is_valid(const double x) {
  return isfinite(x) && x >= 0.0;
}

/* The generated approximation is positive for every finite argument.  Do not
 * replace a caller's finite degeneracy with a tail value: a finite positive
 * polynomial can still be evaluated accurately at large positive arguments.
 * Conversely, an overflow or an underflow to zero means that the required
 * positive intermediate is not representable, so report the rate failure. */
static ghl_error_codes_t evaluate_fermi_dirac_integral(
      const int k, const double z, double *const integral) {
  if(integral == NULL || !isfinite(z))
    return ghl_error_m1_microphysics_failure;
  const ghl_error_codes_t err =
      NRPyLeakage_Fermi_Dirac_integrals(k, z, integral);
  if(err != ghl_success)
    return err;
  if(!isfinite(*integral) || *integral <= 0.0)
    return ghl_error_m1_microphysics_failure;
  return ghl_success;
}

/* Evaluate 1/(exp(x)+1) without forming an overflowing exp(x).  A positive
 * finite x whose reciprocal tail underflows cannot be represented as a
 * nonzero binary64 result; return an explicit microphysics failure instead of
 * silently turning that tail into a different rate. */
static ghl_error_codes_t evaluate_fermi_factor(
      const double x, double *const factor) {
  if(factor == NULL || !isfinite(x))
    return ghl_error_m1_microphysics_failure;
  if(x > 0.0) {
    const double exp_neg_x = exp(-x);
    if(!isfinite(exp_neg_x) || exp_neg_x == 0.0)
      return ghl_error_m1_microphysics_failure;
    *factor = exp_neg_x/(1.0 + exp_neg_x);
  } else {
    const double exp_x = exp(x);
    if(!isfinite(exp_x))
      return ghl_error_m1_microphysics_failure;
    *factor = 1.0/(1.0 + exp_x);
  }
  if(!isfinite(*factor))
    return ghl_error_m1_microphysics_failure;
  return ghl_success;
}

static ghl_error_codes_t validate_raw_species(
      const ghl_m1_nrpyleakage_species_raw_rates *restrict r) {
  if(!isfinite(r->neutrino_degeneracy) ||
     !raw_value_is_valid(r->F2) || !raw_value_is_valid(r->F3) ||
     !raw_value_is_valid(r->F4) || !raw_value_is_valid(r->F5) ||
     !raw_value_is_valid(r->n_eq_cgs) ||
     !raw_value_is_valid(r->J_eq_mev_cgs) ||
     !raw_value_is_valid(r->mean_energy_mev) ||
     !raw_value_is_valid(r->eta_N_beta_cgs) ||
     !raw_value_is_valid(r->eta_N_pair_cgs) ||
     !raw_value_is_valid(r->eta_N_plasmon_cgs) ||
     !raw_value_is_valid(r->eta_N_brems_cgs) ||
     !raw_value_is_valid(r->eta_E_beta_mev_cgs) ||
     !raw_value_is_valid(r->eta_E_pair_mev_cgs) ||
     !raw_value_is_valid(r->eta_E_plasmon_mev_cgs) ||
     !raw_value_is_valid(r->eta_E_brems_mev_cgs) ||
     !raw_value_is_valid(r->kappa_a_N_cc_cgs) ||
     !raw_value_is_valid(r->kappa_a_E_cc_cgs) ||
     !raw_value_is_valid(r->kappa_s_N_neutron_cgs) ||
     !raw_value_is_valid(r->kappa_s_N_proton_cgs) ||
     !raw_value_is_valid(r->kappa_s_E_neutron_cgs) ||
     !raw_value_is_valid(r->kappa_s_E_proton_cgs) ||
     r->F2 == 0.0 || r->mean_energy_mev <= 0.0)
    return ghl_error_m1_microphysics_failure;
  return ghl_success;
}

static ghl_error_codes_t validate_kernel_result(
      const ghl_m1_nrpyleakage_legacy_kernel_result *restrict result) {
  for(int s = 0; s < ghl_m1_nrpyleakage_species_count; ++s) {
    const ghl_error_codes_t err = validate_raw_species(&result->raw.species[s]);
    if(err != ghl_success)
      return err;
  }
  const double *const arrays[] = {
    result->eta_N, result->eta_E, result->kappa_N, result->kappa_E,
    result->n_eq, result->J_eq
  };
  for(size_t a = 0; a < sizeof(arrays)/sizeof(arrays[0]); ++a)
    for(int s = 0; s < ghl_m1_nrpyleakage_species_count; ++s)
      if(!isfinite(arrays[a][s]) || arrays[a][s] < 0.0)
        return ghl_error_m1_microphysics_failure;
  return ghl_success;
}

ghl_error_codes_t ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
      const double rho,
      const double Ye,
      const double T,
      const double muhat,
      const double mu_e,
      const double mu_p,
      const double mu_n,
      const double X_n,
      const double X_p,
      const bool strict_validation,
      ghl_m1_nrpyleakage_thermo_state *restrict thermo) {
  if(thermo == NULL)
    return ghl_error_m1_null_pointer;
  if(!isfinite(rho) || rho <= 0.0 || !isfinite(T) || T <= 0.0 ||
     !isfinite(Ye) || Ye < 0.0 || Ye > 1.0)
    return ghl_error_m1_microphysics_failure;

  ghl_m1_nrpyleakage_thermo_state candidate = {
    .rho = rho, .rho_cgs = rho*NRPyLeakage_units_geom_to_cgs_D,
    .T = T, .Ye = Ye, .muhat = muhat, .mu_e = mu_e,
    .mu_p = mu_p, .mu_n = mu_n, .X_n = X_n, .X_p = X_p
  };
  const double Y_p = Ye;
  const double Y_n = 1.0 - Ye;
  const double metahat_over_T = -candidate.muhat/T;
  /* Evaluate the two composition ratios without forming exp(metahat/T).
   * The direct generated expression becomes inf/inf for a proton-rich,
   * neutron-degenerate tail even though its finite limiting ratio is well
   * defined. */
  if(Ye < 0.5) {
    candidate.Y_np = (2.0*Ye - 1.0)/expm1(metahat_over_T);
    candidate.Y_pn = Y_p;
  } else if(Ye > 0.5) {
    candidate.Y_np = Y_n;
    candidate.Y_pn = metahat_over_T > 0.0
        ? (2.0*Ye - 1.0)/(-expm1(-metahat_over_T))
        : (2.0*Ye - 1.0)*exp(metahat_over_T)/expm1(metahat_over_T);
  } else {
    candidate.Y_np = Y_n;
    candidate.Y_pn = Y_p;
  }
  if(candidate.Y_np <= 0.0) candidate.Y_np = Y_n;
  if(candidate.Y_pn <= 0.0) candidate.Y_pn = Y_p;

  if(strict_validation) {
    const double values[] = {
      candidate.rho, candidate.rho_cgs, candidate.T, candidate.Ye,
      candidate.muhat, candidate.mu_e, candidate.mu_p, candidate.mu_n,
      candidate.X_n, candidate.X_p, candidate.Y_np, candidate.Y_pn
    };
    for(size_t i = 0; i < sizeof(values)/sizeof(values[0]); ++i)
      if(!isfinite(values[i]))
        return ghl_error_m1_microphysics_failure;
    if(candidate.rho_cgs <= 0.0 || candidate.X_n < 0.0 ||
       candidate.X_p < 0.0 || candidate.Y_np < 0.0 || candidate.Y_pn < 0.0)
      return ghl_error_m1_microphysics_failure;
  }
  *thermo = candidate;
  return ghl_success;
}

static ghl_error_codes_t compute_thermo_state_impl(
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double Ye,
      const double T,
      ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const bool strict_validation,
      const bool use_legacy_eos_callback) {
  if(eos == NULL || thermo == NULL)
    return ghl_error_m1_null_pointer;
#ifdef GHL_DISABLE_HDF5
  (void)rho; (void)Ye; (void)T;
  (void)strict_validation; (void)use_legacy_eos_callback;
  return ghl_error_used_disabled_hdf5;
#else
  double muhat, mu_e, mu_p, mu_n, X_n, X_p;
  ghl_error_codes_t err;
  if(use_legacy_eos_callback) {
    if(ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T == NULL)
      return ghl_error_m1_microphysics_failure;
    err = ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T(
        eos, rho, Ye, T, &muhat, &mu_e, &mu_p, &mu_n, &X_n, &X_p);
  }
  else {
    if(eos->eos_type != ghl_eos_tabulated ||
       eos->table_type != ghl_eos_table_stellarcollapse ||
       !isfinite(rho) || rho <= 0.0 || !isfinite(T) || T <= 0.0 ||
       !isfinite(Ye) || Ye < 0.0 || Ye > 1.0)
      return ghl_error_m1_microphysics_failure;
    err = NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T(
        eos, rho, Ye, T, &muhat, &mu_e, &mu_p, &mu_n, &X_n, &X_p);
  }
  if(err != ghl_success)
    return err;
  return ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
      rho, Ye, T, muhat, mu_e, mu_p, mu_n, X_n, X_p,
      strict_validation, thermo);
#endif
}

ghl_error_codes_t ghl_m1_nrpyleakage_compute_thermo_state(
      const ghl_eos_parameters *restrict eos, const double rho,
      const double Ye, const double T,
      ghl_m1_nrpyleakage_thermo_state *restrict thermo) {
  return compute_thermo_state_impl(eos, rho, Ye, T, thermo, true, false);
}

ghl_error_codes_t ghl_m1_nrpyleakage_compute_thermo_state_legacy(
      const ghl_eos_parameters *restrict eos, const double rho,
      const double Ye, const double T,
      ghl_m1_nrpyleakage_thermo_state *restrict thermo) {
  return compute_thermo_state_impl(eos, rho, Ye, T, thermo, false, true);
}

#define EnsureFinite(value) \
  ghl_m1_nrpyleakage_finite_rate((value), strict_rates, &rate_failure)
#define EvaluateFermiDiracOrReturn(out, k, z)                              \
  do {                                                                     \
    const ghl_error_codes_t err =                                          \
        evaluate_fermi_dirac_integral((k), (z), &(out));                   \
    if(err != ghl_success)                                                 \
      return err;                                                          \
  } while(0)
#define EvaluateFermiFactorOrReturn(out, x)                                \
  do {                                                                     \
    const ghl_error_codes_t err = evaluate_fermi_factor((x), &(out));      \
    if(err != ghl_success)                                                 \
      return err;                                                          \
  } while(0)
static ghl_error_codes_t compute_kernel(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double eta[ghl_m1_nrpyleakage_species_count],
      ghl_m1_nrpyleakage_legacy_kernel_result *restrict result,
      const bool strict_rates) {
  if(thermo == NULL || eta == NULL || result == NULL)
    return ghl_error_m1_null_pointer;
  bool rate_failure = false;
  for(int s = 0; s < ghl_m1_nrpyleakage_species_count; ++s)
    if(!isfinite(eta[s]))
      return ghl_error_m1_microphysics_failure;

  const double T = thermo->T;
  const double mu_e = thermo->mu_e;
  const double mu_p = thermo->mu_p;
  const double mu_n = thermo->mu_n;
  const double rho_cgs = thermo->rho_cgs;
  const double Y_e = thermo->Ye;
  const double X_n = thermo->X_n;
  const double X_p = thermo->X_p;
  const double Y_np = thermo->Y_np;
  const double Y_pn = thermo->Y_pn;
  ghl_m1_nrpyleakage_legacy_kernel_result candidate = {0};
  candidate.raw.nux_single_species_multiplicity = 1;

  const double tmp_0 = (1.0/(T));
  const double tmp_1 = mu_e*tmp_0;
  double tmp_2_fd, tmp_3_fd, tmp_16_fd, tmp_18_fd, tmp_19_fd;
  EvaluateFermiDiracOrReturn(tmp_2_fd, 4, tmp_1);
  EvaluateFermiDiracOrReturn(tmp_3_fd, 5, tmp_1);
  EvaluateFermiDiracOrReturn(tmp_16_fd, 3, tmp_1);
  EvaluateFermiDiracOrReturn(tmp_18_fd, 4, -tmp_1);
  EvaluateFermiDiracOrReturn(tmp_19_fd, 3, -tmp_1);
  const double tmp_2 = tmp_2_fd;
  const double tmp_3 = tmp_3_fd/tmp_2;
  const double tmp_6 = eta[ghl_m1_nrpyleakage_nue];
  const double tmp_7 = ((NRPyLeakage_alpha)*(NRPyLeakage_alpha));
  const double tmp_8 = M_PI/NRPyLeakage_hc3;
  const double tmp_10 = 8*NRPyLeakage_N_A*NRPyLeakage_beta*((T)*(T)*(T)*(T)*(T))*rho_cgs*tmp_8*((3.0/8.0)*tmp_7 + 1.0/8.0);
  double tmp_11_factor;
  EvaluateFermiFactorOrReturn(tmp_11_factor, -tmp_3 + tmp_6);
  const double tmp_11 = NRPyLeakage_enable_beta_nue*EnsureFinite(
      Y_pn*tmp_10*tmp_2*tmp_11_factor);
  const double tmp_12 = NRPyLeakage_enable_brems_nui_anui*EnsureFinite(NRPyLeakage_Brems_C1*NRPyLeakage_Brems_zeta*pow(T, 4.5)*rho_cgs*(((X_n)*(X_n)) + (28.0/3.0)*X_n*X_p + ((X_p)*(X_p))));
  const double tmp_15 = eta[ghl_m1_nrpyleakage_anue];
  const double tmp_16 = tmp_16_fd;
  const double tmp_17 = (1.0/(tmp_16));
  const double tmp_18 = tmp_18_fd;
  const double tmp_19 = tmp_19_fd;
  const double tmp_20 = (1.0/(tmp_19));
  const double tmp_21 = -1.0/2.0*tmp_17*tmp_2 - 1.0/2.0*tmp_18*tmp_20;
  const double tmp_22 = ((M_PI)*(M_PI));
  const double tmp_24 = (1.0/((NRPyLeakage_hc3)*(NRPyLeakage_hc3)));
  const double tmp_25 = tmp_16*tmp_22*tmp_24;
  const double tmp_26 = pow(T, 8);
  const double tmp_28 = (16.0/9.0)*NRPyLeakage_beta*tmp_19*tmp_25*tmp_26;
  double tmp_29_factor_anue, tmp_29_factor_nue;
  EvaluateFermiFactorOrReturn(tmp_29_factor_anue, tmp_15 + tmp_21);
  EvaluateFermiFactorOrReturn(tmp_29_factor_nue, tmp_21 + tmp_6);
  const double tmp_29 = NRPyLeakage_enable_pair_nue_anue*EnsureFinite(
      NRPyLeakage_C1pC2_nue_anue*tmp_28*tmp_29_factor_anue*tmp_29_factor_nue);
  const double tmp_31 = (1.0/3.0)*tmp_22 + ((mu_e)*(mu_e))/((T)*(T));
  const double tmp_32 = NRPyLeakage_gamma_0*sqrt(tmp_31);
  const double tmp_34 = ((NRPyLeakage_gamma_0)*(NRPyLeakage_gamma_0))*tmp_31/(tmp_32 + 1);
  const double tmp_35 = -1.0/2.0*tmp_34 - 1;
  const double tmp_36 = (1.0/3.0)*((M_PI)*(M_PI)*(M_PI))*NRPyLeakage_beta*pow(NRPyLeakage_gamma_0, 6)*tmp_24*tmp_26*((tmp_31)*(tmp_31)*(tmp_31))*(tmp_32 + 1)*exp(-tmp_32)/NRPyLeakage_alpha_fs;
  double tmp_37_factor_anue, tmp_37_factor_nue;
  EvaluateFermiFactorOrReturn(tmp_37_factor_anue, tmp_15 + tmp_35);
  EvaluateFermiFactorOrReturn(tmp_37_factor_nue, tmp_35 + tmp_6);
  const double tmp_37 = NRPyLeakage_enable_plasmon_nue_anue*EnsureFinite(
      ((NRPyLeakage_C_V)*(NRPyLeakage_C_V))*tmp_36*tmp_37_factor_anue*
      tmp_37_factor_nue);
  const double tmp_38 = tmp_12 + tmp_29 + tmp_37;
  const double tmp_41 = (1 - Y_e)*((5.0/24.0)*tmp_7 + 1.0/24.0)/((2.0/3.0)*fmax(mu_n*tmp_0, 0) + 1);
  const double tmp_42 = NRPyLeakage_N_A*NRPyLeakage_sigma_0*((T)*(T))*rho_cgs/((NRPyLeakage_m_e_c2)*(NRPyLeakage_m_e_c2));
  double tmp_43_fd, tmp_44_fd, tmp_51_fd, tmp_57_fd, tmp_60_fd, tmp_61_fd, tmp_63_fd, tmp_70_fd, tmp_71_fd, tmp_78_fd, tmp_83_fd, tmp_87_fd4, tmp_87_fd2;
  EvaluateFermiDiracOrReturn(tmp_43_fd, 2, tmp_6);
  EvaluateFermiDiracOrReturn(tmp_44_fd, 4, tmp_6);
  EvaluateFermiDiracOrReturn(tmp_51_fd, 5, tmp_6);
  EvaluateFermiDiracOrReturn(tmp_57_fd, 5, -tmp_1);
  EvaluateFermiDiracOrReturn(tmp_60_fd, 2, tmp_15);
  EvaluateFermiDiracOrReturn(tmp_61_fd, 4, tmp_15);
  EvaluateFermiDiracOrReturn(tmp_63_fd, 5, tmp_15);
  EvaluateFermiDiracOrReturn(tmp_70_fd, 3, 0);
  EvaluateFermiDiracOrReturn(tmp_71_fd, 5, 0);
  EvaluateFermiDiracOrReturn(tmp_78_fd, 3, tmp_6);
  EvaluateFermiDiracOrReturn(tmp_83_fd, 3, tmp_15);
  EvaluateFermiDiracOrReturn(tmp_87_fd4, 4, 0);
  EvaluateFermiDiracOrReturn(tmp_87_fd2, 2, 0);
  const double tmp_43 = tmp_43_fd;
  const double tmp_44 = tmp_44_fd;
  const double tmp_45 = tmp_44/tmp_43;
  const double tmp_47 = ((NRPyLeakage_C_V - 1)*(NRPyLeakage_C_V - 1));
  const double tmp_48 = Y_e*((1.0/6.0)*tmp_47 + (5.0/24.0)*tmp_7)/((2.0/3.0)*fmax(mu_p*tmp_0, 0) + 1);
  const double tmp_49 = tmp_42*tmp_48;
  const double tmp_50 = (3.0/4.0)*tmp_7 + 1.0/4.0;
  const double tmp_51 = tmp_51_fd;
  double tmp_52_factor;
  EvaluateFermiFactorOrReturn(tmp_52_factor, tmp_1 - tmp_51/tmp_44);
  const double tmp_52 = Y_np*tmp_50*tmp_52_factor;
  const double nue_N_p = EnsureFinite(tmp_45*tmp_49);
  const double nue_N_n = EnsureFinite(tmp_41*tmp_42*tmp_45);
  const double nue_N_cc = EnsureFinite(tmp_42*tmp_45*tmp_52);
  const double tmp_53 = nue_N_p + nue_N_n + nue_N_cc;
  const double tmp_55 = 4*((T)*(T)*(T))*tmp_8;
  const double tmp_57 = tmp_57_fd/tmp_18;
  double tmp_58_factor;
  EvaluateFermiFactorOrReturn(tmp_58_factor, tmp_15 - tmp_57);
  const double tmp_58 = NRPyLeakage_enable_beta_anue*EnsureFinite(
      Y_np*tmp_10*tmp_18*tmp_58_factor);
  const double tmp_60 = tmp_60_fd;
  const double tmp_61 = tmp_61_fd;
  const double tmp_62 = tmp_42*tmp_61/tmp_60;
  const double tmp_63 = tmp_63_fd;
  double tmp_64_factor;
  EvaluateFermiFactorOrReturn(tmp_64_factor, -tmp_1 - tmp_63/tmp_61);
  const double tmp_64 = Y_pn*tmp_50*tmp_64_factor;
  const double anue_N_n = EnsureFinite(tmp_41*tmp_62);
  const double anue_N_p = EnsureFinite(tmp_48*tmp_62);
  const double anue_N_cc = EnsureFinite(tmp_62*tmp_64);
  const double tmp_65 = anue_N_n + anue_N_p + anue_N_cc;
  const double tmp_66 = EnsureFinite(NRPyLeakage_Brems_C2*T*tmp_12/NRPyLeakage_Brems_C1);
  const double tmp_67 = 32*pow(T, 9);
  const double tmp_68 = (1.0/64.0)*((NRPyLeakage_hc3)*(NRPyLeakage_hc3))*tmp_17*tmp_20*(tmp_18*tmp_25*tmp_67 + tmp_19*tmp_2*tmp_22*tmp_24*tmp_67)/(tmp_22*tmp_26);
  const double tmp_69 = (1.0/2.0)*T*(tmp_34 + 2);
  const double tmp_70 = tmp_70_fd;
  const double tmp_71 = tmp_71_fd/tmp_70;
  const double nux_E_p = EnsureFinite(tmp_49*tmp_71);
  const double nux_E_n = EnsureFinite(tmp_41*tmp_42*tmp_71);
  const double tmp_73 = nux_E_p + nux_E_n;
  const double tmp_74 = 4*((T)*(T)*(T)*(T))*tmp_8;
  const double pair_E = EnsureFinite(tmp_29*tmp_68);
  const double plasmon_E = EnsureFinite(tmp_37*tmp_69);
  const double tmp_75 = tmp_66 + pair_E + plasmon_E;
  const double beta_E_nue = EnsureFinite(T*tmp_11*tmp_3);
  const double tmp_76 = tmp_75 + beta_E_nue;
  const double tmp_78 = tmp_78_fd;
  const double tmp_79 = tmp_51/tmp_78;
  const double nue_E_p = EnsureFinite(tmp_49*tmp_79);
  const double nue_E_n = EnsureFinite(tmp_41*tmp_42*tmp_79);
  const double nue_E_cc = EnsureFinite(tmp_42*tmp_52*tmp_79);
  const double tmp_81 = nue_E_p + nue_E_n + nue_E_cc;
  const double beta_E_anue = EnsureFinite(T*tmp_57*tmp_58);
  const double tmp_82 = tmp_75 + beta_E_anue;
  const double tmp_83 = tmp_83_fd;
  const double tmp_84 = tmp_63/tmp_83;
  const double anue_E_p = EnsureFinite(tmp_49*tmp_84);
  const double anue_E_n = EnsureFinite(tmp_41*tmp_42*tmp_84);
  const double anue_E_cc = EnsureFinite(tmp_42*tmp_64*tmp_84);
  const double tmp_86 = anue_E_p + anue_E_n + anue_E_cc;
  const double tmp_87 = tmp_87_fd4/tmp_87_fd2;
  const double nux_N_p = EnsureFinite(tmp_49*tmp_87);
  const double nux_N_n = EnsureFinite(tmp_41*tmp_42*tmp_87);
  const double nux_N = nux_N_p + nux_N_n;
  double pair_N_nux_factor, plasmon_N_nux_factor;
  EvaluateFermiFactorOrReturn(pair_N_nux_factor, tmp_21);
  EvaluateFermiFactorOrReturn(plasmon_N_nux_factor, tmp_35);
  const double pair_N_nux = NRPyLeakage_enable_pair_nux_anux*EnsureFinite(
      NRPyLeakage_C1pC2_nux_anux*tmp_28*pair_N_nux_factor*pair_N_nux_factor);
  const double plasmon_N_nux = NRPyLeakage_enable_plasmon_nux_anux*EnsureFinite(
      tmp_36*tmp_47*plasmon_N_nux_factor*plasmon_N_nux_factor);
  const double pair_E_nux = EnsureFinite(tmp_68*pair_N_nux);
  const double plasmon_E_nux = EnsureFinite(tmp_69*plasmon_N_nux);
  const double eta_E_nux = tmp_66 + pair_E_nux + plasmon_E_nux;

  ghl_m1_nrpyleakage_species_raw_rates *const nue = &candidate.raw.species[0];
  ghl_m1_nrpyleakage_species_raw_rates *const anue = &candidate.raw.species[1];
  ghl_m1_nrpyleakage_species_raw_rates *const nux = &candidate.raw.species[2];
  const double F2[3] = {tmp_43, tmp_60, tmp_87_fd2};
  const double F3[3] = {tmp_78, tmp_83, tmp_70};
  const double F4[3] = {tmp_44, tmp_61, tmp_87_fd4};
  const double F5[3] = {tmp_51, tmp_63, tmp_71_fd};
  for(int s = 0; s < 3; ++s) {
    ghl_m1_nrpyleakage_species_raw_rates *const r = &candidate.raw.species[s];
    r->neutrino_degeneracy = eta[s];
    r->F2 = F2[s]; r->F3 = F3[s]; r->F4 = F4[s]; r->F5 = F5[s];
    r->n_eq_cgs = 4*M_PI*((T)*(T)*(T))*F2[s]/NRPyLeakage_hc3;
    r->J_eq_mev_cgs = 4*M_PI*((T)*(T)*(T)*(T))*F3[s]/NRPyLeakage_hc3;
    r->mean_energy_mev = T*F3[s]/F2[s];
    r->eta_N_brems_cgs = tmp_12;
    r->eta_E_brems_mev_cgs = tmp_66;
  }
  nue->eta_N_beta_cgs = tmp_11; anue->eta_N_beta_cgs = tmp_58;
  nue->eta_E_beta_mev_cgs = beta_E_nue;
  anue->eta_E_beta_mev_cgs = beta_E_anue;
  nue->eta_N_pair_cgs = anue->eta_N_pair_cgs = tmp_29;
  nue->eta_E_pair_mev_cgs = anue->eta_E_pair_mev_cgs = pair_E;
  nue->eta_N_plasmon_cgs = anue->eta_N_plasmon_cgs = tmp_37;
  nue->eta_E_plasmon_mev_cgs = anue->eta_E_plasmon_mev_cgs = plasmon_E;
  nux->eta_N_pair_cgs = pair_N_nux; nux->eta_E_pair_mev_cgs = pair_E_nux;
  nux->eta_N_plasmon_cgs = plasmon_N_nux;
  nux->eta_E_plasmon_mev_cgs = plasmon_E_nux;
  nue->kappa_s_N_neutron_cgs = nue_N_n; nue->kappa_s_N_proton_cgs = nue_N_p;
  nue->kappa_a_N_cc_cgs = nue_N_cc; nue->kappa_s_E_neutron_cgs = nue_E_n;
  nue->kappa_s_E_proton_cgs = nue_E_p; nue->kappa_a_E_cc_cgs = nue_E_cc;
  anue->kappa_s_N_neutron_cgs = anue_N_n; anue->kappa_s_N_proton_cgs = anue_N_p;
  anue->kappa_a_N_cc_cgs = anue_N_cc; anue->kappa_s_E_neutron_cgs = anue_E_n;
  anue->kappa_s_E_proton_cgs = anue_E_p; anue->kappa_a_E_cc_cgs = anue_E_cc;
  nux->kappa_s_N_neutron_cgs = nux_N_n; nux->kappa_s_N_proton_cgs = nux_N_p;
  nux->kappa_s_E_neutron_cgs = nux_E_n; nux->kappa_s_E_proton_cgs = nux_E_p;

  candidate.eta_N[0] = tmp_11 + tmp_38;
  candidate.eta_N[1] = tmp_38 + tmp_58;
  candidate.eta_N[2] = tmp_12 + pair_N_nux + plasmon_N_nux;
  candidate.eta_E[0] = tmp_76; candidate.eta_E[1] = tmp_82;
  candidate.eta_E[2] = eta_E_nux;
  candidate.kappa_N[0] = tmp_53; candidate.kappa_N[1] = tmp_65;
  candidate.kappa_N[2] = nux_N;
  candidate.kappa_E[0] = tmp_81; candidate.kappa_E[1] = tmp_86;
  candidate.kappa_E[2] = tmp_73;
  candidate.n_eq[0] = tmp_43*tmp_55;
  candidate.n_eq[1] = tmp_55*tmp_60;
  candidate.n_eq[2] = tmp_87_fd2*tmp_55;
  candidate.J_eq[0] = tmp_74*tmp_78;
  candidate.J_eq[1] = tmp_74*tmp_83;
  candidate.J_eq[2] = tmp_70*tmp_74;
  if(strict_rates) {
    if(rate_failure)
      return ghl_error_m1_microphysics_failure;
    const ghl_error_codes_t err = validate_kernel_result(&candidate);
    if(err != ghl_success)
      return err;
  }
  *result = candidate;
  return ghl_success;
}
#undef EnsureFinite
#undef EvaluateFermiDiracOrReturn
#undef EvaluateFermiFactorOrReturn

ghl_error_codes_t ghl_m1_nrpyleakage_compute_legacy_kernel_from_thermo(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double eta[ghl_m1_nrpyleakage_species_count],
      ghl_m1_nrpyleakage_legacy_kernel_result *restrict result) {
  return compute_kernel(thermo, eta, result, false);
}

ghl_error_codes_t ghl_m1_nrpyleakage_compute_raw_rates_from_thermo_legacy(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double eta[ghl_m1_nrpyleakage_species_count],
      ghl_m1_nrpyleakage_raw_rates *restrict raw) {
  if(raw == NULL) return ghl_error_m1_null_pointer;
  ghl_m1_nrpyleakage_legacy_kernel_result result;
  const ghl_error_codes_t err = compute_kernel(thermo, eta, &result, false);
  if(err == ghl_success) *raw = result.raw;
  return err;
}

ghl_error_codes_t ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(
      const ghl_m1_nrpyleakage_thermo_state *restrict thermo,
      const double eta[ghl_m1_nrpyleakage_species_count],
      ghl_m1_nrpyleakage_raw_rates *restrict raw) {
  if(thermo == NULL || eta == NULL || raw == NULL)
    return ghl_error_m1_null_pointer;
  const double state_values[] = {
    thermo->rho, thermo->rho_cgs, thermo->T, thermo->Ye, thermo->muhat,
    thermo->mu_e, thermo->mu_p, thermo->mu_n, thermo->X_n, thermo->X_p,
    thermo->Y_np, thermo->Y_pn
  };
  for(size_t i = 0; i < sizeof(state_values)/sizeof(state_values[0]); ++i)
    if(!isfinite(state_values[i])) return ghl_error_m1_microphysics_failure;
  if(thermo->rho <= 0.0 || thermo->rho_cgs <= 0.0 || thermo->T <= 0.0 ||
     thermo->Ye < 0.0 || thermo->Ye > 1.0 || thermo->X_n < 0.0 ||
     thermo->X_p < 0.0 || thermo->Y_np < 0.0 || thermo->Y_pn < 0.0)
    return ghl_error_m1_microphysics_failure;
  ghl_m1_nrpyleakage_legacy_kernel_result result;
  ghl_error_codes_t err = compute_kernel(thermo, eta, &result, true);
  if(err != ghl_success) return err;
  for(int s = 0; s < ghl_m1_nrpyleakage_species_count; ++s) {
    err = validate_raw_species(&result.raw.species[s]);
    if(err != ghl_success) return err;
  }
  *raw = result.raw;
  return ghl_success;
}
