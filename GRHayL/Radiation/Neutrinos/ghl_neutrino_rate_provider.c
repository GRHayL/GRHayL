#include "ghl_radiation.h"
#include "ghl_m1_nrpyleakage_kernel.h"
#ifndef GHL_DISABLE_HDF5
#include "ghl_nrpyeos_tabulated.h"
#endif
#include <float.h>
#include <string.h>

/*
 * Shared orchestration for deterministic reference and production
 * NRPyLeakage backends of the frozen aggregate-rate interface.
 *
 * The default initializer preserves the deterministic reference backend. The
 * explicit production initializer consumes the private, tau-free
 * NRPyLeakage-derived Ruffert channel kernel and constructs effective
 * aggregate grey coefficients by detailed balance. Both return an
 * already-summed nu_x bundle and accept only nu_x_multiplicity == 4.
 */

#ifndef CODE_TO_CGS_DENSITY
#define CODE_TO_CGS_DENSITY  6.17714470405638e+17
#endif

#ifndef CGS_TO_CODE_LENGTH
#define CGS_TO_CODE_LENGTH   6.77269222552442e-06
#endif

#define MEV_TO_ERG 1.602176634e-6

/* Detailed-balance aggregation divides channel emissivities by the
 * equilibrium number and energy targets.  Preserve a finite, positive target
 * when a physical tail underflows before reaching that division.  This is far
 * below any resolved production rate and is only active after the raw value
 * has rounded to zero. */
#define GHL_M1_NRPYLEAKAGE_EQUILIBRIUM_FLOOR 1.0e-300

static double clamp_double(const double x, const double lo, const double hi) {
  return fmin(hi, fmax(lo, x));
}

static double safe_exp(const double x) {
  return exp(clamp_double(x, -40.0, 40.0));
}

static double species_lepton_weight(const ghl_m1_neutrino_species_t species) {
  switch(species) {
    case ghl_m1_neutrino_nue:  return  1.0;
    case ghl_m1_neutrino_anue: return -1.0;
    case ghl_m1_neutrino_nux:  return  0.0;
    default:                   return  0.0;
  }
}

static void fill_transparent_rates(
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count],
      const double mean_energy) {

  const double safe_mean = fmax(mean_energy, 1.0);
  for(int s = 0; s < ghl_m1_neutrino_species_count; s++) {
    rates[s] = (ghl_m1_neutrino_rates){0};
    rates[s].species = (ghl_m1_neutrino_species_t)s;
    rates[s].mean_energy = safe_mean;
    rates[s].lepton_weight = species_lepton_weight(rates[s].species);
  }
}

static ghl_error_codes_t publish_recovered_rates(
      const ghl_m1_neutrino_rates candidate[ghl_m1_neutrino_species_count],
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {

  for(int s = 0; s < ghl_m1_neutrino_species_count; ++s) {
    const ghl_error_codes_t error =
        ghl_m1_validate_neutrino_rates(&candidate[s], NULL);
    if(error != ghl_success)
      return error;
  }
  memcpy(rates, candidate,
         sizeof(ghl_m1_neutrino_rates) * ghl_m1_neutrino_species_count);
  return ghl_success;
}

ghl_error_codes_t ghl_neutrino_rate_provider_initialize_default(
      ghl_neutrino_rate_provider_context *restrict provider) {

  if(provider == NULL)
    return ghl_error_m1_null_pointer;

  *provider = (ghl_neutrino_rate_provider_context){0};
  /* The default provider is deliberately table-free in every build mode.
   * Tabulated NRPyLeakage rates are selected only by the explicit
   * ghl_neutrino_rate_provider_initialize_nrpyleakage() initializer. */
  provider->use_tabulated_eos = false;
  provider->channel_mask = ghl_neutrino_rate_channel_charged_current |
                           ghl_neutrino_rate_channel_nucleon_scattering;
  provider->failure_policy = ghl_neutrino_rate_failure_abort;
  provider->table_bounds_policy = ghl_neutrino_rate_table_bounds_abort;
  provider->nu_x_multiplicity = 4.0;
  provider->eos_generation = 0;
  provider->rho_code_to_cgs = 1.0;
  provider->temperature_code_to_mev = 1.0;
  provider->opacity_cgs_to_code = 1.0;
  provider->emissivity_cgs_to_code = 1.0;
  provider->baryon_mass_code = 1.0;
  provider->charged_current_scale = 1.0e-2;
  provider->scattering_scale = 5.0e-3;
  provider->pair_scale = 1.0e-4;
  provider->bremsstrahlung_scale = 1.0e-4;
  provider->plasmon_scale = 1.0e-5;
  provider->min_mean_energy = 1.0e-12;
  provider->backend = ghl_neutrino_rate_backend_reference;
  provider->equilibrium_recovery_rate = 1.0;
  return ghl_success;
}

void ghl_neutrino_rate_provider_cache_initialize(
      ghl_neutrino_rate_provider_cache *restrict cache) {

  if(cache == NULL)
    return;

  *cache = (ghl_neutrino_rate_provider_cache){0};
  /* Keep the cache validity contract explicit even though the compound
   * literal above already initializes both flags to false. */
  cache->thermo_valid = false;
  cache->rates_valid = false;
}

ghl_error_codes_t ghl_neutrino_rate_provider_initialize_nrpyleakage(
      ghl_neutrino_rate_provider_context *restrict provider) {

  if(provider == NULL)
    return ghl_error_m1_null_pointer;
#ifdef GHL_DISABLE_HDF5
  return ghl_error_used_disabled_hdf5;
#else
  ghl_error_codes_t err = ghl_neutrino_rate_provider_initialize_default(provider);
  if(err != ghl_success)
    return err;
  provider->backend = ghl_neutrino_rate_backend_nrpyleakage;
  provider->use_tabulated_eos = true;
  provider->channel_mask = ghl_neutrino_rate_channel_charged_current |
                           ghl_neutrino_rate_channel_nucleon_scattering |
                           ghl_neutrino_rate_channel_pair |
                           ghl_neutrino_rate_channel_bremsstrahlung |
                           ghl_neutrino_rate_channel_plasmon;
  provider->failure_policy = ghl_neutrino_rate_failure_abort;
  provider->table_bounds_policy = ghl_neutrino_rate_table_bounds_abort;
  provider->nu_x_multiplicity = 4.0;
  provider->rho_code_to_cgs = NRPyLeakage_units_geom_to_cgs_D;
  provider->temperature_code_to_mev = 1.0;
  provider->opacity_cgs_to_code = NRPyLeakage_units_geom_to_cgs_L;
  provider->emissivity_cgs_to_code =
      pow(NRPyLeakage_units_geom_to_cgs_L, 3)*
      NRPyLeakage_units_geom_to_cgs_T;
  provider->baryon_mass_code =
      NRPyLeakage_amu/NRPyLeakage_units_geom_to_cgs_M;
  provider->charged_current_scale = 1.0;
  provider->scattering_scale = 1.0;
  provider->pair_scale = 1.0;
  provider->bremsstrahlung_scale = 1.0;
  provider->plasmon_scale = 1.0;
  provider->min_mean_energy = MEV_TO_ERG/
      (NRPyLeakage_units_geom_to_cgs_M*NRPyLeakage_c_light*
       NRPyLeakage_c_light);
  return ghl_success;
#endif
}

static ghl_error_codes_t validate_provider_context(
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_eos_parameters *restrict eos) {

  const int valid_channels = ghl_neutrino_rate_channel_charged_current |
      ghl_neutrino_rate_channel_nucleon_scattering |
      ghl_neutrino_rate_channel_pair |
      ghl_neutrino_rate_channel_bremsstrahlung |
      ghl_neutrino_rate_channel_plasmon;
  const double nonnegative[] = {
    provider->charged_current_scale, provider->scattering_scale,
    provider->pair_scale, provider->bremsstrahlung_scale,
    provider->plasmon_scale
  };
  const double positive[] = {
    provider->temperature_code_to_mev, provider->baryon_mass_code,
    provider->min_mean_energy, provider->equilibrium_recovery_rate
  };

  if(provider->backend < ghl_neutrino_rate_backend_reference ||
     provider->backend > ghl_neutrino_rate_backend_nrpyleakage ||
     !isfinite(provider->nu_x_multiplicity) ||
     provider->nu_x_multiplicity != 4.0 ||
     (provider->channel_mask & ~valid_channels) != 0 ||
     provider->failure_policy < ghl_neutrino_rate_failure_abort ||
     provider->failure_policy > ghl_neutrino_rate_failure_hold_last ||
     provider->table_bounds_policy < ghl_neutrino_rate_table_bounds_abort ||
     provider->table_bounds_policy > ghl_neutrino_rate_table_bounds_clamp)
    return ghl_error_m1_microphysics_failure;
  for(size_t i = 0; i < sizeof(nonnegative)/sizeof(nonnegative[0]); ++i)
    if(!isfinite(nonnegative[i]) || nonnegative[i] < 0.0)
      return ghl_error_m1_microphysics_failure;
  for(size_t i = 0; i < sizeof(positive)/sizeof(positive[0]); ++i)
    if(!isfinite(positive[i]) || positive[i] <= 0.0)
      return ghl_error_m1_microphysics_failure;
#ifdef GHL_DISABLE_HDF5
  if(provider->backend == ghl_neutrino_rate_backend_nrpyleakage)
    return ghl_error_used_disabled_hdf5;
#endif
  if(provider->backend == ghl_neutrino_rate_backend_reference) {
    if(provider->rho_code_to_cgs != 1.0 ||
       provider->opacity_cgs_to_code != 1.0 ||
       provider->emissivity_cgs_to_code != 1.0)
      return ghl_error_m1_microphysics_failure;
  } else {
    const double nrpy_min_mean_energy = MEV_TO_ERG/
        (NRPyLeakage_units_geom_to_cgs_M*NRPyLeakage_c_light*
         NRPyLeakage_c_light);
    if(!provider->use_tabulated_eos || eos == NULL ||
       eos->eos_type != ghl_eos_tabulated ||
       eos->table_type != ghl_eos_table_stellarcollapse ||
       provider->rho_code_to_cgs != NRPyLeakage_units_geom_to_cgs_D ||
       provider->temperature_code_to_mev != 1.0 ||
       provider->opacity_cgs_to_code != NRPyLeakage_units_geom_to_cgs_L ||
       provider->emissivity_cgs_to_code !=
           pow(NRPyLeakage_units_geom_to_cgs_L, 3)*
           NRPyLeakage_units_geom_to_cgs_T ||
       provider->baryon_mass_code !=
           NRPyLeakage_amu/NRPyLeakage_units_geom_to_cgs_M ||
       provider->charged_current_scale != 1.0 ||
       provider->scattering_scale != 1.0 ||
       provider->pair_scale != 1.0 ||
       provider->bremsstrahlung_scale != 1.0 ||
       provider->plasmon_scale != 1.0 ||
       provider->min_mean_energy != nrpy_min_mean_energy)
      return ghl_error_m1_microphysics_failure;
  }
  if(provider->use_tabulated_eos && eos == NULL)
    return ghl_error_m1_microphysics_failure;
#ifdef GHL_DISABLE_HDF5
  if(provider->use_tabulated_eos)
    return ghl_error_m1_microphysics_failure;
#endif
  return ghl_success;
}

static bool same_provider_configuration(
      const ghl_neutrino_rate_provider_context *restrict lhs,
      const ghl_neutrino_rate_provider_context *restrict rhs) {

  return lhs->use_tabulated_eos == rhs->use_tabulated_eos &&
      lhs->channel_mask == rhs->channel_mask &&
      lhs->failure_policy == rhs->failure_policy &&
      lhs->table_bounds_policy == rhs->table_bounds_policy &&
      lhs->nu_x_multiplicity == rhs->nu_x_multiplicity &&
      lhs->eos_generation == rhs->eos_generation &&
      lhs->rho_code_to_cgs == rhs->rho_code_to_cgs &&
      lhs->temperature_code_to_mev == rhs->temperature_code_to_mev &&
      lhs->opacity_cgs_to_code == rhs->opacity_cgs_to_code &&
      lhs->emissivity_cgs_to_code == rhs->emissivity_cgs_to_code &&
      lhs->baryon_mass_code == rhs->baryon_mass_code &&
      lhs->charged_current_scale == rhs->charged_current_scale &&
      lhs->scattering_scale == rhs->scattering_scale &&
      lhs->pair_scale == rhs->pair_scale &&
      lhs->bremsstrahlung_scale == rhs->bremsstrahlung_scale &&
      lhs->plasmon_scale == rhs->plasmon_scale &&
      lhs->min_mean_energy == rhs->min_mean_energy &&
      lhs->backend == rhs->backend &&
      lhs->equilibrium_recovery_rate == rhs->equilibrium_recovery_rate;
}

static bool same_provenance(
      const ghl_neutrino_rate_provider_cache *restrict cache,
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_eos_parameters *restrict eos) {

  return cache != NULL && cache->eos_snapshot == eos &&
      same_provider_configuration(&cache->provider_snapshot, provider);
}

static bool same_thermo_key(
      const ghl_neutrino_rate_provider_cache *restrict cache,
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_eos_parameters *restrict eos,
      const double rho, const double T, const double Ye) {

  return same_provenance(cache, provider, eos) && cache->thermo_valid &&
      cache->thermo_rho == rho && cache->thermo_T == T &&
      cache->thermo_Ye == Ye;
}

static bool same_rates_key(
      const ghl_neutrino_rate_provider_cache *restrict cache,
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_eos_parameters *restrict eos,
      const double rho, const double T, const double Ye) {

  return same_provenance(cache, provider, eos) && cache->rates_valid &&
      cache->rho == rho && cache->T == T && cache->Ye == Ye;
}

static ghl_error_codes_t validate_inputs(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict T,
      double *restrict Ye) {

  if(!isfinite(*rho) || !isfinite(*T) || !isfinite(*Ye) ||
     *rho <= 0.0 || *T <= 0.0 || *Ye < 0.0 || *Ye > 1.0)
    return ghl_error_m1_microphysics_failure;

  if(!provider->use_tabulated_eos || eos == NULL)
    return ghl_success;

  bool out_of_bounds = false;
  const double rho_lo = eos->table_rho_min > 0.0 ? eos->table_rho_min : eos->rho_min;
  const double rho_hi = eos->table_rho_max > 0.0 ? eos->table_rho_max : eos->rho_max;
  const double T_lo   = eos->table_T_min   > 0.0 ? eos->table_T_min   : eos->T_min;
  const double T_hi   = eos->table_T_max   > 0.0 ? eos->table_T_max   : eos->T_max;
  const double Ye_lo  = eos->table_Y_e_min > 0.0 ? eos->table_Y_e_min : eos->Y_e_min;
  const double Ye_hi  = eos->table_Y_e_max > 0.0 ? eos->table_Y_e_max : eos->Y_e_max;

  if(!isfinite(rho_lo) || !isfinite(rho_hi) || rho_lo <= 0.0 || rho_lo >= rho_hi ||
     !isfinite(T_lo) || !isfinite(T_hi) || T_lo <= 0.0 || T_lo >= T_hi ||
     !isfinite(Ye_lo) || !isfinite(Ye_hi) || Ye_lo < 0.0 || Ye_hi > 1.0 ||
     Ye_lo >= Ye_hi)
    return ghl_error_m1_microphysics_failure;

  if(*rho < rho_lo || *rho > rho_hi)
    out_of_bounds = true;
  if(*T < T_lo || *T > T_hi)
    out_of_bounds = true;
  if(*Ye < Ye_lo || *Ye > Ye_hi)
    out_of_bounds = true;

  if(!out_of_bounds)
    return ghl_success;

  if(diagnostics != NULL)
    diagnostics->table_bound_hits++;

  if(provider->table_bounds_policy != ghl_neutrino_rate_table_bounds_clamp)
    return ghl_error_m1_microphysics_failure;

  *rho = clamp_double(*rho, rho_lo, rho_hi);
  *T = clamp_double(*T, T_lo, T_hi);
  *Ye = clamp_double(*Ye, Ye_lo, Ye_hi);
  if(diagnostics != NULL)
    diagnostics->clamped_inputs++;
  return ghl_success;
}

static ghl_error_codes_t compute_thermo(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_cache *restrict cache,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims,
      double *restrict rho,
      double *restrict T,
      double *restrict Ye,
      double *restrict muhat,
      double *restrict mu_e,
      double *restrict mu_p,
      double *restrict mu_n,
      double *restrict X_n,
      double *restrict X_p,
      bool *restrict recovered_key_valid) {

  *recovered_key_valid = false;
  *rho = prims->rho;
  *Ye = prims->Y_e;
  *T = prims->temperature;
  if(provider->use_tabulated_eos) {
#ifdef GHL_DISABLE_HDF5
    return ghl_error_m1_microphysics_failure;
#else
    if(!isfinite(*T) || *T <= 0.0) {
      const double T_lo = eos->table_T_min > 0.0 ? eos->table_T_min : eos->T_min;
      const double T_hi = eos->table_T_max > 0.0 ? eos->table_T_max : eos->T_max;
      if(!isfinite(T_lo) || !isfinite(T_hi) || T_lo <= 0.0 || T_lo >= T_hi ||
         !isfinite(prims->eps))
        return ghl_error_m1_microphysics_failure;
      *T = exp(0.5 * (log(T_lo) + log(T_hi)));
      /* Bounds are applied before the inverse EOS call below. */
      ghl_error_codes_t err = validate_inputs(provider, diagnostics, eos, rho, T, Ye);
      if(err != ghl_success)
        return err;
      err = NRPyEOS_T_from_rho_Ye_eps(eos, *rho, *Ye, prims->eps, T);
      if(err != ghl_success)
        return err;
    }
#endif
  } else if(!isfinite(*T) || *T <= 0.0) {
    if(!isfinite(prims->eps) || prims->eps <= 0.0)
      return ghl_error_m1_microphysics_failure;
    *T = prims->eps;
  }

  ghl_error_codes_t err = validate_inputs(provider, diagnostics, eos, rho, T, Ye);
  if(err != ghl_success)
    return err;
  *recovered_key_valid = true;

  if(same_thermo_key(cache, provider, eos, *rho, *T, *Ye)) {
    *muhat = cache->muhat;
    *mu_e = cache->mu_e;
    *mu_p = cache->mu_p;
    *mu_n = cache->mu_n;
    *X_n = cache->X_n;
    *X_p = cache->X_p;
    return ghl_success;
  }

  if(provider->use_tabulated_eos) {
#ifdef GHL_DISABLE_HDF5
    return ghl_error_m1_microphysics_failure;
#else
    if(eos == NULL)
      return ghl_error_m1_microphysics_failure;
    if(provider->backend == ghl_neutrino_rate_backend_reference) {
      double P = 0.0, eps = 0.0;
      err = NRPyEOS_P_eps_muhat_mue_mup_and_mun_from_rho_Ye_T(
          eos, *rho, *Ye, *T, &P, &eps, muhat, mu_e, mu_p, mu_n);
      if(err != ghl_success)
        return err;
    }
    err = NRPyEOS_muhat_mue_mup_mun_Xn_and_Xp_from_rho_Ye_T(
        eos, *rho, *Ye, *T, muhat, mu_e, mu_p, mu_n, X_n, X_p);
    if(err != ghl_success)
      return err;
#endif
  } else {
    const double ye = clamp_double(*Ye, 1.0e-12, 1.0 - 1.0e-12);
    *X_p = ye;
    *X_n = 1.0 - ye;
    *mu_e = (*T) * log(ye / (1.0 - ye));
    *mu_p = (*T) * log(fmax(*X_p, 1.0e-12));
    *mu_n = (*T) * log(fmax(*X_n, 1.0e-12));
    *muhat = *mu_n - *mu_p;
  }

  if(!isfinite(*muhat) || !isfinite(*mu_e) || !isfinite(*mu_p) ||
     !isfinite(*mu_n) || !isfinite(*X_n) || !isfinite(*X_p) ||
     *X_n < 0.0 || *X_p < 0.0)
    return ghl_error_m1_microphysics_failure;

  if(cache != NULL) {
    cache->thermo_valid = true;
    cache->thermo_rho = *rho;
    cache->thermo_T = *T;
    cache->thermo_Ye = *Ye;
    cache->muhat = *muhat;
    cache->mu_e = *mu_e;
    cache->mu_p = *mu_p;
    cache->mu_n = *mu_n;
    cache->X_n = *X_n;
    cache->X_p = *X_p;
  }
  return ghl_success;
}

static ghl_error_codes_t assemble_nrpyleakage_rates(
      const ghl_neutrino_rate_provider_context *restrict provider,
      const double rho,
      const double T,
      const double Ye,
      const double muhat,
      const double mu_e,
      const double mu_p,
      const double mu_n,
      const double X_n,
      const double X_p,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count],
      double beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_species_count],
      bool beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_species_count]) {

  ghl_m1_nrpyleakage_thermo_state thermo;
  ghl_error_codes_t err =
      ghl_m1_nrpyleakage_build_thermo_state_from_eos_quantities(
          rho, Ye, T, muhat, mu_e, mu_p, mu_n, X_n, X_p, true, &thermo);
  if(err != ghl_success)
    return err;

  const double eta[ghl_m1_nrpyleakage_species_count] = {
    (mu_e - muhat)/T, -(mu_e - muhat)/T, 0.0
  };
  ghl_m1_nrpyleakage_raw_rates raw;
  err = ghl_m1_nrpyleakage_compute_raw_rates_from_thermo(&thermo, eta, &raw);
  if(err != ghl_success)
    return err;
  /* The raw bridge exposes one unsummed heavy species.  The provider owns
   * the sole conversion to the configured four-flavor aggregate below. */
  if(raw.nux_single_species_multiplicity != 1)
    return ghl_error_m1_microphysics_failure;

  const double L0 = NRPyLeakage_units_geom_to_cgs_L;
  const double t0 = NRPyLeakage_units_geom_to_cgs_T;
  const double E0 = NRPyLeakage_units_geom_to_cgs_M*
                    NRPyLeakage_c_light*NRPyLeakage_c_light;
  const double volume = L0*L0*L0;
  const double number_emissivity_conversion = volume*t0;
  const double energy_emissivity_conversion =
      MEV_TO_ERG*volume*t0/E0;
  const int mask = provider->channel_mask;

  for(int s = 0; s < ghl_m1_neutrino_species_count; ++s) {
    const ghl_m1_nrpyleakage_species_raw_rates *const rr = &raw.species[s];
    ghl_m1_neutrino_rates *const out = &rates[s];
    const double multiplicity = s == ghl_m1_neutrino_nux
                              ? provider->nu_x_multiplicity : 1.0;
    beta_kirchhoff_relative_mismatch[s] = 0.0;
    beta_kirchhoff_mismatch_valid[s] = false;
    *out = (ghl_m1_neutrino_rates){0};
    out->species = (ghl_m1_neutrino_species_t)s;
    out->lepton_weight = species_lepton_weight(out->species);
    const double raw_n_eq = multiplicity*rr->n_eq_cgs*volume;
    const double raw_mean_energy = rr->mean_energy_mev*MEV_TO_ERG/E0;
    if(!isfinite(raw_n_eq) || raw_n_eq < 0.0 ||
       !isfinite(raw_mean_energy) || raw_mean_energy <= 0.0)
      return ghl_error_m1_microphysics_failure;
    out->n_eq = fmax(raw_n_eq,
                     GHL_M1_NRPYLEAKAGE_EQUILIBRIUM_FLOOR);
    out->mean_energy = fmax(raw_mean_energy, provider->min_mean_energy);
    out->J_eq = out->n_eq*out->mean_energy;
    if(!isfinite(out->J_eq) || out->J_eq <= 0.0)
      return ghl_error_m1_microphysics_failure;

    if(s != ghl_m1_neutrino_nux) {
      const double beta = rr->eta_N_beta_cgs*number_emissivity_conversion;
      const double kirchhoff = rr->kappa_a_N_cc_cgs*L0*out->n_eq;
      beta_kirchhoff_relative_mismatch[s] = fabs(beta - kirchhoff)/
          fmax(fmax(fabs(beta), fabs(kirchhoff)), DBL_MIN);
      beta_kirchhoff_mismatch_valid[s] = true;
    }

    double kappa_a_N = 0.0;
    double kappa_a_E = 0.0;
    if((mask & ghl_neutrino_rate_channel_charged_current) != 0 &&
       s != ghl_m1_neutrino_nux) {
      out->kappa_a_N_cc = rr->kappa_a_N_cc_cgs*L0;
      kappa_a_N += out->kappa_a_N_cc;
      kappa_a_E += rr->kappa_a_E_cc_cgs*L0;
    }

    if((mask & ghl_neutrino_rate_channel_pair) != 0 &&
       s == ghl_m1_neutrino_nux) {
      const double eta_N_channel = multiplicity * rr->eta_N_pair_cgs
                                 * number_emissivity_conversion;
      const double eta_E_channel = multiplicity * rr->eta_E_pair_mev_cgs
                                 * energy_emissivity_conversion;
      kappa_a_N += eta_N_channel / out->n_eq;
      kappa_a_E += eta_E_channel / out->J_eq;
    }
    if((mask & ghl_neutrino_rate_channel_plasmon) != 0 &&
       s == ghl_m1_neutrino_nux) {
      const double eta_N_channel = multiplicity * rr->eta_N_plasmon_cgs
                                 * number_emissivity_conversion;
      const double eta_E_channel = multiplicity * rr->eta_E_plasmon_mev_cgs
                                 * energy_emissivity_conversion;
      kappa_a_N += eta_N_channel / out->n_eq;
      kappa_a_E += eta_E_channel / out->J_eq;
    }
    if((mask & ghl_neutrino_rate_channel_bremsstrahlung) != 0 &&
       s == ghl_m1_neutrino_nux) {
      const double eta_N_channel = multiplicity * rr->eta_N_brems_cgs
                                 * number_emissivity_conversion;
      const double eta_E_channel = multiplicity * rr->eta_E_brems_mev_cgs
                                 * energy_emissivity_conversion;
      kappa_a_N += eta_N_channel / out->n_eq;
      kappa_a_E += eta_E_channel / out->J_eq;
    }

    /* Electron-flavor emission remains separate from the aggregate
     * charged-current/scattering coefficients.  The raw adapter supplies a
     * common number emissivity for the electron pair; copy both raw number
     * and energy rates without rebuilding either from a grey opacity. */
    if(s != ghl_m1_neutrino_nux) {
      if((mask & ghl_neutrino_rate_channel_pair) != 0) {
        out->eta_N_pair[ghl_m1_neutrino_pair_process_pair] =
            rr->eta_N_pair_cgs * number_emissivity_conversion;
        out->eta_E_pair[ghl_m1_neutrino_pair_process_pair] =
            rr->eta_E_pair_mev_cgs * energy_emissivity_conversion;
      }
      if((mask & ghl_neutrino_rate_channel_plasmon) != 0) {
        out->eta_N_pair[ghl_m1_neutrino_pair_process_plasmon] =
            rr->eta_N_plasmon_cgs * number_emissivity_conversion;
        out->eta_E_pair[ghl_m1_neutrino_pair_process_plasmon] =
            rr->eta_E_plasmon_mev_cgs * energy_emissivity_conversion;
      }
      if((mask & ghl_neutrino_rate_channel_bremsstrahlung) != 0) {
        out->eta_N_pair[ghl_m1_neutrino_pair_process_bremsstrahlung] =
            rr->eta_N_brems_cgs * number_emissivity_conversion;
        out->eta_E_pair[ghl_m1_neutrino_pair_process_bremsstrahlung] =
            rr->eta_E_brems_mev_cgs * energy_emissivity_conversion;
      }
    }

    out->kappa_a_N = kappa_a_N;
    out->kappa_a_E = kappa_a_E;
    if((mask & ghl_neutrino_rate_channel_nucleon_scattering) != 0)
      out->kappa_s = L0*(rr->kappa_s_E_neutron_cgs +
                         rr->kappa_s_E_proton_cgs);
    out->kappa_tr = out->kappa_a_E + out->kappa_s;
    out->eta_N = out->kappa_a_N*out->n_eq;
    out->eta_E = out->kappa_a_E*out->J_eq;
    out->eta_N_cc = out->kappa_a_N_cc*out->n_eq;
  }
  return ghl_success;
}

static void compute_staged_rates(
      const ghl_neutrino_rate_provider_context *restrict provider,
      const double rho,
      const double T,
      const double mu_e,
      const double mu_p,
      const double mu_n,
      const double X_n,
      const double X_p,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {

  const int mask = provider->channel_mask;
  const double mu_nue = mu_e + mu_p - mu_n;
  const double theta = fmax(T * provider->temperature_code_to_mev,
                            provider->min_mean_energy);
  const double thermal_number = rho / fmax(provider->baryon_mass_code, DBL_MIN);

  for(int s = 0; s < ghl_m1_neutrino_species_count; s++) {
    const ghl_m1_neutrino_species_t species = (ghl_m1_neutrino_species_t)s;
    const double lepton_weight = species_lepton_weight(species);
    const double mu_nu = (species == ghl_m1_neutrino_nue)  ? mu_nue :
                         (species == ghl_m1_neutrino_anue) ? -mu_nue : 0.0;
    const double degeneracy = mu_nu / theta;
    const double exp_degeneracy = safe_exp(degeneracy);
    const double occupancy = exp_degeneracy / (1.0 + exp_degeneracy);
    const double multiplicity = (species == ghl_m1_neutrino_nux)
                              ? provider->nu_x_multiplicity : 1.0;
    const double mean_energy =
      fmax(provider->min_mean_energy, theta * (3.1514 + 0.25 * fabs(degeneracy)));

    rates[s] = (ghl_m1_neutrino_rates){0};
    rates[s].species = species;
    rates[s].mean_energy = mean_energy;
    rates[s].lepton_weight = lepton_weight;

    double n_eq = multiplicity * thermal_number * 1.0e-3 * occupancy;
    double kappa_a_N = 0.0;
    double kappa_a_E = 0.0;
    double kappa_cc_N = 0.0;
    double kappa_s = 0.0;

    if((mask & ghl_neutrino_rate_channel_charged_current) != 0) {
      if(species == ghl_m1_neutrino_nue) {
        kappa_cc_N = provider->charged_current_scale * rho * fmax(X_n, 0.0);
        kappa_a_N += kappa_cc_N;
        kappa_a_E += kappa_cc_N * mean_energy * mean_energy;
      } else if(species == ghl_m1_neutrino_anue) {
        kappa_cc_N = provider->charged_current_scale * rho * fmax(X_p, 0.0);
        kappa_a_N += kappa_cc_N;
        kappa_a_E += kappa_cc_N * mean_energy * mean_energy;
      }
    }

    if((mask & ghl_neutrino_rate_channel_nucleon_scattering) != 0) {
      kappa_s += provider->scattering_scale * rho * fmax(X_n + X_p, 0.0)
               * mean_energy * mean_energy;
    }

    if((mask & ghl_neutrino_rate_channel_pair) != 0) {
      if(species == ghl_m1_neutrino_nux) {
        const double pair_kappa = provider->pair_scale * rho * theta * theta;
        kappa_a_N += pair_kappa;
        kappa_a_E += pair_kappa * mean_energy;
      }
      n_eq += multiplicity * thermal_number * 1.0e-4;
    }
    if((mask & ghl_neutrino_rate_channel_bremsstrahlung) != 0 &&
       species == ghl_m1_neutrino_nux) {
      const double brems_kappa = provider->bremsstrahlung_scale * rho * rho;
      kappa_a_N += brems_kappa;
      kappa_a_E += brems_kappa * mean_energy;
    }
    if((mask & ghl_neutrino_rate_channel_plasmon) != 0 &&
       species == ghl_m1_neutrino_nux) {
      const double plasmon_kappa = provider->plasmon_scale * theta * theta * theta;
      kappa_a_N += plasmon_kappa;
      kappa_a_E += plasmon_kappa * mean_energy;
    }

    rates[s].kappa_a_N = fmax(0.0, kappa_a_N);
    rates[s].kappa_a_E = fmax(0.0, kappa_a_E);
    rates[s].kappa_s = fmax(0.0, kappa_s);
    rates[s].kappa_tr = rates[s].kappa_a_E + rates[s].kappa_s;
    rates[s].n_eq = fmax(0.0, n_eq);
    rates[s].J_eq = mean_energy * rates[s].n_eq;
    rates[s].eta_N = rates[s].kappa_a_N * rates[s].n_eq;
    rates[s].kappa_a_N_cc = lepton_weight == 0.0 ? 0.0 : kappa_cc_N;
    rates[s].eta_N_cc = rates[s].kappa_a_N_cc * rates[s].n_eq;
    rates[s].eta_E = rates[s].kappa_a_E * rates[s].J_eq;
  }

  /* The table-free backend has no raw weak-rate bridge.  For every electron
   * process, use the original synthetic number opacity times the geometric
   * mean of the two electron equilibrium number targets.  This gives the
   * required common number emissivity while retaining the original
   * process-specific energy weighting kappa_E,s*J_eq,s.  Heavy flavor keeps
   * the legacy aggregate approximation and zero process arrays. */
  if((mask & (ghl_neutrino_rate_channel_pair |
              ghl_neutrino_rate_channel_plasmon |
              ghl_neutrino_rate_channel_bremsstrahlung)) != 0) {
    const double electron_n_eq_geometric_mean =
        sqrt(fmax(rates[ghl_m1_neutrino_nue].n_eq, 0.0)) *
        sqrt(fmax(rates[ghl_m1_neutrino_anue].n_eq, 0.0));
    if((mask & ghl_neutrino_rate_channel_pair) != 0) {
      const double kappa_N = provider->pair_scale * rho * theta * theta;
      const double kappa_E_nue = kappa_N *
          rates[ghl_m1_neutrino_nue].mean_energy;
      const double kappa_E_anue = kappa_N *
          rates[ghl_m1_neutrino_anue].mean_energy;
      rates[ghl_m1_neutrino_nue].eta_N_pair[
          ghl_m1_neutrino_pair_process_pair] = kappa_N * electron_n_eq_geometric_mean;
      rates[ghl_m1_neutrino_anue].eta_N_pair[
          ghl_m1_neutrino_pair_process_pair] = kappa_N * electron_n_eq_geometric_mean;
      rates[ghl_m1_neutrino_nue].eta_E_pair[
          ghl_m1_neutrino_pair_process_pair] = kappa_E_nue *
          rates[ghl_m1_neutrino_nue].J_eq;
      rates[ghl_m1_neutrino_anue].eta_E_pair[
          ghl_m1_neutrino_pair_process_pair] = kappa_E_anue *
          rates[ghl_m1_neutrino_anue].J_eq;
    }
    if((mask & ghl_neutrino_rate_channel_plasmon) != 0) {
      const double kappa_N = provider->plasmon_scale * theta * theta * theta;
      const double kappa_E_nue = kappa_N *
          rates[ghl_m1_neutrino_nue].mean_energy;
      const double kappa_E_anue = kappa_N *
          rates[ghl_m1_neutrino_anue].mean_energy;
      rates[ghl_m1_neutrino_nue].eta_N_pair[
          ghl_m1_neutrino_pair_process_plasmon] = kappa_N * electron_n_eq_geometric_mean;
      rates[ghl_m1_neutrino_anue].eta_N_pair[
          ghl_m1_neutrino_pair_process_plasmon] = kappa_N * electron_n_eq_geometric_mean;
      rates[ghl_m1_neutrino_nue].eta_E_pair[
          ghl_m1_neutrino_pair_process_plasmon] = kappa_E_nue *
          rates[ghl_m1_neutrino_nue].J_eq;
      rates[ghl_m1_neutrino_anue].eta_E_pair[
          ghl_m1_neutrino_pair_process_plasmon] = kappa_E_anue *
          rates[ghl_m1_neutrino_anue].J_eq;
    }
    if((mask & ghl_neutrino_rate_channel_bremsstrahlung) != 0) {
      const double kappa_N = provider->bremsstrahlung_scale * rho * rho;
      const double kappa_E_nue = kappa_N *
          rates[ghl_m1_neutrino_nue].mean_energy;
      const double kappa_E_anue = kappa_N *
          rates[ghl_m1_neutrino_anue].mean_energy;
      rates[ghl_m1_neutrino_nue].eta_N_pair[
          ghl_m1_neutrino_pair_process_bremsstrahlung] = kappa_N * electron_n_eq_geometric_mean;
      rates[ghl_m1_neutrino_anue].eta_N_pair[
          ghl_m1_neutrino_pair_process_bremsstrahlung] = kappa_N * electron_n_eq_geometric_mean;
      rates[ghl_m1_neutrino_nue].eta_E_pair[
          ghl_m1_neutrino_pair_process_bremsstrahlung] = kappa_E_nue *
          rates[ghl_m1_neutrino_nue].J_eq;
      rates[ghl_m1_neutrino_anue].eta_E_pair[
          ghl_m1_neutrino_pair_process_bremsstrahlung] = kappa_E_anue *
          rates[ghl_m1_neutrino_anue].J_eq;
    }
  }
}

static ghl_error_codes_t recover_failure(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_cache *restrict cache,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict eos,
      const double rho, const double T, const double Ye,
      const bool recovered_key_valid,
      const ghl_error_codes_t error,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {

  if(diagnostics != NULL) {
    diagnostics->failures++;
    diagnostics->last_error = error;
  }

  ghl_m1_neutrino_rates candidate[ghl_m1_neutrino_species_count];

  switch(provider->failure_policy) {
    case ghl_neutrino_rate_failure_transparent:
      fill_transparent_rates(candidate, provider->min_mean_energy);
      if(publish_recovered_rates(candidate, rates) != ghl_success)
        return error;
      if(diagnostics != NULL) {
        diagnostics->transparent_recoveries++;
        diagnostics->last_recovery = ghl_neutrino_rate_recovery_transparent;
      }
      return ghl_success;
    case ghl_neutrino_rate_failure_equilibrium:
      fill_transparent_rates(candidate, provider->min_mean_energy);
      for(int s = 0; s < ghl_m1_neutrino_species_count; s++) {
        candidate[s].n_eq = 1.0e-30;
        candidate[s].J_eq = candidate[s].mean_energy * candidate[s].n_eq;
        candidate[s].kappa_a_N = provider->equilibrium_recovery_rate;
        candidate[s].kappa_a_E = provider->equilibrium_recovery_rate;
        candidate[s].kappa_tr = candidate[s].kappa_a_E;
        candidate[s].eta_N = candidate[s].kappa_a_N * candidate[s].n_eq;
        candidate[s].eta_E = candidate[s].kappa_a_E * candidate[s].J_eq;
        if(candidate[s].lepton_weight != 0.0) {
          candidate[s].kappa_a_N_cc = candidate[s].kappa_a_N;
          candidate[s].eta_N_cc = candidate[s].eta_N;
        }
      }
      if(publish_recovered_rates(candidate, rates) != ghl_success)
        return error;
      if(diagnostics != NULL) {
        diagnostics->equilibrium_recoveries++;
        diagnostics->last_recovery = ghl_neutrino_rate_recovery_equilibrium;
      }
      return ghl_success;
    case ghl_neutrino_rate_failure_hold_last:
      if(!recovered_key_valid ||
         !same_rates_key(cache, provider, eos, rho, T, Ye))
        return error;
      memcpy(rates, cache->rates, sizeof(cache->rates));
      if(diagnostics != NULL) {
        diagnostics->hold_last_recoveries++;
        diagnostics->last_recovery = ghl_neutrino_rate_recovery_hold_last;
      }
      return ghl_success;
    case ghl_neutrino_rate_failure_abort:
    default:
      return error;
  }
}

ghl_error_codes_t ghl_neutrino_rate_provider_compute_cell(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_cache *restrict cache,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {

  if(provider == NULL || prims == NULL || rates == NULL)
    return ghl_error_m1_null_pointer;

  if(diagnostics != NULL)
    diagnostics->last_recovery = ghl_neutrino_rate_recovery_none;

  /* The built-in provider's nu_x bundle is already the summed heavy-lepton
   * sector. A different multiplicity is a provider-boundary configuration
   * error, not a recoverable cell-microphysics failure: reject it before a
   * cache lookup or any rate calculation so a host cannot double count nu_x.
   * External providers own their flavor conversion and return an already-summed
   * bundle through the same frozen-rate interface. */
  const ghl_error_codes_t context_error = validate_provider_context(provider, eos);
  if(context_error != ghl_success) {
    if(diagnostics != NULL) {
      diagnostics->failures++;
      diagnostics->last_error = context_error;
    }
    return context_error;
  }

  if(diagnostics != NULL) {
    diagnostics->active_channel_mask = provider->channel_mask;
    diagnostics->last_error = ghl_success;
  }

  ghl_neutrino_rate_provider_cache staged_cache;
  ghl_neutrino_rate_provider_cache *work_cache = NULL;
  if(cache != NULL) {
    staged_cache = *cache;
    work_cache = &staged_cache;
  }

  double rho = prims->rho;
  double T = prims->temperature;
  double Ye = prims->Y_e;

  double muhat = 0.0, mu_e = 0.0, mu_p = 0.0, mu_n = 0.0;
  double X_n = 0.0, X_p = 0.0;
  bool recovered_key_valid = false;
  ghl_error_codes_t err = compute_thermo(
      provider, work_cache, diagnostics, eos, prims,
      &rho, &T, &Ye, &muhat, &mu_e, &mu_p, &mu_n, &X_n, &X_p,
      &recovered_key_valid);
  if(err != ghl_success)
    return recover_failure(provider, cache, diagnostics, eos,
                           rho, T, Ye, recovered_key_valid, err, rates);

  if(same_rates_key(cache, provider, eos, rho, T, Ye)) {
    memcpy(rates, cache->rates, sizeof(cache->rates));
    if(diagnostics != NULL) {
      diagnostics->cache_hits++;
      memcpy(diagnostics->beta_kirchhoff_relative_mismatch,
             cache->beta_kirchhoff_relative_mismatch,
             sizeof(cache->beta_kirchhoff_relative_mismatch));
      memcpy(diagnostics->beta_kirchhoff_mismatch_valid,
             cache->beta_kirchhoff_mismatch_valid,
             sizeof(cache->beta_kirchhoff_mismatch_valid));
    }
    return ghl_success;
  }

  if(diagnostics != NULL)
    diagnostics->cache_misses++;

  ghl_m1_neutrino_rates candidate_rates[ghl_m1_neutrino_species_count];
  double candidate_mismatch[ghl_m1_neutrino_species_count] = {0};
  bool candidate_mismatch_valid[ghl_m1_neutrino_species_count] = {false};
  if(provider->backend == ghl_neutrino_rate_backend_nrpyleakage) {
    err = assemble_nrpyleakage_rates(
        provider, rho, T, Ye, muhat, mu_e, mu_p, mu_n, X_n, X_p,
        candidate_rates, candidate_mismatch, candidate_mismatch_valid);
    if(err != ghl_success)
      return recover_failure(provider, cache, diagnostics, eos,
                             rho, T, Ye, true, err, rates);
  } else {
    compute_staged_rates(
        provider, rho, T, mu_e, mu_p, mu_n, X_n, X_p, candidate_rates);
  }

  for(int s = 0; s < ghl_m1_neutrino_species_count; s++) {
    err = ghl_m1_validate_neutrino_rates(&candidate_rates[s], NULL);
    if(err != ghl_success)
      return recover_failure(provider, cache, diagnostics, eos,
                             rho, T, Ye, true, err, rates);
  }

  if(cache != NULL) {
    staged_cache.rates_valid = true;
    staged_cache.rho = rho;
    staged_cache.T = T;
    staged_cache.Ye = Ye;
    memcpy(staged_cache.rates, candidate_rates, sizeof(staged_cache.rates));
    memcpy(staged_cache.beta_kirchhoff_relative_mismatch,
           candidate_mismatch, sizeof(candidate_mismatch));
    memcpy(staged_cache.beta_kirchhoff_mismatch_valid,
           candidate_mismatch_valid, sizeof(candidate_mismatch_valid));
    staged_cache.provider_snapshot = *provider;
    staged_cache.eos_snapshot = eos;
    *cache = staged_cache;
  }
  memcpy(rates, candidate_rates, sizeof(candidate_rates));
  if(diagnostics != NULL) {
    memcpy(diagnostics->beta_kirchhoff_relative_mismatch,
           candidate_mismatch, sizeof(candidate_mismatch));
    memcpy(diagnostics->beta_kirchhoff_mismatch_valid,
           candidate_mismatch_valid, sizeof(candidate_mismatch_valid));
  }
  return ghl_success;
}
