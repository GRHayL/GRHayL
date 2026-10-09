#include "ghl_neutrino_rate_provider.h"
#include "ghl_m1_neutrino_rate_backend.h"
#include <math.h>
#include <string.h>

void ghl_neutrino_rate_provider_cache_initialize(
      ghl_neutrino_rate_provider_cache *restrict cache) {
  if(cache != NULL) {
    *cache = (ghl_neutrino_rate_provider_cache){ 0 };
  }
}

ghl_error_codes_t ghl_neutrino_rate_provider_initialize_nrpyleakage(
      ghl_neutrino_rate_provider_context *restrict provider) {
  if(provider == NULL) {
    return ghl_error_m1_null_pointer;
  }
#ifdef GHL_DISABLE_HDF5
  return ghl_m1_neutrino_rate_backend_initialize();
#else
  /* The enabled backend initializer is a pure, unconditional success. */
  *provider = (ghl_neutrino_rate_provider_context){ 0 };
  provider->channel_mask
        = ghl_neutrino_rate_channel_charged_current
          | ghl_neutrino_rate_channel_nucleon_scattering | ghl_neutrino_rate_channel_pair
          | ghl_neutrino_rate_channel_bremsstrahlung | ghl_neutrino_rate_channel_plasmon;
  provider->failure_policy = ghl_neutrino_rate_failure_return_error;
  provider->table_bounds_policy = ghl_neutrino_rate_table_bounds_return_error;
  provider->nu_x_multiplicity = 4.0;
  provider->eos_generation = 0;
  provider->equilibrium_recovery_rate = 1.0;
  return ghl_success;
#endif
}

ghl_error_codes_t ghl_neutrino_rate_provider_initialize_default(
      ghl_neutrino_rate_provider_context *restrict provider) {
  return ghl_neutrino_rate_provider_initialize_nrpyleakage(provider);
}

static ghl_error_codes_t validate_provider_context(
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_eos_parameters *restrict eos) {
  const int valid_channels
        = ghl_neutrino_rate_channel_charged_current
          | ghl_neutrino_rate_channel_nucleon_scattering | ghl_neutrino_rate_channel_pair
          | ghl_neutrino_rate_channel_bremsstrahlung | ghl_neutrino_rate_channel_plasmon;
  if(!isfinite(provider->nu_x_multiplicity) || provider->nu_x_multiplicity != 4.0
     || (provider->channel_mask & ~valid_channels) != 0
     || provider->failure_policy < ghl_neutrino_rate_failure_return_error
     || provider->failure_policy > ghl_neutrino_rate_failure_equilibrium
     || provider->table_bounds_policy < ghl_neutrino_rate_table_bounds_return_error
     || provider->table_bounds_policy > ghl_neutrino_rate_table_bounds_clamp
     || !isfinite(provider->equilibrium_recovery_rate)
     || provider->equilibrium_recovery_rate <= 0.0) {
    return ghl_error_m1_microphysics_failure;
  }
#ifdef GHL_DISABLE_HDF5
  /* Context validation precedes the disabled-backend error; EOS and primitive
   * validation, cache lookup, and recovery never run in this configuration. */
  (void)eos;
  return ghl_m1_neutrino_rate_backend_initialize();
#else
  if(eos == NULL || eos->eos_type != ghl_eos_tabulated
     || eos->table_type == ghl_eos_table_unknown) {
    return ghl_error_m1_microphysics_failure;
  }
  return ghl_success;
#endif
}

#ifndef GHL_DISABLE_HDF5
/* Table-backed cache, thermodynamics, and recovery implementation. */
static bool same_provider_configuration(
      const ghl_neutrino_rate_provider_context *restrict lhs,
      const ghl_neutrino_rate_provider_context *restrict rhs) {
  /* Require the complete provider snapshot to match before cache reuse. */
  return lhs->channel_mask == rhs->channel_mask
         && lhs->failure_policy == rhs->failure_policy
         && lhs->table_bounds_policy == rhs->table_bounds_policy
         && lhs->nu_x_multiplicity == rhs->nu_x_multiplicity
         && lhs->eos_generation == rhs->eos_generation
         && lhs->equilibrium_recovery_rate == rhs->equilibrium_recovery_rate;
}

static bool same_provenance(
      const ghl_neutrino_rate_provider_cache *restrict cache,
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_eos_parameters *restrict eos) {
  return cache != NULL && cache->eos_snapshot == eos
         && same_provider_configuration(&cache->provider_snapshot, provider);
}

static bool same_thermo_key(
      const ghl_neutrino_rate_provider_cache *restrict cache,
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double T,
      const double Ye) {
  /* A thermodynamic cache hit requires provenance, validity, and all keys. */
  return same_provenance(cache, provider, eos) && cache->thermo_valid
         && cache->thermo_rho == rho && cache->thermo_T == T && cache->thermo_Ye == Ye;
}

static bool same_rates_key(
      const ghl_neutrino_rate_provider_cache *restrict cache,
      const ghl_neutrino_rate_provider_context *restrict provider,
      const ghl_eos_parameters *restrict eos,
      const double rho,
      const double T,
      const double Ye) {
  /* Final rates have their own validity flag and primitive keys. */
  return same_provenance(cache, provider, eos) && cache->rates_valid && cache->rho == rho
         && cache->T == T && cache->Ye == Ye;
}

static ghl_error_codes_t validate_inputs(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict eos,
      double *restrict rho,
      double *restrict T,
      double *restrict Ye) {
  if(!isfinite(*rho) || !isfinite(*Ye) || *rho <= 0.0 || *Ye < 0.0 || *Ye > 1.0) {
    return ghl_error_m1_microphysics_failure;
  }
  bool out_of_bounds = false;
  const double rho_lo = eos->table_rho_min > 0.0 ? eos->table_rho_min : eos->rho_min;
  const double rho_hi = eos->table_rho_max > 0.0 ? eos->table_rho_max : eos->rho_max;
  const double T_lo = eos->table_T_min > 0.0 ? eos->table_T_min : eos->T_min;
  const double T_hi = eos->table_T_max > 0.0 ? eos->table_T_max : eos->T_max;
  const double Ye_lo = (eos->table_Y_e_min > 0.0
                        || (eos->table_Y_e_min == 0.0 && eos->table_Y_e_max > 0.0))
                             ? eos->table_Y_e_min
                             : eos->Y_e_min;
  const double Ye_hi = eos->table_Y_e_max > 0.0 ? eos->table_Y_e_max : eos->Y_e_max;
  if(!isfinite(rho_lo) || !isfinite(rho_hi) || rho_lo <= 0.0 || rho_lo >= rho_hi
     || !isfinite(T_lo) || !isfinite(T_hi) || T_lo <= 0.0 || T_lo >= T_hi
     || !isfinite(Ye_lo) || !isfinite(Ye_hi) || Ye_lo < 0.0 || Ye_hi > 1.0
     || Ye_lo >= Ye_hi) {
    return ghl_error_m1_microphysics_failure;
  }
  out_of_bounds = *rho < rho_lo || *rho > rho_hi || *T < T_lo || *T > T_hi || *Ye < Ye_lo
                  || *Ye > Ye_hi;
  if(!out_of_bounds) {
    return ghl_success;
  }
  if(diagnostics != NULL) {
    diagnostics->table_bound_hits++;
  }
  if(provider->table_bounds_policy != ghl_neutrino_rate_table_bounds_clamp) {
    return ghl_error_m1_microphysics_failure;
  }
  *rho = ghl_clamp(*rho, rho_lo, rho_hi);
  *T = ghl_clamp(*T, T_lo, T_hi);
  *Ye = ghl_clamp(*Ye, Ye_lo, Ye_hi);
  if(diagnostics != NULL) {
    diagnostics->clamped_inputs++;
  }
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
      bool *restrict physical_thermo_valid) {
  *physical_thermo_valid = false;
  *rho = prims->rho;
  *Ye = prims->Y_e;
  *T = prims->temperature;
  if(!isfinite(*T) || *T <= 0.0) {
    const double T_lo = eos->table_T_min > 0.0 ? eos->table_T_min : eos->T_min;
    const double T_hi = eos->table_T_max > 0.0 ? eos->table_T_max : eos->T_max;
    if(!isfinite(T_lo) || !isfinite(T_hi) || T_lo <= 0.0 || T_lo >= T_hi
       || !isfinite(prims->eps)) {
      return ghl_error_m1_microphysics_failure;
    }
    *T = exp(0.5 * (log(T_lo) + log(T_hi)));
    ghl_error_codes_t error = validate_inputs(provider, diagnostics, eos, rho, T, Ye);
    if(error != ghl_success) {
      return error;
    }
    error = ghl_m1_neutrino_rate_backend_temperature_from_eps(
          eos, *rho, *Ye, prims->eps, T);
    if(error != ghl_success) {
      return error;
    }
  }
  ghl_error_codes_t error = validate_inputs(provider, diagnostics, eos, rho, T, Ye);
  if(error != ghl_success) {
    return error;
  }

  if(same_thermo_key(cache, provider, eos, *rho, *T, *Ye)) {
    *muhat = cache->muhat;
    *mu_e = cache->mu_e;
    *mu_p = cache->mu_p;
    *mu_n = cache->mu_n;
    *X_n = cache->X_n;
    *X_p = cache->X_p;
    *physical_thermo_valid = true;
    return ghl_success;
  }

  error = ghl_m1_neutrino_rate_backend_thermo_from_T(
        eos, *rho, *Ye, *T, muhat, mu_e, mu_p, mu_n, X_n, X_p);
  if(error != ghl_success) {
    return error;
  }
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
  *physical_thermo_valid = true;
  return ghl_success;
}

static ghl_error_codes_t publish_recovered_rates(
      const ghl_m1_neutrino_rates candidate[ghl_m1_neutrino_species_count],
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    const ghl_error_codes_t error
          = ghl_m1_validate_neutrino_rates(&candidate[species], NULL);
    /* Equilibrium recovery scales the assembled targets by the configured
     * recovery rate; finite factors can still overflow during that scaling. */
    if(error != ghl_success) {
      return error;
    }
  }
  memcpy(
        rates, candidate, sizeof(ghl_m1_neutrino_rates) * ghl_m1_neutrino_species_count);
  return ghl_success;
}

static ghl_error_codes_t recover_failure(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const bool physical_thermo_valid,
      const double rho,
      const double T,
      const double Ye,
      const double muhat,
      const double mu_e,
      const double mu_p,
      const double mu_n,
      const double X_n,
      const double X_p,
      const ghl_error_codes_t original_error,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {
  if(diagnostics != NULL) {
    diagnostics->failures++;
    diagnostics->last_error = original_error;
  }
  if(!physical_thermo_valid
     || provider->failure_policy == ghl_neutrino_rate_failure_return_error) {
    return original_error;
  }

  ghl_m1_neutrino_rates candidate[ghl_m1_neutrino_species_count];
  double mismatch[ghl_m1_neutrino_species_count] = { 0.0 };
  bool mismatch_valid[ghl_m1_neutrino_species_count] = { false };
  const ghl_error_codes_t target_error = ghl_m1_neutrino_rate_backend_assemble(
        0, provider->nu_x_multiplicity, rho, T, Ye, muhat, mu_e, mu_p, mu_n, X_n, X_p,
        candidate, mismatch, mismatch_valid);
  /* A failed primary assembly does not guarantee representable mask-zero
   * targets; underflow in the same thermodynamics can reject recovery too. */
  if(target_error != ghl_success) {
    return original_error;
  }
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    if(ghl_m1_validate_neutrino_rates(&candidate[species], NULL) != ghl_success) {
      return original_error;
    }
  }

  if(provider->failure_policy == ghl_neutrino_rate_failure_equilibrium) {
    for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
      ghl_m1_neutrino_rates *const rate = &candidate[species];
      rate->kappa_a_N = provider->equilibrium_recovery_rate;
      rate->kappa_a_E = provider->equilibrium_recovery_rate;
      rate->kappa_tr = rate->kappa_a_E;
      rate->eta_N = rate->kappa_a_N * rate->n_eq;
      rate->eta_E = rate->kappa_a_E * rate->J_eq;
      if(rate->lepton_weight != 0.0) {
        rate->kappa_a_N_cc = provider->equilibrium_recovery_rate;
        rate->eta_N_cc = rate->kappa_a_N_cc * rate->n_eq;
      }
    }
  }

  /* The configured equilibrium scaling above can produce invalid rates;
   * validate the final recovery candidate before publishing any species. */
  if(publish_recovered_rates(candidate, rates) != ghl_success) {
    return original_error;
  }
  if(diagnostics != NULL) {
    memcpy(diagnostics->beta_kirchhoff_relative_mismatch, mismatch, sizeof(mismatch));
    memcpy(
          diagnostics->beta_kirchhoff_mismatch_valid, mismatch_valid,
          sizeof(mismatch_valid));
    if(provider->failure_policy == ghl_neutrino_rate_failure_transparent) {
      diagnostics->transparent_recoveries++;
      diagnostics->last_recovery = ghl_neutrino_rate_recovery_transparent;
    }
    else {
      diagnostics->equilibrium_recoveries++;
      diagnostics->last_recovery = ghl_neutrino_rate_recovery_equilibrium;
    }
  }
  return original_error;
}

#endif

#ifdef GHL_DISABLE_HDF5
/* Preserve the public failure boundary without compiling unavailable table
 * operations. Required-pointer errors leave diagnostics unchanged; valid
 * pointers reset recovery before context/backend validation. */
ghl_error_codes_t ghl_neutrino_rate_provider_compute_cell(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_cache *restrict cache,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {
  (void)cache;
  if(provider == NULL || prims == NULL || rates == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(diagnostics != NULL) {
    diagnostics->last_recovery = ghl_neutrino_rate_recovery_none;
  }
  const ghl_error_codes_t error = validate_provider_context(provider, eos);
  if(diagnostics != NULL) {
    diagnostics->failures++;
    diagnostics->last_error = error;
  }
  return error;
}
#else
ghl_error_codes_t ghl_neutrino_rate_provider_compute_cell(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_cache *restrict cache,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {
  if(provider == NULL || prims == NULL || rates == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(diagnostics != NULL) {
    diagnostics->last_recovery = ghl_neutrino_rate_recovery_none;
  }
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
  double rho = 0.0, T = 0.0, Ye = 0.0;
  double muhat = 0.0, mu_e = 0.0, mu_p = 0.0, mu_n = 0.0;
  double X_n = 0.0, X_p = 0.0;
  bool physical_thermo_valid = false;
  ghl_error_codes_t error = compute_thermo(
        provider, work_cache, diagnostics, eos, prims, &rho, &T, &Ye, &muhat, &mu_e,
        &mu_p, &mu_n, &X_n, &X_p, &physical_thermo_valid);
  if(error != ghl_success) {
    return recover_failure(
          provider, diagnostics, false, rho, T, Ye, muhat, mu_e, mu_p, mu_n, X_n, X_p,
          error, rates);
  }

  if(same_rates_key(cache, provider, eos, rho, T, Ye)) {
    memcpy(rates, cache->rates, sizeof(cache->rates));
    if(diagnostics != NULL) {
      diagnostics->cache_hits++;
      memcpy(
            diagnostics->beta_kirchhoff_relative_mismatch,
            cache->beta_kirchhoff_relative_mismatch,
            sizeof(cache->beta_kirchhoff_relative_mismatch));
      memcpy(
            diagnostics->beta_kirchhoff_mismatch_valid,
            cache->beta_kirchhoff_mismatch_valid,
            sizeof(cache->beta_kirchhoff_mismatch_valid));
    }
    return ghl_success;
  }
  if(diagnostics != NULL) {
    diagnostics->cache_misses++;
  }

  ghl_m1_neutrino_rates candidate_rates[ghl_m1_neutrino_species_count];
  double candidate_mismatch[ghl_m1_neutrino_species_count] = { 0.0 };
  bool candidate_mismatch_valid[ghl_m1_neutrino_species_count] = { false };
  error = ghl_m1_neutrino_rate_backend_assemble(
        provider->channel_mask, provider->nu_x_multiplicity, rho, T, Ye, muhat, mu_e,
        mu_p, mu_n, X_n, X_p, candidate_rates, candidate_mismatch,
        candidate_mismatch_valid);
  if(error != ghl_success) {
    return recover_failure(
          provider, diagnostics, physical_thermo_valid, rho, T, Ye, muhat, mu_e, mu_p,
          mu_n, X_n, X_p, error, rates);
  }
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    error = ghl_m1_validate_neutrino_rates(&candidate_rates[species], NULL);
    /* Validate converted output and the required rounding mode even after
     * the raw backend accepted its thermodynamic inputs. */
    if(error != ghl_success) {
      return recover_failure(
            provider, diagnostics, physical_thermo_valid, rho, T, Ye, muhat, mu_e, mu_p,
            mu_n, X_n, X_p, error, rates);
    }
  }

  if(cache != NULL) {
    staged_cache.rates_valid = true;
    staged_cache.rho = rho;
    staged_cache.T = T;
    staged_cache.Ye = Ye;
    memcpy(staged_cache.rates, candidate_rates, sizeof(staged_cache.rates));
    memcpy(
          staged_cache.beta_kirchhoff_relative_mismatch, candidate_mismatch,
          sizeof(candidate_mismatch));
    memcpy(
          staged_cache.beta_kirchhoff_mismatch_valid, candidate_mismatch_valid,
          sizeof(candidate_mismatch_valid));
    staged_cache.provider_snapshot = *provider;
    staged_cache.eos_snapshot = eos;
    *cache = staged_cache;
  }
  memcpy(rates, candidate_rates, sizeof(candidate_rates));
  if(diagnostics != NULL) {
    memcpy(
          diagnostics->beta_kirchhoff_relative_mismatch, candidate_mismatch,
          sizeof(candidate_mismatch));
    memcpy(
          diagnostics->beta_kirchhoff_mismatch_valid, candidate_mismatch_valid,
          sizeof(candidate_mismatch_valid));
  }
  return ghl_success;
}

#endif
