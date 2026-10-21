/* Test-only synthetic reference-rate model. This file intentionally contains
 * no production EOS, NRPyLeakage, or provider implementation. */

#include <float.h>
#include <math.h>
#include <string.h>

static double m1_test_reference_safe_exp(const double x) {
  return exp(ghl_clamp(x, -40.0, 40.0));
}

static double m1_test_reference_lepton_weight(const ghl_m1_neutrino_species_t species) {
  static const double weights[ghl_m1_neutrino_species_count]
        = { [ghl_m1_neutrino_nue] = 1.0,
            [ghl_m1_neutrino_anue] = -1.0,
            [ghl_m1_neutrino_nux] = 0.0 };
  return weights[species];
}

static ghl_error_codes_t m1_test_reference_validate_context(
      const m1_test_reference_provider_context *restrict provider) {
  const int valid_channels
        = ghl_neutrino_rate_channel_charged_current
          | ghl_neutrino_rate_channel_nucleon_scattering | ghl_neutrino_rate_channel_pair
          | ghl_neutrino_rate_channel_bremsstrahlung | ghl_neutrino_rate_channel_plasmon;
  const double nonnegative[]
        = { provider->charged_current_scale, provider->scattering_scale,
            provider->pair_scale, provider->bremsstrahlung_scale,
            provider->plasmon_scale };
  if(!isfinite(provider->nu_x_multiplicity) || provider->nu_x_multiplicity != 4.0
     || (provider->channel_mask & ~valid_channels) != 0
     || provider->failure_policy < ghl_neutrino_rate_failure_return_error
     || provider->failure_policy > ghl_neutrino_rate_failure_equilibrium
     || provider->table_bounds_policy < ghl_neutrino_rate_table_bounds_return_error
     || provider->table_bounds_policy > ghl_neutrino_rate_table_bounds_clamp
     || !isfinite(provider->temperature_code_to_mev)
     || provider->temperature_code_to_mev <= 0.0 || !isfinite(provider->baryon_mass_code)
     || provider->baryon_mass_code <= 0.0 || !isfinite(provider->min_mean_energy)
     || provider->min_mean_energy <= 0.0
     || !isfinite(provider->equilibrium_recovery_rate)
     || provider->equilibrium_recovery_rate <= 0.0) {
    return ghl_error_m1_microphysics_failure;
  }
  for(size_t i = 0; i < sizeof(nonnegative) / sizeof(nonnegative[0]); ++i) {
    if(!isfinite(nonnegative[i]) || nonnegative[i] < 0.0) {
      return ghl_error_m1_microphysics_failure;
    }
  }
  return ghl_success;
}

static bool m1_test_reference_same_configuration(
      const m1_test_reference_provider_context *restrict lhs,
      const m1_test_reference_provider_context *restrict rhs) {
  return lhs->channel_mask == rhs->channel_mask
         && lhs->failure_policy == rhs->failure_policy
         && lhs->table_bounds_policy == rhs->table_bounds_policy
         && lhs->nu_x_multiplicity == rhs->nu_x_multiplicity
         && lhs->eos_generation == rhs->eos_generation
         && lhs->temperature_code_to_mev == rhs->temperature_code_to_mev
         && lhs->baryon_mass_code == rhs->baryon_mass_code
         && lhs->charged_current_scale == rhs->charged_current_scale
         && lhs->scattering_scale == rhs->scattering_scale
         && lhs->pair_scale == rhs->pair_scale
         && lhs->bremsstrahlung_scale == rhs->bremsstrahlung_scale
         && lhs->plasmon_scale == rhs->plasmon_scale
         && lhs->min_mean_energy == rhs->min_mean_energy
         && lhs->equilibrium_recovery_rate == rhs->equilibrium_recovery_rate;
}

static bool m1_test_reference_same_provenance(
      const m1_test_reference_provider_cache *restrict cache,
      const m1_test_reference_provider_context *restrict provider,
      const ghl_eos_parameters *restrict eos) {
  return cache != NULL && cache->eos_snapshot == eos
         && m1_test_reference_same_configuration(&cache->provider_snapshot, provider);
}

static ghl_error_codes_t m1_test_reference_compute_thermo(
      const m1_test_reference_provider_context *restrict provider,
      m1_test_reference_provider_cache *restrict cache,
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
      double *restrict X_p) {
  *rho = prims->rho;
  *T = prims->temperature;
  *Ye = prims->Y_e;
  if(!isfinite(*rho) || *rho <= 0.0 || !isfinite(*Ye) || *Ye < 0.0 || *Ye > 1.0) {
    return ghl_error_m1_microphysics_failure;
  }
  if(!isfinite(*T) || *T <= 0.0) {
    if(!isfinite(prims->eps) || prims->eps <= 0.0) {
      return ghl_error_m1_microphysics_failure;
    }
    *T = prims->eps;
  }
  if(!isfinite(*T) || *T <= 0.0) {
    return ghl_error_m1_microphysics_failure;
  }
  if(cache != NULL && m1_test_reference_same_provenance(cache, provider, eos)
     && cache->thermo_valid && cache->thermo_rho == *rho && cache->thermo_T == *T
     && cache->thermo_Ye == *Ye) {
    *muhat = cache->muhat;
    *mu_e = cache->mu_e;
    *mu_p = cache->mu_p;
    *mu_n = cache->mu_n;
    *X_n = cache->X_n;
    *X_p = cache->X_p;
    return ghl_success;
  }

  const double ye = ghl_clamp(*Ye, 1.0e-12, 1.0 - 1.0e-12);
  *X_p = ye;
  *X_n = 1.0 - ye;
  *mu_e = *T * log(ye / (1.0 - ye));
  *mu_p = *T * log(fmax(*X_p, 1.0e-12));
  *mu_n = *T * log(fmax(*X_n, 1.0e-12));
  *muhat = *mu_n - *mu_p;
  if(!isfinite(*muhat) || !isfinite(*mu_e) || !isfinite(*mu_p) || !isfinite(*mu_n)) {
    return ghl_error_m1_microphysics_failure;
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
  return ghl_success;
}

static void m1_test_reference_compute_rates(
      const m1_test_reference_provider_context *restrict provider,
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
  const double theta
        = fmax(T * provider->temperature_code_to_mev, provider->min_mean_energy);
  const double thermal_number = rho / fmax(provider->baryon_mass_code, DBL_MIN);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    const ghl_m1_neutrino_species_t index = (ghl_m1_neutrino_species_t)species;
    const double mu_nu = index == ghl_m1_neutrino_nue    ? mu_nue
                         : index == ghl_m1_neutrino_anue ? -mu_nue
                                                         : 0.0;
    const double degeneracy = mu_nu / theta;
    const double exp_degeneracy = m1_test_reference_safe_exp(degeneracy);
    const double occupancy = exp_degeneracy / (1.0 + exp_degeneracy);
    const double multiplicity
          = index == ghl_m1_neutrino_nux ? provider->nu_x_multiplicity : 1.0;
    const double mean_energy
          = fmax(provider->min_mean_energy, theta * (3.1514 + 0.25 * fabs(degeneracy)));
    const double lepton_weight = m1_test_reference_lepton_weight(index);
    double n_eq = multiplicity * thermal_number * 1.0e-3 * occupancy;
    double kappa_a_N = 0.0, kappa_a_E = 0.0, kappa_cc_N = 0.0, kappa_s = 0.0;

    rates[species] = (ghl_m1_neutrino_rates){ 0 };
    rates[species].species = index;
    rates[species].mean_energy = mean_energy;
    rates[species].lepton_weight = lepton_weight;
    if((mask & ghl_neutrino_rate_channel_charged_current) != 0) {
      if(index == ghl_m1_neutrino_nue) {
        kappa_cc_N = provider->charged_current_scale * rho * fmax(X_n, 0.0);
      }
      else if(index == ghl_m1_neutrino_anue) {
        kappa_cc_N = provider->charged_current_scale * rho * fmax(X_p, 0.0);
      }
      kappa_a_N += kappa_cc_N;
      kappa_a_E += kappa_cc_N * mean_energy * mean_energy;
    }
    if((mask & ghl_neutrino_rate_channel_nucleon_scattering) != 0) {
      kappa_s = provider->scattering_scale * rho * fmax(X_n + X_p, 0.0) * mean_energy
                * mean_energy;
    }
    if((mask & ghl_neutrino_rate_channel_pair) != 0) {
      if(index == ghl_m1_neutrino_nux) {
        const double pair_kappa = provider->pair_scale * rho * theta * theta;
        kappa_a_N += pair_kappa;
        kappa_a_E += pair_kappa * mean_energy;
      }
      n_eq += multiplicity * thermal_number * 1.0e-4;
    }
    if((mask & ghl_neutrino_rate_channel_bremsstrahlung) != 0
       && index == ghl_m1_neutrino_nux) {
      const double brems_kappa = provider->bremsstrahlung_scale * rho * rho;
      kappa_a_N += brems_kappa;
      kappa_a_E += brems_kappa * mean_energy;
    }
    if((mask & ghl_neutrino_rate_channel_plasmon) != 0 && index == ghl_m1_neutrino_nux) {
      const double plasmon_kappa = provider->plasmon_scale * theta * theta * theta;
      kappa_a_N += plasmon_kappa;
      kappa_a_E += plasmon_kappa * mean_energy;
    }
    rates[species].kappa_a_N = fmax(0.0, kappa_a_N);
    rates[species].kappa_a_E = fmax(0.0, kappa_a_E);
    rates[species].kappa_s = fmax(0.0, kappa_s);
    rates[species].kappa_tr = rates[species].kappa_a_E + rates[species].kappa_s;
    rates[species].n_eq = fmax(0.0, n_eq);
    rates[species].J_eq = mean_energy * rates[species].n_eq;
    rates[species].eta_N = rates[species].kappa_a_N * rates[species].n_eq;
    rates[species].kappa_a_N_cc = lepton_weight == 0.0 ? 0.0 : kappa_cc_N;
    rates[species].eta_N_cc = rates[species].kappa_a_N_cc * rates[species].n_eq;
    rates[species].eta_E = rates[species].kappa_a_E * rates[species].J_eq;
  }

  if((mask
      & (ghl_neutrino_rate_channel_pair | ghl_neutrino_rate_channel_plasmon
         | ghl_neutrino_rate_channel_bremsstrahlung))
     != 0) {
    const double electron_n_eq_geometric_mean
          = sqrt(fmax(rates[ghl_m1_neutrino_nue].n_eq, 0.0))
            * sqrt(fmax(rates[ghl_m1_neutrino_anue].n_eq, 0.0));
    const int process_masks[ghl_m1_neutrino_pair_process_count]
          = { ghl_neutrino_rate_channel_pair, ghl_neutrino_rate_channel_plasmon,
              ghl_neutrino_rate_channel_bremsstrahlung };
    const double process_scales[ghl_m1_neutrino_pair_process_count]
          = { provider->pair_scale * rho * theta * theta,
              provider->plasmon_scale * theta * theta * theta,
              provider->bremsstrahlung_scale * rho * rho };
    for(int process = 0; process < ghl_m1_neutrino_pair_process_count; ++process) {
      if((mask & process_masks[process]) == 0) {
        continue;
      }
      const double kappa_N = process_scales[process];
      for(int index = ghl_m1_neutrino_nue; index <= ghl_m1_neutrino_anue; ++index) {
        rates[index].eta_N_pair[process] = kappa_N * electron_n_eq_geometric_mean;
        rates[index].eta_E_pair[process]
              = kappa_N * rates[index].mean_energy * rates[index].J_eq;
      }
    }
  }
}

static ghl_error_codes_t m1_test_reference_recover(
      const m1_test_reference_provider_context *restrict provider,
      const ghl_primitive_quantities *restrict prims,
      const bool physical_thermo_valid,
      const double rho,
      const double T,
      const double mu_e,
      const double mu_p,
      const double mu_n,
      const double X_n,
      const double X_p,
      const ghl_error_codes_t original_error,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {
  if(diagnostics != NULL) {
    diagnostics->failures++;
    diagnostics->last_error = original_error;
  }
  if(!physical_thermo_valid
     || provider->failure_policy == ghl_neutrino_rate_failure_return_error) {
    return original_error;
  }
  m1_test_reference_provider_context target_provider = *provider;
  target_provider.channel_mask = 0;
  ghl_m1_neutrino_rates candidate[ghl_m1_neutrino_species_count];
  m1_test_reference_compute_rates(
        &target_provider, rho, T, mu_e, mu_p, mu_n, X_n, X_p, candidate);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    if(ghl_m1_validate_neutrino_rates(&candidate[species], NULL) != ghl_success) {
      return original_error;
    }
  }
  if(provider->failure_policy == ghl_neutrino_rate_failure_equilibrium) {
    for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
      candidate[species].kappa_a_N = provider->equilibrium_recovery_rate;
      candidate[species].kappa_a_E = provider->equilibrium_recovery_rate;
      candidate[species].kappa_tr = candidate[species].kappa_a_E;
      candidate[species].eta_N = candidate[species].kappa_a_N * candidate[species].n_eq;
      candidate[species].eta_E = candidate[species].kappa_a_E * candidate[species].J_eq;
      if((provider->channel_mask & ghl_neutrino_rate_channel_charged_current) != 0
         && candidate[species].lepton_weight != 0.0) {
        candidate[species].kappa_a_N_cc = provider->equilibrium_recovery_rate;
        candidate[species].eta_N_cc
              = candidate[species].kappa_a_N_cc * candidate[species].n_eq;
      }
    }
  }
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    if(ghl_m1_validate_neutrino_rates(&candidate[species], NULL) != ghl_success) {
      return original_error;
    }
  }
  memcpy(rates, candidate, sizeof(candidate));
  if(diagnostics != NULL) {
    if(provider->failure_policy == ghl_neutrino_rate_failure_transparent) {
      diagnostics->transparent_recoveries++;
      diagnostics->last_recovery = ghl_neutrino_rate_recovery_transparent;
    }
    else {
      diagnostics->equilibrium_recoveries++;
      diagnostics->last_recovery = ghl_neutrino_rate_recovery_equilibrium;
    }
  }
  (void)prims;
  return original_error;
}

ghl_error_codes_t ghl_neutrino_rate_provider_initialize_reference(
      ghl_neutrino_rate_provider_context *restrict provider) {
  if(provider == NULL) {
    return ghl_error_m1_null_pointer;
  }
  *provider = (ghl_neutrino_rate_provider_context){ 0 };
  provider->channel_mask = ghl_neutrino_rate_channel_charged_current
                           | ghl_neutrino_rate_channel_nucleon_scattering;
  provider->failure_policy = ghl_neutrino_rate_failure_return_error;
  provider->table_bounds_policy = ghl_neutrino_rate_table_bounds_return_error;
  provider->nu_x_multiplicity = 4.0;
  provider->temperature_code_to_mev = 1.0;
  provider->baryon_mass_code = 1.0;
  provider->charged_current_scale = 1.0e-2;
  provider->scattering_scale = 5.0e-3;
  provider->pair_scale = 1.0e-4;
  provider->bremsstrahlung_scale = 1.0e-4;
  provider->plasmon_scale = 1.0e-5;
  provider->min_mean_energy = 1.0e-12;
  provider->equilibrium_recovery_rate = 1.0;
  return ghl_success;
}

void ghl_neutrino_rate_provider_cache_initialize(
      ghl_neutrino_rate_provider_cache *restrict cache) {
  if(cache != NULL) {
    *cache = (ghl_neutrino_rate_provider_cache){ 0 };
  }
}

ghl_error_codes_t ghl_neutrino_rate_provider_compute_cell(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_cache *restrict cache,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict unused_eos,
      const ghl_primitive_quantities *restrict prims,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]) {
  if(provider == NULL || prims == NULL || rates == NULL) {
    return ghl_error_m1_null_pointer;
  }
  if(diagnostics != NULL) {
    diagnostics->last_recovery = ghl_neutrino_rate_recovery_none;
  }
  const ghl_error_codes_t context_error = m1_test_reference_validate_context(provider);
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

  m1_test_reference_provider_cache staged_cache;
  m1_test_reference_provider_cache *work_cache = NULL;
  if(cache != NULL) {
    staged_cache = *cache;
    work_cache = &staged_cache;
  }
  double rho = 0.0, T = 0.0, Ye = 0.0;
  double muhat = 0.0, mu_e = 0.0, mu_p = 0.0, mu_n = 0.0;
  double X_n = 0.0, X_p = 0.0;
  ghl_error_codes_t error = m1_test_reference_compute_thermo(
        provider, work_cache, unused_eos, prims, &rho, &T, &Ye, &muhat, &mu_e, &mu_p,
        &mu_n, &X_n, &X_p);
  if(error != ghl_success) {
    return m1_test_reference_recover(
          provider, prims, false, rho, T, mu_e, mu_p, mu_n, X_n, X_p, error, diagnostics,
          rates);
  }

  if(cache != NULL && m1_test_reference_same_provenance(cache, provider, unused_eos)
     && cache->rates_valid && cache->rho == rho && cache->T == T && cache->Ye == Ye) {
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
  m1_test_reference_compute_rates(
        provider, rho, T, mu_e, mu_p, mu_n, X_n, X_p, candidate_rates);
  for(int species = 0; species < ghl_m1_neutrino_species_count; ++species) {
    error = ghl_m1_validate_neutrino_rates(&candidate_rates[species], NULL);
    if(error != ghl_success) {
      return m1_test_reference_recover(
            provider, prims, true, rho, T, mu_e, mu_p, mu_n, X_n, X_p, error,
            diagnostics, rates);
    }
  }
  if(cache != NULL) {
    staged_cache.rates_valid = true;
    staged_cache.rho = rho;
    staged_cache.T = T;
    staged_cache.Ye = Ye;
    memcpy(staged_cache.rates, candidate_rates, sizeof(staged_cache.rates));
    staged_cache.provider_snapshot = *provider;
    staged_cache.eos_snapshot = unused_eos;
    *cache = staged_cache;
  }
  memcpy(rates, candidate_rates, sizeof(candidate_rates));
  if(diagnostics != NULL) {
    memset(
          diagnostics->beta_kirchhoff_relative_mismatch, 0,
          sizeof(diagnostics->beta_kirchhoff_relative_mismatch));
    memset(
          diagnostics->beta_kirchhoff_mismatch_valid, 0,
          sizeof(diagnostics->beta_kirchhoff_mismatch_valid));
  }
  return ghl_success;
}
