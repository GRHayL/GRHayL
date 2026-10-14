#ifndef GHL_NEUTRINO_RATE_PROVIDER_H_
#define GHL_NEUTRINO_RATE_PROVIDER_H_

#include "ghl_m1.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @addtogroup Radiation
 *  @{ */

/*
 * The provider owns EOS/table lookup and channel-specific microphysics.
 * Radiation consumes frozen ghl_m1_neutrino_rates bundles and validates their
 * contract; it does not treat this interface as a direct production-rate
 * formula source. Every initializer declared here selects the production
 * Ruffert backend, and all of them accept only nu_x_multiplicity == 4 and
 * return an already-summed nu_x bundle; hosts must not multiply the returned
 * heavy-flavor state or exchange a second time. For electron flavors, scalar
 * absorption and emissivity fields contain charged-current contributions only;
 * isoenergetic nucleon scattering is separate in kappa_s, and kappa_tr is the
 * transport sum kappa_a_E + kappa_s. Electron-flavor pair, plasmon, and
 * bremsstrahlung emission is carried in the process-indexed fields of
 * ghl_m1_neutrino_rates. For nu_x, the scalar coefficients include the
 * enabled aggregate pair, plasmon, and bremsstrahlung contributions. For
 * electron-flavor face transport, add the partner-dependent inverse pair
 * energy opacity to scalar kappa_tr; see docs/raw/Radiation_pair_source_model.md.
 */

typedef enum {
  ghl_neutrino_rate_channel_charged_current = 1 << 0,
  ghl_neutrino_rate_channel_nucleon_scattering = 1 << 1,
  ghl_neutrino_rate_channel_pair = 1 << 2,
  ghl_neutrino_rate_channel_bremsstrahlung = 1 << 3,
  ghl_neutrino_rate_channel_plasmon = 1 << 4
} ghl_neutrino_rate_channel_t;

typedef enum {
  ghl_neutrino_rate_failure_return_error = 0,
  ghl_neutrino_rate_failure_transparent,
  ghl_neutrino_rate_failure_equilibrium
} ghl_neutrino_rate_failure_policy_t;

typedef enum {
  ghl_neutrino_rate_recovery_none = 0,
  ghl_neutrino_rate_recovery_transparent,
  ghl_neutrino_rate_recovery_equilibrium
} ghl_neutrino_rate_recovery_status_t;

typedef enum {
  ghl_neutrino_rate_table_bounds_return_error = 0,
  ghl_neutrino_rate_table_bounds_clamp
} ghl_neutrino_rate_table_bounds_policy_t;

/** @ingroup m1_rates */
typedef struct ghl_neutrino_rate_provider_context {
  /**
   * Enabled weak-interaction channels. Production equilibrium moments and
   * beta/Kirchhoff diagnostics are computed independently of this mask;
   * disabled channels contribute zero rates.
   */
  int channel_mask;
  ghl_neutrino_rate_failure_policy_t failure_policy;
  ghl_neutrino_rate_table_bounds_policy_t table_bounds_policy;
  double nu_x_multiplicity;
  /**
   * Host-managed identity for mutable EOS/table contents.  Initialize this
   * field through a provider initializer and advance
   * it after every in-place mutation of the EOS object or its table storage.
   */
  uint64_t eos_generation;
  /** Recovery absorption rate in inverse code-time units. */
  double equilibrium_recovery_rate;
} ghl_neutrino_rate_provider_context;

/** @ingroup m1_rates */
typedef struct ghl_neutrino_rate_provider_cache {
  bool thermo_valid;
  bool rates_valid;
  /* Recovered, validated thermodynamic key. */
  double thermo_rho;
  double thermo_T;
  double thermo_Ye;
  /* Final-rate key, intentionally distinct from the thermo key. */
  double rho;
  double T;
  double Ye;
  double muhat;
  double mu_e;
  double mu_p;
  double mu_n;
  double X_n;
  double X_p;
  ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count];
  ghl_neutrino_rate_provider_context provider_snapshot;
  const ghl_eos_parameters *eos_snapshot;
  double beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_species_count];
  bool beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_species_count];
} ghl_neutrino_rate_provider_cache;

/** @ingroup m1_rates */
typedef struct ghl_neutrino_rate_provider_diagnostics {
  int failures;
  int table_bound_hits;
  int cache_hits;
  int cache_misses;
  int clamped_inputs;
  int transparent_recoveries;
  int equilibrium_recoveries;
  int active_channel_mask;
  ghl_error_codes_t last_error;
  /** Recovery used by the most recent call; none for an ordinary result. */
  ghl_neutrino_rate_recovery_status_t last_recovery;
  /** Last successfully published production-call consistency diagnostic. */
  double beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_species_count];
  bool beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_species_count];
} ghl_neutrino_rate_provider_diagnostics;

/** @ingroup m1_rates
 *
 * Initialize the default provider, which is the production, table-backed
 * NRPyLeakage provider; see ghl_neutrino_rate_provider_initialize_nrpyleakage().
 * It is an error, not a silent substitution, to request it without HDF5.
 *
 * @param provider Caller-owned context to initialize. Existing contents are
 *        replaced on success.
 * @return @c ghl_success when the table-backed context is available;
 *         @c ghl_error_m1_null_pointer for NULL @p provider; or
 *         @c ghl_error_used_disabled_hdf5 when HDF5 support is disabled.
 */
ghl_error_codes_t ghl_neutrino_rate_provider_initialize_default(
      ghl_neutrino_rate_provider_context *restrict provider);

/** @ingroup m1_rates
 *
 * Initialize a caller-owned provider cache before its first use.
 *
 * The initializer clears all cached state, sets thermo_valid and rates_valid
 * false, and clears the cached EOS/provider provenance. A cache may then be
 * reused across calls to ghl_neutrino_rate_provider_compute_cell(); the
 * provider context, EOS pointer/generation, and primitive keys determine
 * whether a cached value is reusable. Reinitialize the cache when its
 * contents are intentionally discarded. The cache must not be shared by
 * concurrent calls without external synchronization. Passing NULL is a
 * no-op.
 *
 * @param cache Caller-owned cache to clear, or NULL to do nothing.
 */
void ghl_neutrino_rate_provider_cache_initialize(
      ghl_neutrino_rate_provider_cache *restrict cache);

/** @ingroup m1_rates
 *
 * Initialize the production, table-backed NRPyLeakage rate provider.
 *
 * @param provider Caller-owned context to initialize. Existing contents are
 *        replaced on success.
 * @return @c ghl_success when the table-backed context is available;
 *         @c ghl_error_m1_null_pointer for NULL @p provider; or
 *         @c ghl_error_used_disabled_hdf5 when HDF5 support is disabled.
 */
ghl_error_codes_t ghl_neutrino_rate_provider_initialize_nrpyleakage(
      ghl_neutrino_rate_provider_context *restrict provider);

/** @ingroup m1_rates
 *
 * Compute one frozen rate bundle per evolved species.
 *
 * The production provider accepts only nu_x_multiplicity == 4 and returns an
 * already-summed heavy-lepton bundle. Electron-flavor pair, plasmon, and
 * bremsstrahlung number and energy emissivities are retained independently in
 * the process-indexed fields of each electron-flavor result. Any other finite
 * or nonfinite multiplicity is rejected before cache lookup, thermodynamic
 * reconstruction, or rate calculation; recovery policies do not mask this
 * boundary error.
 *
 * The caller-owned cache is one transactional record with separate recovered
 * thermodynamic and final-rate keys. Cache hits require an exact primitive
 * key, complete provider-context snapshot, EOS pointer, and eos_generation
 * match. No validity flag or cached value is committed until a complete rate
 * bundle has passed validation.
 *
 * @param provider Initialized provider context; it is read-only for this call.
 * @param cache Optional caller-owned cache. Pass NULL to compute without cache
 *        reuse; a non-NULL cache must not be shared concurrently.
 * @param diagnostics Optional caller-owned diagnostics record. Counters and
 *        last-call status are updated when non-NULL.
 * @param eos Tabulated EOS parameters with registered dispatch pointers for
 *        the configured table type.
 * @param prims Frozen cell primitives. The provider reads density,
 *        temperature, and electron fraction and does not modify them.
 * @param rates Output array indexed by
 *        ghl_m1_neutrino_species_t; all three entries are published only
 *        after the complete bundle passes validation.
 * @return @c ghl_success only on a validated bundle; otherwise an error such
 *         as @c ghl_error_m1_null_pointer,
 *         @c ghl_error_m1_microphysics_failure, or
 *         @c ghl_error_used_disabled_hdf5. When a recovery policy publishes a
 *         bundle into @p rates, the original error is still returned and
 *         diagnostics->last_recovery identifies the recovery; a host that
 *         accepts recovery must check that field. If physical thermodynamic
 *         targets are unavailable, @p rates and cached values remain
 *         unchanged.
 */
ghl_error_codes_t ghl_neutrino_rate_provider_compute_cell(
      const ghl_neutrino_rate_provider_context *restrict provider,
      ghl_neutrino_rate_provider_cache *restrict cache,
      ghl_neutrino_rate_provider_diagnostics *restrict diagnostics,
      const ghl_eos_parameters *restrict eos,
      const ghl_primitive_quantities *restrict prims,
      ghl_m1_neutrino_rates rates[ghl_m1_neutrino_species_count]);

#ifdef __cplusplus
}
#endif

/** @} */

#endif // GHL_NEUTRINO_RATE_PROVIDER_H_
