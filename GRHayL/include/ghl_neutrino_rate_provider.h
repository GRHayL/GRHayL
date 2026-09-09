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
 * formula source. The default backend is a deterministic reference/test
 * provider. The explicit NRPyLeakage initializer selects the production
 * Ruffert backend. Both accept only nu_x_multiplicity == 4 and return an
 * already-summed nu_x bundle; hosts must not multiply the returned
 * heavy-flavor state or exchange a second time. Electron-flavor pair,
 * plasmon, and bremsstrahlung emission is carried in the process-indexed
 * fields of ghl_m1_neutrino_rates, while its scalar absorption/emissivity
 * fields contain charged-current and scattering contributions only.
 * For electron-flavor face transport, add the partner-dependent inverse pair
 * energy opacity to scalar kappa_tr; see Radiation/PAIR_SOURCE_MODEL.md.
 */

typedef enum {
  ghl_neutrino_rate_channel_charged_current    = 1 << 0,
  ghl_neutrino_rate_channel_nucleon_scattering = 1 << 1,
  ghl_neutrino_rate_channel_pair               = 1 << 2,
  ghl_neutrino_rate_channel_bremsstrahlung     = 1 << 3,
  ghl_neutrino_rate_channel_plasmon            = 1 << 4
} ghl_neutrino_rate_channel_t;

typedef enum {
  ghl_neutrino_rate_failure_abort = 0,
  ghl_neutrino_rate_failure_transparent,
  ghl_neutrino_rate_failure_equilibrium,
  ghl_neutrino_rate_failure_hold_last
} ghl_neutrino_rate_failure_policy_t;

typedef enum {
  ghl_neutrino_rate_recovery_none = 0,
  ghl_neutrino_rate_recovery_transparent,
  ghl_neutrino_rate_recovery_equilibrium,
  ghl_neutrino_rate_recovery_hold_last
} ghl_neutrino_rate_recovery_status_t;

typedef enum {
  ghl_neutrino_rate_table_bounds_abort = 0,
  ghl_neutrino_rate_table_bounds_clamp
} ghl_neutrino_rate_table_bounds_policy_t;

typedef enum {
  ghl_neutrino_rate_backend_reference = 0,
  ghl_neutrino_rate_backend_nrpyleakage = 1
} ghl_neutrino_rate_backend_t;

typedef struct {
  bool use_tabulated_eos;
  int channel_mask;
  ghl_neutrino_rate_failure_policy_t failure_policy;
  ghl_neutrino_rate_table_bounds_policy_t table_bounds_policy;
  double nu_x_multiplicity;
  /**
   * Host-managed identity for mutable EOS/table contents.  Initialize this
   * field through ghl_neutrino_rate_provider_initialize_default() and advance
   * it after every in-place mutation of the EOS object or its table storage.
   */
  uint64_t eos_generation;
  /* Backend configuration is validated exactly. The reference backend uses
   * its documented analytic scales; the production backend fixes every scale
   * to 1 because the raw kernel already supplies physical channel rates. */
  double rho_code_to_cgs;
  double temperature_code_to_mev;
  double opacity_cgs_to_code;
  double emissivity_cgs_to_code;
  double baryon_mass_code;
  double charged_current_scale;
  double scattering_scale;
  double pair_scale;
  double bremsstrahlung_scale;
  double plasmon_scale;
  double min_mean_energy;
  /** Backend selection. */
  ghl_neutrino_rate_backend_t backend;
  /** Degraded equilibrium-recovery rate in inverse code-time units. */
  double equilibrium_recovery_rate;
} ghl_neutrino_rate_provider_context;

typedef struct {
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

typedef struct {
  int failures;
  int table_bound_hits;
  int cache_hits;
  int cache_misses;
  int clamped_inputs;
  int transparent_recoveries;
  int equilibrium_recoveries;
  int hold_last_recoveries;
  int active_channel_mask;
  ghl_error_codes_t last_error;
  /** Recovery used by the most recent call; none for an ordinary result. */
  ghl_neutrino_rate_recovery_status_t last_recovery;
  /** Last successfully published production-call consistency diagnostic. */
  double beta_kirchhoff_relative_mismatch[ghl_m1_neutrino_species_count];
  bool beta_kirchhoff_mismatch_valid[ghl_m1_neutrino_species_count];
} ghl_neutrino_rate_provider_diagnostics;

/**
 * Initialize the deterministic, table-free reference provider. This choice
 * is independent of whether HDF5 support was compiled in; use the explicit
 * NRPyLeakage initializer for the production, table-backed provider.
 *
 * @param provider Caller-owned context to initialize. Existing contents are
 *        replaced on success.
 * @return @c ghl_success on initialization, or
 *         @c ghl_error_m1_null_pointer when @p provider is NULL.
 */
ghl_error_codes_t ghl_neutrino_rate_provider_initialize_default(
      ghl_neutrino_rate_provider_context *restrict provider);

/**
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

/**
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

/**
 * Compute one frozen rate bundle per evolved species.
 *
 * Both built-in backends accept only nu_x_multiplicity == 4 and return an
 * already-summed heavy-lepton bundle. Electron-flavor pair, plasmon, and
 * bremsstrahlung number and energy emissivities are retained independently in
 * the process-indexed fields of each electron-flavor result. Any other finite
 * or nonfinite multiplicity is rejected before cache lookup, thermodynamic
 * reconstruction, or rate calculation; recovery policies do not mask this
 * boundary error.
 *
 * The caller-owned cache is one transactional record with separate recovered
 * thermodynamic and final-rate keys. Cache hits and hold-last recovery require
 * an exact primitive key, complete provider-context snapshot, EOS pointer, and
 * eos_generation match. No validity flag or cached value is committed until a
 * complete rate bundle has passed validation.
 *
 * @param provider Initialized provider context; it is read-only for this call.
 * @param cache Optional caller-owned cache. Pass NULL to compute without cache
 *        reuse; a non-NULL cache must not be shared concurrently.
 * @param diagnostics Optional caller-owned diagnostics record. Counters and
 *        last-call status are updated when non-NULL.
 * @param eos EOS parameters used by the selected backend. The reference
 *        backend may run without an EOS object; the NRPyLeakage backend
 *        requires a compatible tabulated EOS.
 * @param prims Frozen cell primitives. The provider reads density,
 *        temperature, and electron fraction and does not modify them.
 * @param rates Output array indexed by
 *        ghl_m1_neutrino_species_t; all three entries are published only
 *        after the complete bundle passes validation.
 * @return @c ghl_success on a validated bundle or an allowed recovery;
 *         otherwise an error such as @c ghl_error_m1_null_pointer,
 *         @c ghl_error_m1_microphysics_failure, or
 *         @c ghl_error_used_disabled_hdf5. On an error that does not publish
 *         a recovery, @p rates and cached values remain unchanged.
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
