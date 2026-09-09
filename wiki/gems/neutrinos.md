# Neutrinos Gem

## Purpose

This page routes the legacy NRPyLeakage interfaces and the neutrino M1 rate
provider. It records source, header, build, HDF5, cache, and local-test
boundaries without replacing the repository files that define them.

## Read Order

1. [Physics And EOS Contract](neutrinos/physics-and-eos-contract.md) for species,
   chemical potentials, free-nucleon fractions, units, and source conventions.
2. [Generator Provenance](neutrinos/generator-provenance.md) for the ancestral
   notebooks and the boundary between generated expressions and checked-in C.
3. [API And Data](neutrinos/api-and-data.md) for public structs, declarations,
   errors, constants, and HDF5/EOS dependencies.
4. [Implementation Flow](neutrinos/implementation-flow.md) for the five legacy
   source files, raw-rate bridge, and neutrino M1 provider flow.
5. [Tests And Fixtures](neutrinos/tests-and-fixtures.md) for the checked-in local
   tests and fixture setup.
6. [NRPyLeakage `bns_nurates`-Class Future Work](../future_work/index.md)
   for explicitly future, not-implemented enhancements.

## Ground Truth

- Legacy source manifest: `GRHayL/Neutrinos/NRPyLeakage/make.code.defn`
  lists five C files.
- Neutrino M1 source manifest: `GRHayL/Radiation/Neutrinos/make.code.defn`
  lists thirteen C files; shared E/F helpers are listed by
  `GRHayL/Radiation/make.code.defn`.
- Public aggregate and legacy interfaces: `GRHayL/include/ghl_radiation.h`
  and `GRHayL/include/ghl_nrpyleakage.h`.
- M1 and provider interfaces: `GRHayL/include/ghl_m1.h` and
  `GRHayL/include/ghl_neutrino_rate_provider.h`.
- Legacy leakage tests are
  `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`,
  `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`, and
  `Unit_Tests/unit_test_nrpyleakage_luminosities.c`.
- Doxygen source: `docs/raw/Radiation.dox` owns the neutrino M1 group, while
  `ghl_radiation.h` declares the legacy `Neutrinos` group.

If this page conflicts with current source, headers, manifests, or tests, trust
the repo-local files and narrow the claim. A declaration or local test does not
by itself prove dispatch, workflow execution, downstream integration, or an
external verification result.

## Provider Boundary

`ghl_neutrino_rate_provider_compute_cell` is the boundary between microphysics
and transport. The provider owns EOS/table lookup, unit conversion, channel
selection, equilibrium targets, and production weak-rate calculations. It
returns a validated, frozen `ghl_m1_neutrino_rates` bundle for `nue`, `anue`,
and already-summed `nux`. The transport operators consume that bundle and do
not perform the production weak-rate calculation.

The exact provider names are: `ghl_neutrino_rate_backend_reference` for the
deterministic reference/test backend selected by
`ghl_neutrino_rate_provider_initialize_default` and
`ghl_neutrino_rate_backend_nrpyleakage` for the table-backed production path
selected by `ghl_neutrino_rate_provider_initialize_nrpyleakage`.
Both paths require `nu_x_multiplicity == 4.0` and return an already-summed
heavy-lepton bundle; callers must not multiply it again.

## HDF5, Cache, And Diagnostics

`ghl_nrpyeos_tabulated.h` includes `<hdf5.h>` only when HDF5 is enabled, so
`ghl_radiation.h`, `ghl_m1.h`, and the provider header remain includable in a
`GHL_DISABLE_HDF5` build. Table-backed EOS/provider operations still return
`ghl_error_used_disabled_hdf5`; the production initializer cannot provide its
normal path without HDF5. The reference provider can be used with
`use_tabulated_eos = false` for table-free calls.

The cache and diagnostics objects are caller-owned. Initialize the cache with
`ghl_neutrino_rate_provider_cache_initialize`; diagnostics have no initializer
function and should be zero-initialized (for example,
`ghl_neutrino_rate_provider_diagnostics diagnostics = {0}`). Both objects may
be omitted by passing `NULL` where the API permits. Diagnostics counters
accumulate until the caller resets the object. Cache reuse and hold-last recovery
require an exact primitive/context snapshot, matching EOS pointer, and matching
caller-managed `eos_generation`; increment that generation after every in-place
EOS/table mutation. Provider output and cache state are committed only after a
complete validated bundle succeeds.

## Legacy Public Surface

- Data: `ghl_neutrino_luminosities`, `ghl_neutrino_opacities`, and the
  `ghl_neutrino_optical_depths` alias.
- Functions: `NRPyLeakage_Fermi_Dirac_integrals`,
  `NRPyLeakage_compute_neutrino_opacities`,
  `NRPyLeakage_compute_neutrino_luminosities`,
  `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`, and
  `NRPyLeakage_optical_depths_PathOfLeastResistance`.

The M1 provider's private thermo builder requires an initialized tabulated EOS
and HDF5. Its private
`ghl_m1_nrpyleakage_compute_raw_rates_from_thermo` adapter consumes an already
validated thermo state and does not perform a table lookup. The five-file
`GRHayL/Neutrinos/NRPyLeakage/make.code.defn` manifest is the legacy build
boundary.

## Scope Notes

The current contribution documented through this route is neutrino-only and
library-level. Shared E/F M1 helpers are dependencies of the neutrino path, not
evidence of a separate non-neutrino implementation. Nothing here establishes
downstream host integration or an external verification campaign.

## Future Work

The [NRPyLeakage `bns_nurates`-Class Future Work](../future_work/index.md)
folder remains a to-do catalog. It does not add dependencies, source formulas,
build entries, or validated capabilities to the current five-file legacy model
or the current neutrino M1 provider boundary.

## Ground Truth References

- [Original Tabulated_EOS_IllinoisGRMHD repository](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD)
