# Radiation M1 API And Build Boundary

Use this page to distinguish the installed neutrino M1/provider surface from
private helpers and source files that are not part of the production manifest.
The headers, source files, and make definitions named below are the ground
truth.

## Installed public headers

`GRHayL/include/make.code.defn` installs these Radiation-facing headers:

| Header | Role |
| --- | --- |
| `GRHayL/include/ghl_m1.h` | Shared M1 data, closure/repair/moment/flux/source helpers, and neutrino state, rates, exchange, number, lepton, and local implicit APIs. |
| `GRHayL/include/ghl_neutrino_rate_provider.h` | Provider context, channel and failure policies, cache, diagnostics, backend initialization, and frozen cell-rate computation. |
| `GRHayL/include/ghl_radiation.h` | Aggregate header for core, M1, provider, and legacy NRPyLeakage radiation declarations. |
| `GRHayL/include/ghl_nrpyleakage.h` | Legacy leakage data and public leakage routines. The M1 raw-rate adapter is private to Radiation. |

`ghl_m1_utils.h` and `ghl_m1_neutrino_implicit.h` are private build inputs
listed through `INCS`; they are not installed as public headers.

The installed `ghl_radiation.h` and `ghl_neutrino_rate_provider.h` surfaces do
not pull in the tabulated-EOS/HDF5 header. External code that needs the
tabulated-EOS API must include `ghl_nrpyeos_tabulated.h` explicitly. The
HDF5-enabled provider translation unit includes that header privately for its
`NRPyEOS_*` calls; no internal build macro or transitive public-header
declaration is required.

## Compiled source lists

`GRHayL/make.code.defn` includes the `Radiation` gem. The shared manifest
`GRHayL/Radiation/make.code.defn` contains the shared closure, moments,
source-geometry, wave-speed, prepared-face transport, stress-energy,
diagnostics, Newton, and utility implementations used by the neutrino path.

The nested manifest `GRHayL/Radiation/Neutrinos/make.code.defn` contains the
provider, rate validation, repair, number flux, prepared-face wrapper, source,
lepton-exchange, and local implicit-solve files.

The component-wise Rusanov arithmetic required by the prepared neutrino
transport operation is compiled by `GRHayL/Flux_Source/make.code.defn` from
`GRHayL/Flux_Source/ghl_calculate_Rusanov_flux.c`.

The Radiation/Neutrinos manifest includes the private
`ghl_m1_nrpyleakage_kernel.c` adapter. It provides the EOS-to-thermodynamic-
state boundary and optical-depth-independent Ruffert channel rates for the
M1 provider. The legacy NRPyLeakage manifest and its generated source files
remain unchanged.

## Ownership and call boundary

- The rate provider owns EOS/table lookup, channel policy, caching, unit
  conversion, and validation. It publishes frozen
  `ghl_m1_neutrino_rates` bundles.
- M1 pointwise kernels consume the frozen bundle and do not perform an EOS
  lookup. They operate on undensitized states and leave host mesh/storage
  policy to the caller.
- The private M1 NRPyLeakage adapter owns the thermodynamic state and
  channel-rate calculation used by the provider. Those types and calls are not
  part of the public NRPyLeakage API.
- The legacy leakage APIs retain their existing opacity, luminosity,
  optical-depth, and GRMHD source contracts.

## Production boundary

Only the files listed by the active make manifests define the supported
production build surface described here. The neutrino M1 implementation has one
primary pointwise numerical path: metric light-cone speeds, the linear
`E^2(1-epsilon)/F^2` realizability rescale, the full four-dimensional closure,
and four-point blended Rusanov transport without a separate diffusion
correction. A finite full-four-dimensional closure candidate that fails the
physical PSD check is replaced by the built-in Eulerian Minerbo admissibility
fallback and marked in the returned closure; this is not a user-selectable
alternative. Its local implicit source update uses the existing
finite-difference/Newton solver. Host mesh loops, reconstruction, boundaries,
schedules, AMR, matter recovery, and the coupled limiter are intentionally
outside these manifests. Legacy compatibility fields and rejected policy
arguments do not promote discarded numerical alternatives. Callback and
opacity experiments, discrete-operator artifacts, downstream host
implementations, and unit-test content are not silently promoted by this
boundary.

The installed-header manifest already lists `GRHayL/include/ghl_m1.h`; no
additional private M1 header is installed for these routes.

## Ground truth references

- `GRHayL/make.code.defn`
- `GRHayL/Radiation/make.code.defn`
- `GRHayL/Radiation/Neutrinos/make.code.defn`
- `GRHayL/Flux_Source/make.code.defn`
- `GRHayL/Flux_Source/ghl_calculate_Rusanov_flux.c`
- `GRHayL/Neutrinos/NRPyLeakage/make.code.defn`
- `GRHayL/Radiation/Neutrinos/ghl_m1_nrpyleakage_kernel.c`
- `GRHayL/Radiation/Neutrinos/ghl_m1_nrpyleakage_kernel.h`
- `GRHayL/include/make.code.defn`
- `GRHayL/include/ghl_m1.h`
- `GRHayL/include/ghl_neutrino_rate_provider.h`
- `GRHayL/include/ghl_nrpyleakage.h`
