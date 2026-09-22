# EOS Initialization And Dispatch

## Purpose

This page maps how GRHayL initializes `ghl_eos_parameters`, installs EOS
function pointers, and routes EOS-dependent calls. It is a guide into the
repo-local source; source, headers, Doxygen source, tests, and `configure`
remain authoritative.

## Public Wrapper Versus Low-Level Initialization

Use the public wrappers when setting up an EOS for normal GRHayL use:

- `ghl_initialize_simple_eos_functions_and_params`
- `ghl_initialize_hybrid_eos_functions_and_params`
- `ghl_initialize_tabulated_eos_functions_and_params`

Those wrappers call `ghl_initialize_eos_functions` first, then call the matching
low-level parameter initializer. The low-level functions
`ghl_initialize_simple_eos`, `ghl_initialize_hybrid_eos`, and
`ghl_initialize_tabulated_eos` populate `ghl_eos_parameters` fields but do not,
by themselves, install the global function-pointer dispatch layer.
Low-level callers must therefore install the matching dispatch family first;
parameter construction and tabulated rollback use those callbacks.

Function-pointer dispatch is process-wide global state defined through
`GRHayL/include/ghl_eos_functions_declaration.h`; it is not stored per EOS
object. Initializing a different EOS family replaces general pointers such as
`ghl_compute_h_and_cs2` and `ghl_con2prim_multi_method` for all callers.

Sources: `GRHayL/GRHayL_Core/initialize_eos.c`, `GRHayL/include/ghl.h`,
`docs/raw/EOS.dox`, `docs/raw/GRHayL_Core.dox`.

## EOS Parameter Setup

`ghl_eos_parameters` stores the selected `eos_type`, optional `table_type`,
atmosphere/floor/ceiling values, hybrid piecewise-polytrope arrays, tabulated
table state, and root-finding settings. The enum choices are
`ghl_eos_simple`, `ghl_eos_hybrid`, and `ghl_eos_tabulated`.

Simple EOS setup uses the hybrid implementation path with a one-piece
ideal-fluid setup. `ghl_initialize_simple_eos` sets `eos_type` to
`ghl_eos_simple`, sets `neos = 1`, uses the input `Gamma` for both
`Gamma_th` and `Gamma_ppoly[0]`, sets one-piece hybrid constants, stores
pressure atmosphere/floor/ceiling values from inputs or defaults, computes
epsilon values from the ideal-fluid pressure relation, computes entropy through
the hybrid entropy helper, and sets `tau_atm = rho_atm * eps_atm`.

Simple and hybrid initializers set the unused-family `Y_e_atm` and `T_atm`
placeholders to zero. `ghl_set_prims_to_constant_atm` therefore copies defined
composition and temperature values for every successfully initialized family;
callers may still override them afterward.

Hybrid EOS setup sets `eos_type` to `ghl_eos_hybrid`, stores the requested
piece count and piece arrays, computes derived `K_ppoly`,
`eps_integ_const`, `p_ppoly`, atmosphere values, floors, and ceilings through
hybrid helpers.

Tabulated EOS setup sets `eos_type` to `ghl_eos_tabulated`, stores
`table_type` and `clean_sound_speed`, reads the table, clamps requested
`rho`, `Y_e`, and `T` bounds to table bounds, initializes atmosphere pressure,
energy, and entropy through tabulated interpolation, sets table-derived
pressure/energy/entropy bounds, initializes beta-equilibrium arrays to `NULL`,
and sets `root_finding_precision = 1e-10`.

Simple and hybrid setup construct a local candidate and leave the destination
unchanged on error. Tabulated setup instead requires a destination with no live
owned allocations: success transfers the completed candidate; failure releases
it and publishes an empty aggregate tagged `ghl_eos_tabulated`. A live table
must be cleaned before reinitialization or an EOS-family switch. Global dispatch
selection is process-wide and is not rolled back when parameter setup fails.

Tabulated initialization sets `tau_atm = rho_min * eps_min`, unlike the
simple/hybrid `rho_atm * eps_atm` assignment. `ghl.h` describes `_atm` fields
as atmosphere values and Con2Prim consumes `tau_atm` as a floor, but repo-local
evidence does not resolve whether this family difference is intended.

The current tabulated public wrapper hard-codes `ghl_eos_table_stellarcollapse`
and both `clean_sound_speed = false` and `enable_neural_net_c2p = false` before
calling `ghl_initialize_tabulated_eos`.

## Caller-Visible Validation Boundary

- Simple/hybrid atmosphere density must be finite and positive. Simple
  atmosphere pressure must be finite and nonnegative. Scalar inputs and used
  derived atmosphere/bound fields must be finite. Cold-piece gamma values reject
  exact `0` and `1`; thermal gamma rejects exact `1`. Other finite sub-unity and
  negative values remain permitted by the Core API.
- Negative density minima normalize to zero. Simple keeps its independently
  configured pressure floor and uses `-DBL_MAX` for energy/entropy minima at a
  zero density floor. Hybrid uses `-DBL_MAX` for pressure, energy, and entropy
  minima there. Negative density or simple-pressure maxima normalize to the
  `1e300` disabled marker; simple derived maxima use `DBL_MAX` if either marker
  is present, while hybrid uses `DBL_MAX` when density is unbounded.
- Hybrid setup enforces `1 <= neos <= MAX_EOS_PARAMS`, requires one finite
  nonsingular gamma per piece, and requires `neos - 1` finite, positive,
  strictly increasing breakpoints. The breakpoint pointer may be `NULL` for a
  one-piece EOS. Constructed `K_ppoly` values and integration constants must be
  finite. A zero `K_ppoly0` intentionally produces a zero cold-pressure curve;
  otherwise every constructed `K_ppoly` must remain nonzero. The zero curve has
  no unique cold-pressure density inverse. Consumed pressure breakpoints must
  be finite.
- Tabulated initialization rejects nonfinite atmosphere and requested bounds
  before reading the table. Requested min/max values are then clamped to
  table bounds, but atmosphere `(rho,Y_e,T)` is not clamped before its direct
  interpolation call and therefore must already lie inside table bounds. A
  final finite-state gate covers published atmosphere values, bounds, derived
  extrema, and `tau_atm`.
- `ghl_initialize_eos_functions` returns `void` and has no final unknown-type
  branch. Pass only a declared `ghl_eos_t`; an invalid value still installs the
  hybrid helper family, leaves `ghl_con2prim_multi_method` at its prior value,
  and does not report an error.

These are observed checks, not an endorsement of missing validation.

Sources: `GRHayL/include/ghl.h`,
`GRHayL/GRHayL_Core/initialize_eos.c`.

## Dispatch Outcomes

`ghl_initialize_eos_functions` always initializes the hybrid EOS pointer family
through `NRPyEOS_initialize_hybrid_functions`. When HDF5 is enabled, it also
initializes the tabulated pointer family through
`NRPyEOS_initialize_tabulated_functions`.

For `ghl_eos_simple` and `ghl_eos_hybrid`,
`ghl_initialize_eos_functions` routes:

- `ghl_compute_h_and_cs2` to `NRPyEOS_hybrid_compute_enthalpy_and_cs2`
- `ghl_con2prim_multi_method` to `ghl_con2prim_hybrid_multi_method`

For `ghl_eos_tabulated` with HDF5 enabled, it routes:

- `ghl_compute_h_and_cs2` to `NRPyEOS_tabulated_compute_enthalpy_and_cs2`
- `ghl_con2prim_multi_method` to `ghl_con2prim_tabulated_multi_method`

The hybrid pointer initializer installs the hybrid helper family used for
piece lookup, cold pressure and energy, entropy, epsilon, rho bounds, and
enthalpy/sound-speed routing. The tabulated pointer initializer installs the
table read/free routines, interpolation families, table-bound enforcement,
beta-equilibrium rho-map helpers, beta-equilibrium cleanup, and tabulated
enthalpy/sound-speed routing. Use
[tabulated interpolator catalog](tabulated-interpolator-catalog.md) for the
exact registry seam.

Concrete Flux_Source routines are split into hybrid, hybrid-entropy,
tabulated, and tabulated-entropy families. Unsuffixed generic HLLE globals
remain only as uninitialized compatibility storage; new callers select a direct
family/direction/entropy variant. Route flux questions to [Flux Source](../flux-source.md) and
`GRHayL/include/ghl_flux_source.h`.

Sources: `GRHayL/GRHayL_Core/initialize_eos.c`,
`GRHayL/include/ghl_eos_functions.h`,
`GRHayL/include/ghl_eos_functions_declaration.h`,
`GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c`,
`GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`,
`GRHayL/include/ghl_flux_source.h`.

## HDF5 Build And Runtime Contract

Default configured builds expect HDF5 support for tabulated EOS. Passing
`--disable-hdf5` to `configure` defines `GHL_DISABLE_HDF5` and filters most
tabulated/HDF5 implementation and test paths out of the generated build. It
retains a narrow documented set of low-level helpers, disabled-feature stubs,
and direct flux symbols. [Build and CI](../../build-and-ci.md) owns the exact
retained/excluded source and test surface. Downstream or manual builds that
bypass `configure` must reproduce that actual filtered source set, not assume
every path containing `Tabulated` disappears.

When `GHL_DISABLE_HDF5` is defined, tabulated EOS runtime paths are disabled.
`ghl_initialize_tabulated_eos` and
`ghl_initialize_tabulated_eos_functions_and_params` return
`ghl_error_used_disabled_hdf5` and publish an empty tagged aggregate for a
valid destination; `ghl_initialize_eos_functions` on
`ghl_eos_tabulated` calls the disabled-HDF5 error macro; HDF5-only code-error
tests are skipped in no-HDF5 builds.

Sources: `README.md`, `configure`,
`GRHayL/GRHayL_Core/initialize_eos.c`,
`GRHayL/include/ghl_nrpyeos_tabulated.h`,
`Unit_Tests/unit_test_code_error.c`.

## Dependent Areas

Keep dependent content routed, not expanded here:

- [Con2Prim recovery flow](../con2prim/recovery-flow.md) owns recovery order,
  solver selection, and backup behavior.
- [Flux Source](../flux-source.md) owns characteristic speeds and HLLE flux
  implementations.
- [Neutrinos](../neutrinos.md) owns tabulated-EOS leakage dependencies.
- `implementations/GRHayLib/` owns downstream integration details.

## Source-Of-Truth Paths

- `AGENTS.md`
- `wiki/index.md`
- `wiki/catalog.md`
- `wiki/gems/eos.md`
- `GRHayL/GRHayL_Core/initialize_eos.c`
- `GRHayL/include/ghl.h`
- `GRHayL/include/ghl_eos_functions.h`
- `GRHayL/include/ghl_eos_functions_declaration.h`
- `GRHayL/include/ghl_flux_source.h`
- `GRHayL/include/ghl_nrpyeos_tabulated.h`
- `GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c`
- `GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`
- `GRHayL/EOS/Hybrid/make.code.defn`
- `GRHayL/EOS/Tabulated/make.code.defn`
- `README.md`
- `configure`
- `docs/raw/EOS.dox`
- `docs/raw/GRHayL_Core.dox`
- `Unit_Tests/unit_test_code_error.c`
