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

The current tabulated public wrapper hard-codes `ghl_eos_table_stellarcollapse`
and `clean_sound_speed = false` before calling `ghl_initialize_tabulated_eos`.

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
beta-equilibrium rho-map helpers, and tabulated enthalpy/sound-speed routing.

EOS-specific HLLE flux pointer declarations live beside the EOS function
pointers as `ghl_calculate_HLLE_fluxes_dirn0`,
`ghl_calculate_HLLE_fluxes_dirn1`, and `ghl_calculate_HLLE_fluxes_dirn2`.
Concrete Flux_Source routines are split into hybrid, hybrid-entropy,
tabulated, and tabulated-entropy families. Current core EOS initialization does
not assign those flux pointers; route flux questions to [Flux Source](../flux-source.md)
and `GRHayL/include/ghl_flux_source.h`.

Sources: `GRHayL/GRHayL_Core/initialize_eos.c`,
`GRHayL/include/ghl_eos_functions.h`,
`GRHayL/include/ghl_eos_functions_declaration.h`,
`GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c`,
`GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`,
`GRHayL/include/ghl_flux_source.h`.

## HDF5 Build And Runtime Contract

Default configured builds expect HDF5 support for tabulated EOS. Passing
`--disable-hdf5` to `configure` defines `GHL_DISABLE_HDF5` and filters
tabulated/HDF5 source and test paths out of the generated build. Downstream or
manual builds that bypass `configure` must mirror both parts of that contract:
define `GHL_DISABLE_HDF5` and omit tabulated/HDF5 implementation sources.

When `GHL_DISABLE_HDF5` is defined, tabulated EOS runtime paths are disabled.
`ghl_initialize_tabulated_eos` and
`ghl_initialize_tabulated_eos_functions_and_params` return
`ghl_error_used_disabled_hdf5`; `ghl_initialize_eos_functions` on
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
- `plan_synth.md`
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
