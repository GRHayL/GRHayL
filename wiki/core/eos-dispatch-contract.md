# Core EOS Dispatch Contract

## Purpose

This page routes Core-owned EOS initialization wrappers and dispatch behavior.
It separates the `GRHayL/GRHayL_Core/initialize_eos.c` contract from EOS gem
implementation details under `GRHayL/EOS/`; use
[EOS initialization and dispatch](../gems/eos/initialization-and-dispatch.md)
and child EOS pages for implementation internals.

Ground truth:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c),
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h),
[`docs/raw/GRHayL_Core.dox`](../../docs/raw/GRHayL_Core.dox), and
[`docs/raw/EOS.dox`](../../docs/raw/EOS.dox).

## Public Wrapper Boundary

Core owns the public wrapper layer:
`ghl_initialize_simple_eos_functions_and_params`,
`ghl_initialize_hybrid_eos_functions_and_params`, and
`ghl_initialize_tabulated_eos_functions_and_params` are declared in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h) and implemented in
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).
The simple and hybrid wrappers install dispatch first through
`ghl_initialize_eos_functions`, then call the matching parameter initializer.
The tabulated wrapper does the same only in an HDF5-capable build; its
no-HDF5 branch returns before dispatch setup.

Core Doxygen describes EOS initialization as Core setup for
`ghl_eos_parameters` and function pointers, while EOS Doxygen describes the
simple, hybrid, and tabulated EOS schemes. Keep equation, table, interpolation,
and helper implementation details routed to EOS pages.
Sources: [`docs/raw/GRHayL_Core.dox`](../../docs/raw/GRHayL_Core.dox),
[`docs/raw/EOS.dox`](../../docs/raw/EOS.dox).

## Function Pointer Dispatch

Dispatch uses process-wide mutable function-pointer globals, not fields inside
`ghl_eos_parameters`. `ghl_con2prim_multi_method` is declared in
[`GRHayL/include/ghl_con2prim.h`](../../GRHayL/include/ghl_con2prim.h) and
has storage in
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).
The EOS/flux pointer family is declared `extern` in
[`GRHayL/include/ghl_eos_functions.h`](../../GRHayL/include/ghl_eos_functions.h)
and receives storage definitions when Core includes
[`GRHayL/include/ghl_eos_functions_declaration.h`](../../GRHayL/include/ghl_eos_functions_declaration.h).
Installed public header routing is listed in
[`GRHayL/include/make.code.defn`](../../GRHayL/include/make.code.defn).

Storage gives these globals static zero initialization; it does not give them a
callable target. The file named `ghl_eos_functions_declaration.h` contains
definitions rather than `extern` declarations, despite its name. Core includes
it in one built translation unit. Installation of that header is not evidence
that consumers should include it or that every pointer is assigned.

`ghl_initialize_eos_functions` performs assignment in three stages:

1. It always calls `NRPyEOS_initialize_hybrid_functions`, assigning hybrid
   pointers and `ghl_compute_h_and_cs2`.
2. In HDF5 builds it then calls `NRPyEOS_initialize_tabulated_functions`,
   assigning tabulated pointers and temporarily overwriting
   `ghl_compute_h_and_cs2`.
3. For a recognized EOS type it selects `ghl_con2prim_multi_method` and resets
   `ghl_compute_h_and_cs2` to the requested family.

Sources:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c),
[`GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c`](../../GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c),
and
[`GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`](../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c).

For `ghl_eos_simple` and `ghl_eos_hybrid`, Core routes
`ghl_con2prim_multi_method` to `ghl_con2prim_hybrid_multi_method` and
`ghl_compute_h_and_cs2` to `NRPyEOS_hybrid_compute_enthalpy_and_cs2`.
Source:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).

For `ghl_eos_tabulated`, Core routes those pointers to tabulated Con2Prim and
tabulated enthalpy/sound-speed functions only when HDF5-backed tabulated code is
available. The tabulated branch in `ghl_initialize_eos_functions` is guarded by
`GHL_DISABLE_HDF5`.
Sources:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c),
[`wiki/build-and-ci.md`](../build-and-ci.md).

An unrecognized `ghl_eos_t` value has no `else` error path and does not return
`ghl_error_unknown_eos_type`. The EOS-selection branch is skipped, so
`ghl_con2prim_multi_method` retains its prior value (initially null). However,
family initialization has already overwritten `ghl_compute_h_and_cs2`: it
points to the hybrid implementation in a no-HDF5 build and to the tabulated
implementation in an HDF5 build. The function returns `void` with this
mode-dependent partial dispatch state
([`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c)).

Repository-wide assignment/call search finds no assignment or call for
`ghl_tabulated_free_beq_quantities` or the three generic
`ghl_calculate_HLLE_fluxes_dirn*` pointer globals. Their declarations and Core
storage make symbols available, but current repo evidence does not establish
initialized or supported dispatch through them. Direct named Flux_Source and
tabulated functions are separate APIs; route product intent for these four
globals to EOS/Flux_Source maintainers.

## Ordering, Mutation, And Failure

The direct parameter initializers are not dispatch-independent:
`ghl_initialize_simple_eos` and `ghl_initialize_hybrid_eos` call hybrid pointer
globals, while `ghl_initialize_tabulated_eos` calls tabulated pointer globals.
Call `ghl_initialize_eos_functions` first or use a matching
`*_functions_and_params` wrapper. Calling a direct initializer against only
zero-initialized pointer storage is not a checked error path.

Wrappers install global dispatch before parameter validation. If simple or
hybrid parameter initialization returns an error, the process-wide pointers
remain assigned; there is no rollback. Tabulated initialization also mutates
the EOS struct and can allocate/read table state before later atmosphere/bound
validation returns an error. Review table cleanup with the EOS lifecycle owner;
the Core wrapper is not transactional.

The tabulated convenience wrapper hard-codes
`ghl_eos_table_stellarcollapse`, `clean_sound_speed = false`, and neural-network
Con2Prim disabled before calling the configurable direct initializer. The
direct API is required for other supported option values.

## HDF5 Gate

With `GHL_DISABLE_HDF5`, `ghl_initialize_tabulated_eos` and
`ghl_initialize_tabulated_eos_functions_and_params` return
`ghl_error_used_disabled_hdf5`.
Source:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).

By contrast, calling `ghl_initialize_eos_functions(ghl_eos_tabulated)`
directly in a no-HDF5 build reaches `GHL_HDF5_ERROR_IF_USED`, which expands to
the terminating `ghl_error` macro. It exits rather than returning
`ghl_error_used_disabled_hdf5`. Keep this direct-dispatch behavior distinct
from the two error-returning parameter entry points
([`GRHayL/include/ghl_nrpyeos_tabulated.h`](../../GRHayL/include/ghl_nrpyeos_tabulated.h),
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c)).

`ghl_error_used_disabled_hdf5` is part of `ghl_error_codes_t` in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h). No-HDF5 configure behavior
defines `GHL_DISABLE_HDF5` and applies a source-selection predicate that omits
most HDF5/tabulated sources while retaining the tabulated primitive-guess and
neural-network helpers documented in
[`wiki/build-and-ci.md`](../build-and-ci.md).

## Atmosphere And Bounds Routing

At Core routing level, simple EOS initialization validates density and pressure
atmosphere inputs, gives negative density or pressure min/max inputs default
floor/ceiling behavior, and returns min-greater-than-max errors. Hybrid EOS
initialization validates density atmosphere/min/max inputs and computes
pressure/energy/entropy atmosphere and bounds through hybrid helpers. Source:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).

Those checks are not complete domain validation. Simple setup accepts zero
`rho_atm`, zero floors, and any `Gamma`; later divisions by density and
`Gamma-1` can therefore produce non-finite `eps` fields. Hybrid setup does not
validate `neos` against `1..MAX_EOS_PARAMS`, the input array extents, EOS
coefficients, or output pointer. These are unchecked caller preconditions, not
documented error-return paths.

Tabulated EOS initialization validates `rho_atm`, `Y_e_atm`, and `T_atm`,
clamps requested `rho`, `Y_e`, and `T` min/max values to table bounds, checks
min/max ordering, then computes atmosphere pressure, energy, and entropy
through table-backed interpolation. Source:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).
For atmosphere reset behavior after EOS setup, use
[`wiki/gems/atmosphere/prescription-contract.md`](../gems/atmosphere/prescription-contract.md).

`tau_atm` has family-dependent current initialization: simple and hybrid use
`rho_atm * eps_atm`, while tabulated setup uses `rho_min * eps_min`. The public
struct comment describes `_atm` fields as atmosphere values, and Con2Prim later
uses `tau_atm` as a conservative floor. Do not assume these formulas are
equivalent; intended tabulated semantics need maintainer resolution
([`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h),
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c),
[`GRHayL/Con2Prim/apply_conservative_limits.c`](../../GRHayL/Con2Prim/apply_conservative_limits.c)).

## Tests And Coverage

The Core suite tests simple EOS default density/pressure floor and ceiling
behavior, initializes hybrid EOS atmosphere data, and checks
`ghl_set_prims_to_constant_atm` for simple and hybrid EOS cases. Source:
[`Unit_Tests/unit_test_grhayl_core_test_suite.c`](../../Unit_Tests/unit_test_grhayl_core_test_suite.c).

Do not claim tabulated default or tabulated atmosphere reset coverage from the
Core suite: its tabulated setup and table-default checks are commented/TODO, and
the atmosphere loop stops before the tabulated case. Source:
[`Unit_Tests/unit_test_grhayl_core_test_suite.c`](../../Unit_Tests/unit_test_grhayl_core_test_suite.c).

Error-path coverage for invalid EOS atmosphere and min/max inputs lives in
[`Unit_Tests/unit_test_code_error.c`](../../Unit_Tests/unit_test_code_error.c),
including HDF5-only skip behavior when `GHL_DISABLE_HDF5` is set.

## Source Of Truth

- [`wiki/gems/eos/initialization-and-dispatch.md`](../gems/eos/initialization-and-dispatch.md)
- [`wiki/build-and-ci.md`](../build-and-ci.md)
- [`docs/raw/GRHayL_Core.dox`](../../docs/raw/GRHayL_Core.dox)
- [`docs/raw/EOS.dox`](../../docs/raw/EOS.dox)
- [`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h)
- [`GRHayL/include/ghl_con2prim.h`](../../GRHayL/include/ghl_con2prim.h)
- [`GRHayL/include/ghl_eos_functions.h`](../../GRHayL/include/ghl_eos_functions.h)
- [`GRHayL/include/ghl_eos_functions_declaration.h`](../../GRHayL/include/ghl_eos_functions_declaration.h)
- [`GRHayL/include/make.code.defn`](../../GRHayL/include/make.code.defn)
- [`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c)
- [`Unit_Tests/unit_test_grhayl_core_test_suite.c`](../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
- [`Unit_Tests/unit_test_code_error.c`](../../Unit_Tests/unit_test_code_error.c)
