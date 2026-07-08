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
Each wrapper installs dispatch first through `ghl_initialize_eos_functions`,
then calls the matching parameter initializer.

Core Doxygen describes EOS initialization as Core setup for
`ghl_eos_parameters` and function pointers, while EOS Doxygen describes the
simple, hybrid, and tabulated EOS schemes. Keep equation, table, interpolation,
and helper implementation details routed to EOS pages.
Sources: [`docs/raw/GRHayL_Core.dox`](../../docs/raw/GRHayL_Core.dox),
[`docs/raw/EOS.dox`](../../docs/raw/EOS.dox).

## Function Pointer Dispatch

`ghl_initialize_eos_functions` writes global dispatch pointers in
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).
The assigned pointers include `ghl_con2prim_multi_method` and
`ghl_compute_h_and_cs2`. `ghl_con2prim_multi_method` is declared in
[`GRHayL/include/ghl_con2prim.h`](../../GRHayL/include/ghl_con2prim.h) and
defined in
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).
`ghl_compute_h_and_cs2` is declared in
[`GRHayL/include/ghl_eos_functions.h`](../../GRHayL/include/ghl_eos_functions.h)
and defined by
[`GRHayL/include/ghl_eos_functions_declaration.h`](../../GRHayL/include/ghl_eos_functions_declaration.h).
Installed public header routing is listed in
[`GRHayL/include/make.code.defn`](../../GRHayL/include/make.code.defn).

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

## HDF5 Gate

With `GHL_DISABLE_HDF5`, `ghl_initialize_tabulated_eos` and
`ghl_initialize_tabulated_eos_functions_and_params` return
`ghl_error_used_disabled_hdf5`.
Source:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).

`ghl_error_used_disabled_hdf5` is part of `ghl_error_codes_t` in
[`GRHayL/include/ghl.h`](../../GRHayL/include/ghl.h). No-HDF5 configure behavior
defines `GHL_DISABLE_HDF5` and omits HDF5/tabulated sources from generated build
lists; route build questions to
[`wiki/build-and-ci.md`](../build-and-ci.md).

## Atmosphere And Bounds Routing

At Core routing level, simple EOS initialization validates density and pressure
atmosphere inputs, gives negative density or pressure min/max inputs default
floor/ceiling behavior, and returns min-greater-than-max errors. Hybrid EOS
initialization validates density atmosphere/min/max inputs and computes
pressure/energy/entropy atmosphere and bounds through hybrid helpers. Source:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).

Tabulated EOS initialization validates `rho_atm`, `Y_e_atm`, and `T_atm`,
clamps requested `rho`, `Y_e`, and `T` min/max values to table bounds, checks
min/max ordering, then computes atmosphere pressure, energy, and entropy
through table-backed interpolation. Source:
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).
For atmosphere reset behavior after EOS setup, use
[`wiki/gems/atmosphere/prescription-contract.md`](../gems/atmosphere/prescription-contract.md).

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
