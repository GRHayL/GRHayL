# Tabulated Interpolator Catalog

## Purpose

This page catalogs built tabulated EOS interpolation wrappers by input tuple and
routes shared helper code. It complements
[tabulated interpolation and bounds](tabulated-interpolation-and-bounds.md) by
making wrapper ownership explicit without copying interpolation formulas.

Public function-pointer names are declared in
[`GRHayL/include/ghl_eos_functions.h`](../../../GRHayL/include/ghl_eos_functions.h).
NRPyEOS declarations, table keys, index macros, and table-bound helpers are in
[`GRHayL/include/ghl_nrpyeos_tabulated.h`](../../../GRHayL/include/ghl_nrpyeos_tabulated.h).
Assignments from public pointers to NRPyEOS implementations are in
[`GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c).

## Built Wrapper Groups

The built interpolator source list is
[`GRHayL/EOS/Tabulated/interpolators/make.code.defn`](../../../GRHayL/EOS/Tabulated/interpolators/make.code.defn).
All files listed there live under
[`GRHayL/EOS/Tabulated/interpolators/`](../../../GRHayL/EOS/Tabulated/interpolators/).

### rho_Ye_T

Wrappers whose input tuple is `(rho,Y_e,T)` call the direct shared helper
`NRPyEOS_from_rho_Ye_T_interpolate_n_quantities`. Built routes include
`P`, `eps`, `cs2`, `P+eps`, `P+eps+S`, `P+eps+cs2`,
`P+eps+S+cs2`, `P+eps+depsdT`, chemical-potential bundles, and neutron/proton
fraction bundles. Direct interpolation performs table-bound checks, chooses
interpolation spots, interpolates requested table keys, and applies source-owned
log/energy-shift output conversions for pressure, enthalpy, and internal energy.

### rho_Ye_eps

Wrappers whose input tuple is `(rho,Y_e,eps)` call
`NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities` with
`NRPyEOS_eps_key`, recover temperature using `eos->root_finding_precision`, and
then reuse the direct `(rho,Y_e,T)` helper. Built routes include `T`, `P+T`,
`P+S+T`, `P+cs2+T`, and `P+S+depsdT+T`.

### rho_Ye_P

Wrappers whose input tuple is `(rho,Y_e,P)` use the same auxiliary-temperature
helper with `NRPyEOS_press_key`. Built routes include `eps+T`,
`eps+S+T`, and `eps+cs2+T`.

### rho_Ye_S

Wrappers whose input tuple is `(rho,Y_e,S)` use the auxiliary-temperature
helper with `NRPyEOS_entropy_key`. Built routes include `P+T` and
`P+eps+T`.

### rho_Ye_h

Wrappers whose input tuple is `(rho,Y_e,h)` use the auxiliary-temperature
helper with `NRPyEOS_enthalpy_key`. Built routes include `P+eps+S+T` and
`P+eps+dPdrho+dPdeps+T`.

## Shared Interpolation Helpers

Generic interpolation, bounds checks, interpolation-spot lookup, linear
interpolation, temperature recovery, and bisection helpers are source-local to
[`GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h`](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h).
The direct helper is
[`GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_T_interpolate_n_quantities.c`](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_T_interpolate_n_quantities.c).
The auxiliary-input helper is
[`GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities.c`](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities.c).

Table index routing currently exposes `ghl_tabulated_get_index_T`, assigned to
`NRPyEOS_tabulated_get_index_T` in
[`GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c).
The target implementation is
[`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_get_index.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_get_index.c).
Header prototypes for `NRPyEOS_tabulated_get_index_rho` and
`NRPyEOS_tabulated_get_index_Ye` exist, but repo-local source currently shows
only the `NRPyEOS_tabulated_get_index_T` implementation. Treat this as
prototype-only API drift, not a new page.

Bounds enforcement helpers live in
[`GRHayL/EOS/Tabulated/NRPyEOS_enforce_table_bounds.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_enforce_table_bounds.c).
They mutate caller-provided `(rho,Y_e,T)`, `(rho,Y_e,eps)`, `(rho,Y_e,S)`, or
`(rho,Y_e,P)` variables into configured `ghl_eos_parameters` bounds. Hard-error
checks in `NRPyEOS_tabulated_helpers.h` are separate and return table-bound
error codes without clamping.

## Beta-Equilibrium And Rho Maps

Beta-equilibrium/rho-map helpers are compactly routed here instead of split into
a separate page. The source is
[`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c).

`NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T` builds cached
`Ye_of_lr` by finding the electron fraction where `munu` changes sign at a
constant temperature. `NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T`
extends that cache with `lp_of_lr`, `le_of_lr`, and `lh_of_lr`. Consumer routes
include `Ye_of_rho`, `P_from_rho`, `rho_from_P`, `eps_from_rho`,
`dP_drho_from_rho`, and `deps_dP_from_rho`. Memory cleanup is
`NRPyEOS_tabulated_free_beq_quantities`.

## Enthalpy And Sound Speed

Derived enthalpy tabulation is owned by
[`GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c).
It fills the enthalpy table slot and `table_logh` using already converted table
state and the stored energy shift.

Sound-speed normalization and optional cleaning are owned by
[`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c).
That pass converts non-relativistic table sound speeds using enthalpy, counts
non-finite, superluminal, and negative values, and only clamps values when
`eos->clean_sound_speed` is true.

Runtime `compute_h_and_cs2` routing for tabulated EOS is
[`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c).
It enforces `(rho,Y_e,T)` bounds on `ghl_primitive_quantities`, computes
pressure, internal energy, and `cs2` through the `(rho,Y_e,T)` wrapper, and
returns enthalpy to the caller.

## Source-Of-Truth Paths

- [`GRHayL/include/ghl_nrpyeos_tabulated.h`](../../../GRHayL/include/ghl_nrpyeos_tabulated.h)
- [`GRHayL/include/ghl_eos_functions.h`](../../../GRHayL/include/ghl_eos_functions.h)
- [`GRHayL/EOS/Tabulated/interpolators/`](../../../GRHayL/EOS/Tabulated/interpolators/)
- [`GRHayL/EOS/Tabulated/interpolators/make.code.defn`](../../../GRHayL/EOS/Tabulated/interpolators/make.code.defn)
- [`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c)
