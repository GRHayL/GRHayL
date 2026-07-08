# Tabulated EOS Table Contract

## Purpose

This page routes agents through the tabulated EOS table lifecycle: HDF5 table
loading, stellar-collapse adapter conversion, GRHayL table state, energy shift,
enthalpy tabulation, sound-speed handling, and memory cleanup. It is a routing
contract only; exact enums, datasets, macros, and signatures stay in source.

## Current Adapter

The current table adapter path is `ghl_eos_table_stellarcollapse`, and it is
HDF5-backed. The public wrapper
`ghl_initialize_tabulated_eos_functions_and_params` currently hard-codes
`ghl_eos_table_stellarcollapse` and `clean_sound_speed=false` before calling
`ghl_initialize_tabulated_eos`; see
[`GRHayL/GRHayL_Core/initialize_eos.c`](../../../GRHayL/GRHayL_Core/initialize_eos.c).
The table type enum lives in [`GRHayL/include/ghl.h`](../../../GRHayL/include/ghl.h).

## Table Keys And Quantities

`NRPyEOS_keys` in
[`GRHayL/include/ghl_nrpyeos_tabulated.h`](../../../GRHayL/include/ghl_nrpyeos_tabulated.h)
is the exact table-key authority. At category level, the GRHayL table stores
thermodynamic values, sound speed, derivatives, chemical potentials,
composition fractions, nuclear averages, adiabatic index, and the derived
enthalpy slot. Do not copy the full enum into KB pages.

The stellar-collapse adapter has its own quantity enum and HDF5 dataset-name
mapping in
[`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h)
and
[`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c).
The conversion map from stellar-collapse quantities into `NRPyEOS_keys` belongs
to
[`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c).

## Read Flow

1. Public init: `ghl_initialize_tabulated_eos_functions_and_params` initializes
   tabulated function pointers, selects the current default table adapter, and
   passes `clean_sound_speed=false` to the low-level initializer in
   [`GRHayL/GRHayL_Core/initialize_eos.c`](../../../GRHayL/GRHayL_Core/initialize_eos.c).
2. Low-level init: `ghl_initialize_tabulated_eos` records EOS type, table type,
   and `clean_sound_speed`, then calls the table read function and derives EOS
   atmosphere, floor, ceiling, and root-finding state.
3. Table read: `NRPyEOS_read_table_set_EOS_params` validates `ghl_eos_tabulated`,
   dispatches on `eos->table_type`, and currently routes
   `ghl_eos_table_stellarcollapse` through the stellar-collapse reader in
   [`GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c).
4. HDF5 adapter: `NRPyEOS_stellarcollapse_read_table` opens the HDF5 file,
   reads scalar dimensions, `energy_shift`, grid arrays, table-backed
   quantities, and `have_rel_cs2` when present; HDF5 dataset helpers live in
   [`GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c)
   and
   [`GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.h`](../../../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.h).
5. Conversion: `NRPyEOS_stellarcollapse_to_ghl` allocates `ghl_eos_parameters`
   table arrays, converts grid units, copies mapped table quantities, and stores
   the converted energy shift.
6. Post-read processing: `NRPyEOS_read_table_set_EOS_params` converts pressure,
   energy, sound-speed, and derivative units; fills `table_eps`; calls
   `NRPyEOS_tabulate_enthalpy`; calls
   `NRPyEOS_tabulated_adjust_sound_speed`; computes interpolation stride
   inverses and table bounds.
7. Cleanup: the temporary stellar-collapse table is freed after conversion, and
   GRHayL table memory is later released through `NRPyEOS_free_memory`; see
   [`GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c).

## Units And Energy Shift

Unit conversion constants are owned by
[`GRHayL/include/ghl_nrpyeos_tabulated.h`](../../../GRHayL/include/ghl_nrpyeos_tabulated.h).
Use that header for exact `CODE_TO_CGS_*` and `CGS_TO_CODE_*` values rather
than duplicating macro blocks. The stellar-collapse conversion applies density,
temperature-grid, and energy shift conversions before storing data in
`ghl_eos_parameters`; the table read path then converts pressure, specific
energy, sound speed, and derivative quantities to code units.

The stellar-collapse table stores an energy shift. GRHayL converts that shift,
stores it as `eos->energy_shift`, tabulates `table_eps` from shifted internal
energy data, and subtracts the shift when computing table energy bounds. The
enthalpy tabulation uses the same shift when forming the derived enthalpy table;
see
[`GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c).

## Sound Speed

The HDF5 adapter checks `have_rel_cs2`. If the flag is absent, the table treats
sound speed as non-relativistic; otherwise the dataset value decides. The
adjustment pass converts non-relativistic sound speed through enthalpy and then
counts non-finite, superluminal, and negative values. When
`clean_sound_speed` is true, those values are clamped; with the current public
wrapper default `clean_sound_speed=false`, the pass reports them without fixing
the table values. The implementation is
[`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c).

## HDF5 Build Gate

Tabulated EOS runtime table support depends on HDF5. The default `configure`
path builds with HDF5, while `./configure --disable-hdf5` defines
`GHL_DISABLE_HDF5` and omits tabulated/HDF5 implementation sources from the
generated build. The EOS source list for HDF5-enabled tabulated builds is
[`GRHayL/EOS/Tabulated/make.code.defn`](../../../GRHayL/EOS/Tabulated/make.code.defn);
build behavior is routed through [`configure`](../../../configure) and
[`README.md`](../../../README.md). Under `GHL_DISABLE_HDF5`, HDF5-backed entry
points return or raise disabled-HDF5 errors instead of loading tables.

## Tests And Fixtures

Direct tabulated EOS read/interpolation coverage is in
[`Unit_Tests/unit_test_tabulated_eos.c`](../../../Unit_Tests/unit_test_tabulated_eos.c).
The repo-local sample table fixture and generator live under
[`Unit_Tests/sample_table/`](../../../Unit_Tests/sample_table/), including
`sample_table/table_info.txt` for fixture dimensions. CI routes include
`download_test_data EOS/simple_table.h5`, `./test/unit_test_tabulated_eos
simple_table.h5`, and larger downloaded stellar-collapse tables in
[`.github/run_tests.sh`](../../../.github/run_tests.sh).

## Source-Of-Truth Paths

- [`GRHayL/include/ghl_nrpyeos_tabulated.h`](../../../GRHayL/include/ghl_nrpyeos_tabulated.h)
- [`GRHayL/include/ghl.h`](../../../GRHayL/include/ghl.h)
- [`GRHayL/GRHayL_Core/initialize_eos.c`](../../../GRHayL/GRHayL_Core/initialize_eos.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.h`](../../../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.h)
- [`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c)
- [`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c)
- [`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h)
- [`GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_adjust_sound_speed.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_free_memory.c)
- [`GRHayL/EOS/Tabulated/make.code.defn`](../../../GRHayL/EOS/Tabulated/make.code.defn)
- [`README.md`](../../../README.md)
- [`configure`](../../../configure)
- [`.github/run_tests.sh`](../../../.github/run_tests.sh)
- [`docs/raw/EOS.dox`](../../../docs/raw/EOS.dox)
- [`Unit_Tests/unit_test_tabulated_eos.c`](../../../Unit_Tests/unit_test_tabulated_eos.c)
- [`Unit_Tests/sample_table/`](../../../Unit_Tests/sample_table/)
