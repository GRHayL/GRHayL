# Stellarcollapse Table Adapter

## Purpose

This page routes the HDF5-backed stellar-collapse table adapter for tabulated
EOS. It focuses on `table_type` dispatch, HDF5 dataset reads, temporary adapter
ownership, unit conversion, and conversion into `ghl_eos_parameters`. For
general tabulated lifecycle context, see
[tabulated table contract](tabulated-table-contract.md).

## Dispatch Boundary

`ghl_eos_table_stellarcollapse` is the only concrete table adapter currently
handled by
[`GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c).
That entry point validates that the caller supplied a tabulated EOS, switches
on `eos->table_type`, and routes the stellar-collapse case through a local
reader before applying shared post-read processing. Table type values and
`ghl_eos_parameters` fields live in
[`GRHayL/include/ghl.h`](../../../GRHayL/include/ghl.h); tabulated keys,
indexing macros, unit macros, and public NRPyEOS declarations live in
[`GRHayL/include/ghl_nrpyeos_tabulated.h`](../../../GRHayL/include/ghl_nrpyeos_tabulated.h).

## HDF5 Helper Boundary

Generic HDF5 reads belong to
[`GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c)
and
[`GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.h`](../../../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.h).
Those helpers open named datasets from an already-open file, verify the total
dataset size against the expected table size, allocate the destination array,
read either integer or double data, and close HDF5 handles on the way out. They
do not know stellar-collapse quantity semantics; callers own dataset names,
expected sizes, and mapping into EOS fields.

The stellar-collapse reader in
[`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c)
owns those semantics. It opens the HDF5 file, reads scalar grid dimensions and
`energy_shift`, reads grid arrays, reads the stellar-collapse quantity datasets,
and checks `have_rel_cs2` when present. The exact dataset-name list is source
owned by that file, while the adapter struct and stellar-collapse quantity enum
are declared in
[`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h).

## Ownership And Conversion Flow

The stellar-collapse table is a temporary adapter object. Its allocation,
dataset pointers, and cleanup are owned by
`NRPyEOS_stellarcollapse_read_table` and
`NRPyEOS_stellarcollapse_free_table`; GRHayL-owned table storage begins only
after
[`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c)
allocates `ghl_eos_parameters` arrays.

The reader uses zero-initialized temporary storage and frees all successfully
read datasets on a read error. Conversion separately allocates six GRHayL
arrays. On conversion allocation failure it frees those allocations and
returns `ghl_error_out_of_memory`, but does not reset corresponding EOS fields
to `NULL`; this is partial-state evidence, not a safe retry/cleanup contract.

Conversion responsibilities split this way:

- `NRPyEOS_stellarcollapse_to_ghl` copies table dimensions into
  `ghl_eos_parameters`, allocates `table_logrho`, `table_logT`, `table_Y_e`,
  `table_all`, `table_eps`, and `table_logh`, converts density-grid and
  energy-shift units, and maps stellar-collapse quantity slots into
  `NRPyEOS_keys`.
- `NRPyEOS_read_table_set_EOS_params` applies the remaining post-read code-unit
  conversions for pressure, internal energy, sound speed, and derivative slots;
  fills `table_eps`; builds derived enthalpy; adjusts sound speed; computes
  interpolation stride inverses; and stores table min/max bounds on
  `ghl_eos_parameters`.

Do not duplicate unit constants or dataset maps in KB pages. Use
[`GRHayL/include/ghl_nrpyeos_tabulated.h`](../../../GRHayL/include/ghl_nrpyeos_tabulated.h)
and the stellar-collapse source files for exact names and values.

## Offline CompOSE Regularization

[`tools/compose/compose_to_grhayl.py`](../../../tools/compose/compose_to_grhayl.py)
produces this same StellarCollapse schema from one fixed CompOSE-generated
SRO(SLy4) SNA table-141 schema. The production identity is
`sro-sly4-sna-141-regularized-v1`. It is a physics-changing surrogate, not an
unmodified CompOSE artifact or a native table backend. Exact input controls,
regularization equations, distortion disclosure, and commands belong to the
[`tools/compose` README](../../../tools/compose/README.md).

No runtime dispatch or public C surface changes for this route. The output sets
`have_rel_cs2=1`, so the loader retains its relativistic sound speed; callers
must leave optional runtime sound-speed cleaning disabled. The reader consumes
the standard axes, scalars, and 19 fields and ignores the converter's diagnostic
`grhayl_compose/manifest_json` group.

Composition mass and charge closure are hard node gates. Independent
trilinear interpolation of `Abar`, `Zbar`, and `Xh` cannot preserve the
nonlinear heavy-charge product off grid. Current NRPyLeakage calls instead use
the six-value `muhat,mu_e,mu_p,mu_n,Xn,Xp` wrapper; do not generalize this into
an off-grid equilibrium-composition claim. The auto-discovered
[`unit_test_tabulated_eos_compose.c`](../../../Unit_Tests/unit_test_tabulated_eos_compose.c)
checks every serialized node and all 19 mappings, storage-space interpolation
and inversions, the six-value ABI order and range failures, no-fallback
Palenzuela recovery, analytic characteristic-speed, HLLE/entropy-flux, and
source goldens, all eight combined-leakage outputs against an independent
high-precision Ruffert-equation oracle, and cleanup.

## Build Gate

The adapter is an HDF5 source/build boundary. With `GHL_DISABLE_HDF5`, the
tabulated table read path returns a disabled-HDF5 error instead of opening
files, and HDF5-only code is excluded by preprocessor guards. The configured
no-HDF5 build path is controlled by `--disable-hdf5` in
[`configure`](../../../configure), which defines `GHL_DISABLE_HDF5` and omits
this adapter and other HDF5-dependent EOS sources from generated builds.
Con2Prim retains specific non-HDF5 tabulated helper/NN exceptions; none makes
the stellar-collapse adapter available. This is a source and build contract,
not a Doxygen task.

Built stellar-collapse adapter sources are listed in
[`GRHayL/EOS/Tabulated/stellarcollapse/make.code.defn`](../../../GRHayL/EOS/Tabulated/stellarcollapse/make.code.defn).
The parent tabulated source list is
[`GRHayL/EOS/Tabulated/make.code.defn`](../../../GRHayL/EOS/Tabulated/make.code.defn).

## Evidence Links

Sample tables and tests are evidence only. Use
[`Unit_Tests/unit_test_tabulated_eos.c`](../../../Unit_Tests/unit_test_tabulated_eos.c),
[`Unit_Tests/unit_test_tabulated_eos_compose.c`](../../../Unit_Tests/unit_test_tabulated_eos_compose.c),
[`Unit_Tests/compose/`](../../../Unit_Tests/compose/),
[`Unit_Tests/sample_table/`](../../../Unit_Tests/sample_table/), and
[`.github/run_tests.sh`](../../../.github/run_tests.sh) to verify how repo-local
fixtures and CI exercise table loading. Do not treat fixture contents as adapter
API authority.

## Source-Of-Truth Paths

- [`GRHayL/include/ghl_nrpyeos_tabulated.h`](../../../GRHayL/include/ghl_nrpyeos_tabulated.h)
- [`GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_read_table_set_EOS_params.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c`](../../../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c)
- [`GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.h`](../../../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.h)
- [`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c)
- [`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h)
- [`GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c`](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c)
- [`GRHayL/EOS/Tabulated/stellarcollapse/make.code.defn`](../../../GRHayL/EOS/Tabulated/stellarcollapse/make.code.defn)
- [`tools/compose/compose_to_grhayl.py`](../../../tools/compose/compose_to_grhayl.py)
- [`tools/compose/README.md`](../../../tools/compose/README.md)
