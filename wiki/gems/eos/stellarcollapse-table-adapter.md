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

## Producer Energy Conventions

Reading a StellarCollapse-format file does not establish its producer's energy
reference or qualify a leakage blocking approximation. The actual legacy
SLy4 fixture and the regularized SRO-141 input have the following conventions;
these are not defaults for arbitrary files accepted by the loader.

| Input / component | Individual nucleon chemical potentials | Rest energies and consequence |
| --- | --- | --- |
| `SLy4_3335_rho391_temp163_ye66.h5`, SNA component | Both are referenced to the free-neutron rest energy. The proton mean field includes the negative neutron–proton rest gap; interactions remain in both potentials. | The embedded input/source selects `m_n=939.56540`, `m_p=938.27204 MeV`, giving `Q=1.29336 MeV`. Removing the proton rest offset does not remove its mean field. |
| Same legacy artifact, NSE component | The writer subtracts one neutron reference from both full potentials using its nuclear partition data. The delivered artifact blends SNA and NSE. | Neutron reference `939.5654133 MeV`; nominal proton mass `938.2720813 MeV` minus its isotope's `0.00006 MeV` binding entry gives `938.2720213 MeV`, hence `Q=1.293392 MeV`. The blend has no single exact component-independent rest gap. |
| CompOSE SRO SLy4 SNA table 141, regularized profile | Both potentials use the common reference `939.5654 MeV`; `muhat=-mu_q` retains the full difference. The electron potential includes its rest energy. | Header/profile masses are `939.5654,938.2754 MeV`, formally a `1.2900 MeV` gap. Dilute raw proton effective-mass data approaches a kinetic mass near `938.27204 MeV`; the approximately `-0.00336 MeV` residual is in the mean-field convention. Header-gap subtraction alone does not produce kinetic degeneracies. |

The legacy convention follows the actual HDF5-embedded producers:
`SNA-skyrme.in` selects non-default constants;
`SNA-src.tar.gz!src/read_input/read_skyrme_coefficients.f90` sets the masses;
`SNA-src.tar.gz!src/bulk_matter/Skyrme_bulk.f90` constructs the proton field;
`NSE-src.tar.gz!src/read_input/nse_read_nuclear_data.F90` reads nuclear masses;
and `NSE-src.tar.gz!src/make_table/write_to_table.F90` applies the common
output reference. The [SRO distribution](https://stellarcollapse.org/SROEOS.html)
distinguishes SNA, NSE and merged products. Do not infer these masses by
fitting a single `muhat(Ye=.5)` value.

For [CompOSE table 141](https://compose.obspm.fr/eos/141), the
[published raw archive](https://zenodo.org/records/14811397),
[table description](https://compose.obspm.fr/download/3D/SRO/SLy4/eos.pdf) and
[manual, sections 3.5 and 4.2.2](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf)
establish the header and microscopic conventions. In
[`compose_to_grhayl.py`](../../../tools/compose/compose_to_grhayl.py),
`raw["mu_b"]` is the generated, already rest-shifted quantity, not the full
thermodynamic baryon potential. The converter restores `profile.m_n_mev`
before subtracting `profile.m_ref_mu_mev`, then writes
`mu_p=mu_n-muhat`. Raw microscopic identifiers `10040/11040` supply
effective-mass ratios and `10050/11050` supply mean fields. The converter
does not serialize these into the runtime chemistry contract. Legacy HDF5
also contains effective-mass data not returned by the six-output callback.

The legacy internal-energy writer adds `20 MeV/baryon` before conversion to
`erg/g` and logarithmic storage; its `energy_shift` is approximately
`1.91312955e19 erg/g`. The regularized SRO conversion selects a storage
shift of `30.051798127348633 MeV/baryon`. Neither shift is a chemical-potential
reference. The reader and conversion source above do not infer a chemical
reference or a microscopic spectrum from it.

These established producer conventions do not remove the generic loader's
missing-reference boundary: EOSmaker's SFHo atomic-mass reference is a
counterexample to a universal neutron reference; see the
[CompOSE chemical mapping](../neutrinos/compose-eos-adapter-how-to.md#nucleon-chemical-potentials-and-muhat).
The separate consequences for rest corrections, mean fields, effective masses,
free fractions, and the installed reference-invariant approximation belong to
the [Neutrinos Physics And EOS Contract](../neutrinos/physics-and-eos-contract.md#chemical-potentials-and-density-derived-blocking).

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
and inversions from distinct valid initial guesses, the six-value ABI order
and range failures, no-fallback
Palenzuela recovery, analytic characteristic-speed, HLLE/entropy-flux, and
source goldens, all eight combined-leakage outputs against implementation-derived
regression goldens for the analytic and regularized SRO-141 inputs, and
cleanup. Those goldens lock the selected implementation after external rate
and thin-gas qualification; they do not independently qualify its physics.
The correction uses `rho`, `T`, `Xn`, and `Xp`, rather than raw
`mu_n` and `mu_p`, because the free fractions determine the available nucleon
populations without depending on the producer's chemical-energy zero. Its
interacting-EOS and spectral-pairing limits, and the external qualification's
reproducibility boundary, are documented in the
[nucleon-blocking conventions leaf](../neutrinos/nucleon-blocking-and-eos-conventions.md).

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

## Ground Truth References

- [SRO EOS distribution](https://stellarcollapse.org/SROEOS.html): legacy SNA, NSE and merged products with embedded producer evidence.
- [CompOSE SRO SLy4 table 141](https://compose.obspm.fr/eos/141), [raw archive](https://zenodo.org/records/14811397) and [table description](https://compose.obspm.fr/download/3D/SRO/SLy4/eos.pdf): actual input identity, header and microscopic data.
- [CompOSE manual](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf): thermodynamic, chemical-potential and microscopic quantity conventions.
