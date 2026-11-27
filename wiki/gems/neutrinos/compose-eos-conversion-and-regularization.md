# CompOSE EOS Conversion And Regularization

> Integration routes from CompOSE EOS tables to GRHayL, the `sro-sly4-sna-141-regularized-v1` offline profile, table preflight checks, and raw vs. generated CompOSE decoding. · Status: confirmed
> Up: [Neutrinos](index.md)

## Summary

Supplying CompOSE equations of state to GRHayL's NRPyLeakage routines follows one of three integration routes: offline conversion to the full 19-field StellarCollapse HDF5 format, a native CompOSE-backed GRHayL EOS backend, or a direct host port. GRHayL implements the fixed offline regularized profile `sro-sly4-sna-141-regularized-v1` in `tools/compose/compose_to_grhayl.py`, which produces a physics-changing surrogate table for CompOSE table 141. Raw `eos.thermo` files and generated `eos.table`/`eoscompose.h5` outputs use different normalizations and index indirections that must never be mixed.

## Detail

### Authority Order And Integration Routes

Use three authorities in order: (1) repo-local GRHayL source and headers for callback signatures, units, and error paths, (2) the CompOSE manual and selected table data sheet for axes, fields, particle indices, and thermodynamic conventions, and (3) O'Connor-Ott `EOSmaker` source for the common-reference StellarCollapse convention.

| Route | When to select | Required implementation | Primary constraint |
| --- | --- | --- | --- |
| Offline StellarCollapse HDF5 conversion | Keep GRHayL's existing tabulated reader unchanged | Build compatible regular grids, convert units and chemical-potential zero points, derive/aggregate all required fields, write every dataset read by `NRPyEOS_stellarcollapse_read_table` | Reader requires the complete 19-field StellarCollapse dataset family, not just the six leakage outputs |
| Native CompOSE GRHayL backend | GRHayL loads CompOSE tables directly | Add table storage, loader, state mapping, interpolation, callback registration, and GRHayL errors | EOS-dependent leakage routines retain an HDF5 compile guard even if the backend avoids HDF5 |
| Direct NRPyLeakage port to another host | Host provides its own GRMHD, EOS, interpolation, and errors | Implement the six-output callback and translate GRHayL types, errors, units, and HDF5 guard | Reproducing current GRHayL behavior requires keeping documented conversion conventions explicit |

CompOSE's `eoscompose.h5` is incompatible with `NRPyEOS_stellarcollapse_read_table`. The StellarCollapse format requires 19 three-dimensional datasets, three 1D grid arrays, three grid-size scalars, and `energy_shift` — including pressure, entropy, sound speed, internal energy and shift, thermodynamic derivatives, adiabatic-index convention, `munu`, and nuclear composition (`Abar`, `Zbar`, `Xh`, `Xa`). If those non-leakage fields cannot be constructed consistently, build a native backend instead of writing placeholders; grouped nuclear fields never replace free `X_n` and `X_p`.

Evidence (descriptive behavior): `GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c::NRPyEOS_stellarcollapse_read_table`; corroboration: `tools/compose/compose_to_grhayl.py::main`.

### Implemented Offline Profile `sro-sly4-sna-141-regularized-v1`

`tools/compose/compose_to_grhayl.py` (documented in `tools/compose/README.md`) implements the fixed profile `sro-sly4-sna-141-regularized-v1`, reading the official generated-HDF5 output for CompOSE table 141 and emitting GRHayL's StellarCollapse HDF5 format without adding a native C backend, resampling engine, or new public C API.

This profile intentionally produces a physics-changing regularized surrogate: strict temperature-ray regression, entropy reconstruction, joint causal derivative projection, chemical canonicalization, and constrained composition projection repair invalid raw fields. Never treat its output as official unmodified table 141 or as an oracle for raw CompOSE values; inspect the embedded manifest and tool report for exact policies and measured deltas. Independent trilinear interpolation of `Abar`, `Zbar`, and `Xh` does not preserve nonlinear charge closure away from grid nodes (recorded diagnostic; NRPyLeakage consumes `X_n` and `X_p`, not `Abar`/`Zbar`/`Xh`). `Unit_Tests/compose/` and `Unit_Tests/unit_test_tabulated_eos_compose.c` verify serialization, runtime interpolation, temperature inversions from distinct guesses, six-output callback order, and range errors with runtime sound-speed cleaning disabled.

Evidence (descriptive behavior): `tools/compose/README.md::sro-sly4-sna-141-regularized-v1`; corroboration: `Unit_Tests/unit_test_tabulated_eos_compose.c::main`.

### Prerequisites And Table Preflight Checklist

Every route requires the table's web page, data sheet, parameter/thermodynamic files, and individual free-nucleon composition data. Generating or validating against official `compose` output requires a matching `compose` release and its documented HDF5 toolchain (build defaults can differ between the manual and current source). Current EOS-dependent GRHayL leakage calls require an HDF5-enabled build. Preflight every candidate CompOSE table before coding:

| Check | Required evidence | Stop condition |
| --- | --- | --- |
| Parameter space | Table page, data sheet, `eos.t`, `eos.nb`, `eos.yq` cover intended positive `T`, `n_b`, `Y_q` | Cold, beta-equilibrated, 1D, or truncated table cannot represent evolved states |
| Thermodynamic metadata | First row of raw `eos.thermo` gives model `m_n`, `m_p`, and lepton flag `I_l` | Metadata missing or detached from generated files with no provenance |
| Lepton model | `I_l`, particle list, and data sheet specify electrons, muons, and any $\mu_{l_e}/\mu_{l_\mu}$ relation | Electron chemical potential cannot be reconstructed unambiguously |
| Composition | `eos.compo` supplies individual net neutron/proton indices `10` and `11`; negligible antinucleons | Only grouped nuclei, total nucleons, absent free nucleons, or appreciable antinucleons |
| Charge coordinate | Hydrodynamic `Y_e` maps deterministically to CompOSE `Y_q` under table neutrality rules | Evolved `Y_e` does not determine `Y_q` |
| Source layer | Reader knows whether input is raw `Q3`-`Q5` or generated physical MeV outputs | Double scaling or lost neutron-mass shift |
| Grid and interpolation | Exact axis arrays, ordering, bounds, and interpolation order are known | Implicit extrapolation or assumed uniform spacing not in source table |

### Acquiring Data And Decoding Raw Vs. Generated CompOSE Formats

CompOSE provides raw ASCII tables (`eos.t`, `eos.nb`, `eos.yq`, required `eos.thermo`, optional `eos.compo`/`eos.micro`) and official `compose` software producing `eos.quantities`, `eos.parameters`, and `eos.table` or `eoscompose.h5` (no public REST API exists; script against downloaded files or official code). In `compose` task 1, explicitly select thermodynamic identifiers `3`, `4`, `5` and composition-pair identifiers `10`, `11` (plus every extra field required for offline StellarCollapse conversion); unselected identifiers are omitted from generated outputs.

Decode **exactly one** input layer:

- **Raw `eos.thermo`:** Header row 1 provides model $m_n$, $m_p$ (MeV), and $I_l$. Raw columns store $Q_3 = \mu_b/m_n - 1$, $Q_4 = \mu_q/m_n$, and $Q_5 = \mu_l/m_n$; recover relativistic rest-mass-inclusive conserved-charge potentials via $\mu_b = m_n(1 + Q_3)$, $\mu_q = m_n Q_4$, $\mu_l = m_n Q_5$. Similarly, raw $Q_7 = e/(n_b m_n) - 1$ gives total energy per baryon $E = m_n(1 + Q_7)$ including rest mass; never relabel raw $Q_7$ as hydrodynamic specific internal energy `eps` because $m_n n_b$ is not total mass density.
- **Generated `eos.table` or `eoscompose.h5`:** Parse ASCII `eos.table` columns (`T`, `n_b`, `Y_q`, then selected quantities) from `eos.quantities`, never hard-coded column positions. In `eoscompose.h5` (per official `hdf5compose.f90`, whose group hierarchy is omitted from manual tables), axes reside at `Parameters/nb`, `Parameters/t`, `Parameters/yq` with attributes `pointsnb`, `pointst`, `pointsyq`, and `tabulation_scheme`.
- **Generated index lookup:** Identifier `3` stores `bshift` $= \mu_b - m_n$ (MeV), `4` stores $\mu_q$ (MeV), and `5` stores $\mu_l$ (MeV) — never multiply them by $m_n$. Locate thermodynamic slots dynamically via `Thermo_qty/index_thermo` into `Thermo_qty/thermo` (attribute `pointsqty`, group omitted if zero selected), recover $\mu_b = \texttt{bshift} + m_n$, and locate free neutrons/protons (`10`, `11`) via `Composition_pairs/index_yi` into `Composition_pairs/yi` (attribute `pointspairs`). Abort if required groups or indices are absent, and preserve raw-header $m_n$, $m_p$, and $I_l$ alongside generated HDF5 because the manual's generated dataset list does not guarantee all three.

Evidence (descriptive behavior): `tools/compose/compose_to_grhayl.py::index_thermo`; corroboration: `tools/compose/compose_to_grhayl.py::index_yi`.

## Sources

- [tools/compose/compose_to_grhayl.py](../../../tools/compose/compose_to_grhayl.py) - `main`, `index_thermo`, `index_yi`, fixed SRO-141 regularized converter
- [tools/compose/README.md](../../../tools/compose/README.md) - `sro-sly4-sna-141-regularized-v1` policy and invocation
- [GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c) - `NRPyEOS_stellarcollapse_read_table` 19-field reader
- [Unit_Tests/unit_test_tabulated_eos_compose.c](../../../Unit_Tests/unit_test_tabulated_eos_compose.c) - runtime integration checks for converted table
- [CompOSE Reference Manual Version 3.01](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf) - sections 3.1, 3.3-3.5, 4.1, 4.2, 7.3, 7.5, tables 3.2, 3.3, 7.1, 7.7, 7.8, appendices A-C
- [CompOSE Quick Guide For Users](https://compose.obspm.fr/download/pdf/CompOSE_Quick_Guide_for_Users.pdf) - task-1 workflow
- [CompOSE Table Catalog](https://compose.obspm.fr/table/) - table pages and data sheets
- [CompOSE Software Page](https://compose.obspm.fr/software/) - official `compose` releases
- [Official CompOSE Code Repository](https://gitlab.obspm.fr/data_and_software_compose/code-compose/-/tree/master) - `compose` source
- [Official CompOSE HDF5 Writer (hdf5compose.f90)](https://gitlab.obspm.fr/data_and_software_compose/code-compose/-/raw/master/hdf5compose.f90) - `Parameters`, `Thermo_qty`, `Composition_pairs` HDF5 hierarchy

## See Also

- Depends on: [CompOSE EOS Leakage Callback](compose-eos-leakage-callback.md)
- Depends on: [StellarCollapse Table Adapter](../eos/stellarcollapse-table-adapter.md)
- See also: [EOS Inputs And Chemical Potentials](eos-inputs-and-chemical-potentials.md)
- Validated by: [EOS Tests And Reference Data](../eos/tests-and-reference-data.md)
