# StellarCollapse Table Adapter

> HDF5 dataset reading and conversion for `ghl_eos_table_stellarcollapse`, producer chemical-potential energy references, and the offline CompOSE regularized-table interface. · Status: confirmed
> Up: [Equation Of State (EOS)](index.md)

## Summary

`GRHayL/EOS/Tabulated/stellarcollapse/` implements the sole native tabulated EOS reader (`ghl_eos_table_stellarcollapse`), reading HDF5 datasets via `NRPyEOS_hdf5_helpers.*` into a temporary adapter struct and mapping them into `ghl_eos_parameters`. Reading a StellarCollapse HDF5 file does **not** establish the producer's nucleon chemical-potential energy reference (`SLy4` legacy SNA/NSE, CompOSE SRO-141, and EOSmaker SFHo differ). Offline `tools/compose/compose_to_grhayl.py` emits this same HDF5 format without changing C runtime dispatch.

## Detail

### HDF5 Helper Layer And Temporary Adapter Conversion

- **Generic HDF5 helpers (`NRPyEOS_hdf5_helpers.c`, `.h`):** open named datasets on an open HDF5 handle, verify exact element counts, allocate buffers, read integer or double arrays, and close handles without knowing StellarCollapse semantics.
- **Reader & cleanup (`stellarcollapse/NRPyEOS_stellarcollapse.c`, `.h`):** defines the dataset-name list, initializes output to `NULL`, frees partial reads on error, reads scalar dims, `energy_shift`, optional 1-element `have_rel_cs2`, uniform increasing grid axes, and table quantities into a temporary struct freed by `NRPyEOS_stellarcollapse_free_table`.
- **Conversion (`stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c`):** rechecks dimension/count representability, allocates the six GRHayL arrays (freeing and nulling all six on allocation failure), converts density-grid and `energy_shift` units, and maps StellarCollapse slots into `NRPyEOS_keys`. Remaining unit conversions, enthalpy tabulation, and sound-speed adjustment run in `NRPyEOS_read_table_set_EOS_params.c`. Excluded under `configure --disable-hdf5`.

Evidence (descriptive behavior): `GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c::NRPyEOS_stellarcollapse_read_table`; corroboration: `GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c` and `GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c`.

### Producer Chemical-Potential Energy References

| Table / component | Nucleon chemical-potential reference | Rest-energy details & consequence |
| --- | --- | --- |
| Legacy `SLy4_3335_rho391_temp163_ye66.h5` (SNA) | Both `mu_n`, `mu_p` referenced to free-neutron rest energy; proton mean field includes negative `n - p` rest gap plus interactions. | Embedded `SNA-skyrme.in` / `read_skyrme_coefficients.f90` / `Skyrme_bulk.f90` set `m_n = 939.56540`, `m_p = 938.27204 MeV` (`Q = 1.29336 MeV`). |
| Legacy `SLy4_3335_rho391_temp163_ye66.h5` (NSE, blended with SNA) | `write_to_table.F90` subtracts one neutron reference (`939.5654133 MeV`) from both full potentials. | Proton mass `938.2720813 - 0.00006 = 938.2720213 MeV` (`Q = 1.293392 MeV`); blended file has no single exact component-independent gap. |
| CompOSE SRO SLy4 table 141 (`sro-sly4-sna-141-regularized-v1`) | Both use common reference `939.5654 MeV`; `muhat = -mu_q` keeps full difference; `mu_e` includes `m_e`. | Header masses `939.5654, 938.2754 MeV` (`1.2900 MeV` gap), while dilute effective mass approaches `938.27204 MeV` (`-0.00336 MeV` residual in mean field). |

In `tools/compose/compose_to_grhayl.py`, `raw["mu_b"]` is already rest-shifted by CompOSE (`manual_v3.00.pdf` Secs. 3.5, 4.2.2): the converter adds `profile.m_n_mev`, subtracts `profile.m_ref_mu_mev`, and writes `mu_p = mu_n - muhat`. Neither raw CompOSE microscopic IDs (`10040/11040` effective-mass ratios, `10050/11050` mean fields) nor legacy HDF5 effective-mass datasets are exposed by NRPyLeakage's six-output callback. Storage energy shifts (`20 MeV/baryon`, `~1.91312955e19 erg/g` in legacy SLy4; `30.051798127348633 MeV/baryon` in regularized SRO-141) are unrelated to chemical-potential references, and EOSmaker's SFHo atomic-mass reference confirms no universal neutron reference exists. For leakage blocking consequences, see [EOS Inputs And Chemical Potentials](../neutrinos/eos-inputs-and-chemical-potentials.md) and [Nucleon Blocking And EOS Conventions](../neutrinos/nucleon-blocking-and-eos-conventions.md).

### Offline CompOSE Regularized-Table Interface

`tools/compose/compose_to_grhayl.py` produces the fixed surrogate profile `sro-sly4-sna-141-regularized-v1` in StellarCollapse HDF5 layout (`have_rel_cs2 = 1`, requiring `clean_sound_speed = false`) with no C runtime dispatch changes; GRHayL's reader loads the standard axes, scalars, and 19 fields while ignoring `grhayl_compose/manifest_json`. Mass/charge closure holds at grid nodes, whereas off-grid trilinear interpolation of `Abar`, `Zbar`, and `Xh` cannot preserve the nonlinear heavy-charge product; NRPyLeakage instead uses the six-output `(muhat, mu_e, mu_p, mu_n, Xn, Xp)` callback and density-derived blocking. See [CompOSE EOS Conversion And Regularization](../neutrinos/compose-eos-conversion-and-regularization.md) and [CompOSE EOS Leakage Callback](../neutrinos/compose-eos-leakage-callback.md).

Evidence (descriptive behavior): `tools/compose/compose_to_grhayl.py`; corroboration: `Unit_Tests/unit_test_tabulated_eos_compose.c`.

## Sources

- [GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_hdf5_helpers.c) - generic HDF5 dataset readers
- [GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c) - StellarCollapse reader and dataset names
- [GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c](../../../GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c) - conversion into `ghl_eos_parameters`
- [GRHayL/EOS/Tabulated/stellarcollapse/make.code.defn](../../../GRHayL/EOS/Tabulated/stellarcollapse/make.code.defn) - adapter manifest
- [tools/compose/compose_to_grhayl.py](../../../tools/compose/compose_to_grhayl.py) - offline SRO-141 regularized converter
- [tools/compose/README.md](../../../tools/compose/README.md) - converter controls and distortion disclosures
- [Unit_Tests/unit_test_tabulated_eos_compose.c](../../../Unit_Tests/unit_test_tabulated_eos_compose.c) - C integration test for converted tables
- [SRO EOS distribution](https://stellarcollapse.org/SROEOS.html) - legacy SNA, NSE, and merged products
- [CompOSE SRO SLy4 table 141](https://compose.obspm.fr/eos/141) - CompOSE table 141 catalog entry
- [CompOSE table 141 raw archive](https://zenodo.org/records/14811397) - published raw table archive
- [CompOSE SRO SLy4 datasheet](https://compose.obspm.fr/download/3D/SRO/SLy4/eos.pdf) - table 141 masses and bounds
- [CompOSE Manual v3.00](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf) - Secs. 3.5 and 4.2.2 chemical-potential conventions

## See Also

- Depends on: [Tabulated Table Layout](tabulated-table-layout.md)
- See also: [CompOSE EOS Conversion And Regularization](../neutrinos/compose-eos-conversion-and-regularization.md)
- See also: [EOS Inputs And Chemical Potentials](../neutrinos/eos-inputs-and-chemical-potentials.md)
- See also: [Nucleon Blocking And EOS Conventions](../neutrinos/nucleon-blocking-and-eos-conventions.md)
