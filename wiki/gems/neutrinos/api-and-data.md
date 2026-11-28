# Neutrinos API And Data

> Public radiation containers, `NRPyLeakage_*` entry points, input/unit preconditions, error codes, and HDF5/EOS boundaries. · Status: confirmed
> Up: [Neutrinos](index.md)

## Summary

`GRHayL/include/ghl_radiation.h` defines the three-species radiation containers (`ghl_neutrino_luminosities`, `ghl_neutrino_opacities`, `ghl_neutrino_optical_depths`) and includes `GRHayL/include/ghl_nrpyleakage.h`, which exports the `NRPyLeakage_*` entry points, constants, unit macros, and finite-value classifiers. No `ghl_*` wrapper family exists for leakage calls. EOS-dependent routines require an HDF5-enabled build and an initialized six-output tabulated-EOS chemistry callback, whereas the Fermi-Dirac helper and optical-depth update do not call the EOS or HDF5.

## Detail

### Radiation Data Containers

`GRHayL/include/ghl_radiation.h` includes `ghl.h`, declares the radiation containers, and then includes `ghl_nrpyleakage.h`:

- `ghl_neutrino_luminosities` holds scalar `nue`, `anue`, and `nux` fields in that declaration order.
- `ghl_neutrino_opacities` holds two-element arrays `nue[2]`, `anue[2]`, and `nux[2]`.
- `ghl_neutrino_optical_depths` is a `typedef` alias of `ghl_neutrino_opacities`.

Public headers do not assign named labels to slots `[0]` and `[1]`. Implementation formulas pair slot `[0]` with `R_source` and equilibrium-degeneracy interpolation and pair slot `[1]` with `Q_source` and luminosity suppression; keep `[0]`/`[1]` indexing at caller boundaries. Binary luminosity reference files serialize both `ghl_neutrino_optical_depths` (unperturbed input records) and `ghl_neutrino_luminosities` as raw structs, making field declaration order in both `ghl_neutrino_opacities` and `ghl_neutrino_luminosities` part of luminosity reference layout, whereas constant-density-sphere opacity and optical-depth reference files store twelve separate arrays in species-then-slot order.

Evidence (scientific/interface requirement): `ghl_radiation.h::ghl_neutrino_opacities`; corroboration: `ghl_radiation.h::ghl_neutrino_luminosities`.

### Public NRPyLeakage Entry Points And Portability Boundary

`GRHayL/include/ghl_nrpyleakage.h` exports:

- `NRPyLeakage_Fermi_Dirac_integrals` (returns `ghl_error_codes_t`)
- `NRPyLeakage_compute_neutrino_opacities` (returns `ghl_error_codes_t`)
- `NRPyLeakage_compute_neutrino_luminosities` (returns `ghl_error_codes_t`)
- `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms` (returns `ghl_error_codes_t`)
- `NRPyLeakage_optical_depths_PathOfLeastResistance` (returns `void`)

Every implementation file includes `ghl_radiation.h`, binding signatures and bodies to `ghl_eos_parameters`, the radiation structs, `ghl_error_codes_t` (`ghl_error_used_disabled_hdf5`, `ghl_error_invalid_fermi_dirac_integral_key`, `ghl_error_nrpyleakage_blocking`, `ghl_error_nrpyleakage_nonfinite_output`), `ghl_nrpyleakage.h` unit/constant/finiteness macros (`NRPYLEAKAGE_FD_OR_RETURN`), and `-lm`. The opacity, combined source/opacity, and luminosity files also include `ghl_nrpyeos_tabulated.h`, which pulls in `hdf5.h` when HDF5 is enabled even though leakage reads EOS quantities solely through the chemistry callback.

Porting NRPyLeakage to another host requires: (1) supplying the six EOS outputs in the units specified in [EOS Inputs And Chemical Potentials](eos-inputs-and-chemical-potentials.md), (2) preserving or translating the three-species/two-slot arrays, (3) replacing GRHayL return codes and `NRPYLEAKAGE_FD_OR_RETURN` without altering early-exit writeback behavior, (4) keeping density, opacity, source, luminosity, and cgs suppression-ratio conversions mutually consistent, (5) carrying the private blocking and rate headers (or equivalent behavior), and (6) revalidating finite-value classification, metric conventions, and path lengths.

Evidence (descriptive behavior): `ghl_nrpyleakage.h::NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`; corroboration: `ghl_nrpyleakage.h::NRPYLEAKAGE_FD_OR_RETURN`.

### Input Preconditions, Constants, And Floating-Point Discipline

Routines do not perform general null-pointer checks. Caller preconditions established by the declarations and formulas are:

- `rho` / `rho_b` is in GRHayL geometric density units and is converted internally via `NRPyLeakage_units_geom_to_cgs_D`.
- `Y_e` and positive nonzero `T` (in MeV) must lie within the initialized tabulated EOS domain; formulas divide directly by `T`.
- Opacity outputs are multiplied by `NRPyLeakage_units_geom_to_cgs_L` (inverse geometric length) to pair with geometric coordinate path lengths in dimensionless optical depths.
- `R_source` and `Q_source` are signed fluid-frame scalar rates per proper volume and proper time, already carrying the baryon mass factor in `R_source` and the factor of four for heavy leptons in `Q_source`; metric, lapse, Lorentz, and velocity projections remain caller responsibilities (see [Physics Model And Sources](physics-model-and-sources.md)).
- Luminosities carry `W * alpha^2 * sqrt(det(gamma_ij))` and require finite positive lapse, finite `W >= 1`, and a finite positive-definite spatial metric (a nonnegative determinant alone is insufficient). Each returned luminosity is a pointwise coordinate-volume integrand; total luminosity integrates over coordinate volume without a second `sqrt(gamma)`.
- `NRPyLeakage_optical_depths_PathOfLeastResistance` expects `dxx[3]` in `(x, y, z)` order, metric stencils ordered `[minus-one, center, plus-one]`, finite positive cell widths and face-averaged diagonal metrics, and finite nonnegative opacities and optical depths.

`GRHayL/include/ghl_nrpyleakage.h` defines all physical constants, switches, unit macros, `robust_isnan`, `robust_isfinite`, and `NRPYLEAKAGE_FD_OR_RETURN`. On IEEE binary64 systems, `robust_isnan` and `robust_isfinite` use alias-safe `memcpy` to `uint64_t` with exponent/fraction bit masks, falling back to C99 predicates on non-binary64 targets. Because unsafe floating-point flags can corrupt preceding arithmetic, `configure` rejects `-ffast-math`, `-Ofast`, `-ffinite-math-only`, `-funsafe-math-optimizations`, `-fassociative-math`, `-freciprocal-math`, `-fno-signed-zeros`, `-fno-trapping-math`, Intel `-ftz`, and Intel fast floating-point models, probes predefined fast/finite-only macros and gradual underflow, and appends `-fp-model=precise -no-ftz` for Intel LLVM (which direct builds and executable links must retain).

Evidence (normative rule): `configure::no-ftz`; corroboration: `ghl_nrpyleakage.h::robust_isfinite`.

### Error Behavior And Writeback Guarantees

- **Fermi-Dirac helper:** `NRPyLeakage_Fermi_Dirac_integrals` zeroes `*Fermi_Dirac_integral` before branch dispatch (`z > 1e-3` vs. `z <= 1e-3`), supports keys `0` through `5`, and returns `ghl_error_invalid_fermi_dirac_integral_key` for unsupported keys (leaving output zeroed). `Unit_Tests/unit_test_code_error.c` checks invalid keys on both branches.
- **HDF5 and EOS errors:** Under `GHL_DISABLE_HDF5`, EOS-dependent routines return `ghl_error_used_disabled_hdf5` before table access. With HDF5 enabled, they call `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T` and propagate any EOS error unchanged before writing outputs.
- **Nucleon-blocking and shifted-moment errors:** Source-private `NRPyLeakage_nucleon_blocking.h` requires finite positive cgs density and temperature and finite `X_n`, `X_p` in `[0, 1]`, normalizing only endpoint excursions within the eight-corner interpolation forward-error bound. Both-zero composition succeeds neutrally; an exactly-one-zero state preserves occupied-species scattering and applies the analytic zero limit to charged-current products without forming divergent degeneracy differences. Larger excursions, nonfinite kinetic-degeneracy inversions or overlaps, or invalid/nonfinite shifted beta-moment inputs or results (including overflow from finite reaction shifts) return `ghl_error_nrpyleakage_blocking` before public output writeback.
- **Nonfinite-output sanitization:** After ordinary writeback, every EOS-dependent routine sanitizes public outputs, replacing nonfinite opacities with the converted cgs inverse-length floor and nonfinite luminosities or signed sources with `0.0`. Any replacement returns `ghl_error_nrpyleakage_nonfinite_output` while keeping finite sibling outputs (including finite negative values) intact. Intermediate `EnsureFinite` guards can substitute values without raising this error, so `ghl_success` does not guarantee absence of intermediate fallbacks.
- **Optical-depth update:** Returns `void` with no validation or error channel.

No-HDF5 builds still compile all guarded NRPyLeakage implementation files so `ghl_error_used_disabled_hdf5` paths remain linkable; `configure` excludes only the table-backed NRPyLeakage tests while keeping `unit_test_nrpyleakage_physics.c` and `unit_test_nrpyleakage_classifier_fallback.c`. Test-only EOS setup calls (`ghl_initialize_tabulated_eos_functions_and_params`, `ghl_tabulated_compute_eps_from_T`, `ghl_tabulated_compute_T_from_eps`) live in `Unit_Tests/nrpyleakage_main.h` and `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`.

Evidence (descriptive behavior): `NRPyLeakage_Fermi_Dirac_integrals.c::NRPyLeakage_Fermi_Dirac_integrals`; corroboration: `ghl.h::ghl_error_nrpyleakage_blocking`.

## Sources

- [GRHayL/include/ghl_radiation.h](../../../GRHayL/include/ghl_radiation.h) - `ghl_neutrino_luminosities`, `ghl_neutrino_opacities`, `ghl_neutrino_optical_depths`
- [GRHayL/include/ghl_nrpyleakage.h](../../../GRHayL/include/ghl_nrpyleakage.h) - public declarations, unit/physical constants, `robust_isnan`, `robust_isfinite`, `NRPYLEAKAGE_FD_OR_RETURN`
- [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) - `ghl_eos_parameters`, `ghl_error_codes_t`, `ghl_initialize_tabulated_eos_functions_and_params`
- [GRHayL/include/ghl_nrpyeos_tabulated.h](../../../GRHayL/include/ghl_nrpyeos_tabulated.h) - `hdf5.h` inclusion and tabulated EOS definitions
- [GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h) - EOS function-pointer dispatch
- [GRHayL/include/ghl_eos_functions_declaration.h](../../../GRHayL/include/ghl_eos_functions_declaration.h) - `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c) - Fermi-Dirac branch dispatch and invalid-key zeroing
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h) - private blocking validation and error checks
- [GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c) - tabulated chemistry callback registration
- [configure](../../../configure) - floating-point flag rejection, gradual-underflow probe, `--disable-hdf5` filtering
- [Unit_Tests/nrpyleakage_main.h](../../../Unit_Tests/nrpyleakage_main.h) - test-only EOS table initialization
- [Unit_Tests/unit_test_code_error.c](../../../Unit_Tests/unit_test_code_error.c) - invalid Fermi-Dirac key checks

## See Also

- Depends on: [EOS Inputs And Chemical Potentials](eos-inputs-and-chemical-potentials.md)
- Depends on: [Physics Model And Sources](physics-model-and-sources.md)
- Implements: [Implementation Flow](implementation-flow.md)
- Validated by: [Neutrinos Tests And Reference Data](tests-and-reference-data.md)
- See also: [CompOSE EOS Leakage Callback](compose-eos-leakage-callback.md)
