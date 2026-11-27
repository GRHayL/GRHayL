# Neutrinos Generator Provenance

> Ancestral NRPy+ notebooks in `Tabulated_EOS_IllinoisGRMHD`, `tmp_*` common-subexpression origin, hand-applied C corrections, and ILEAS/FDINT attribution. · Status: confirmed
> Up: [Neutrinos](index.md)

## Summary

No Python or NRPy+ generator for NRPyLeakage resides under `GRHayL/Neutrinos/`; ancestral symbolic notebooks and generated outputs live in Leo Werneck's external `Tabulated_EOS_IllinoisGRMHD` repository. Current checked-in GRHayL C, headers, manifests, and tests govern all compiled behavior, including hand-applied physics and numerical corrections not synchronized back to those notebooks.

## Detail

### Authority Order And External Notebook Roles

Treat sources in strict precedence: (1) repo-local GRHayL C, headers, manifests, and tests govern compiled behavior; (2) the external implementation notebook provides provenance for symbolic construction and initial generation; (3) the formulation and validation notebooks provide derivation cross-checks; and (4) external checked-in C is historical comparison material only, never a drop-in replacement. Because the upstream `master` branch can move and carries no pinned Python lockfile, durable facts are recorded here:

- **`Tutorial-Neutrino_leakage-Formulation.ipynb`:** Derives beta processes, pair annihilation, transverse plasmon decay, bremsstrahlung, opacities, optical depths, diffusion times, effective rates, and `R`/`Q` sources. Also discusses neutrino pressure, which is outside GRHayL's current build.
- **`Tutorial-Leakage_Scheme-Implementation.ipynb`:** Imports vendored `nrpy_core`, SymPy, and Astropy constants; builds symbolic Fermi-Dirac fits of orders `0` through `5`, free rates, densities, opacities, diffusion times, and sources; emits C via `outputC.outputC`, `outputC.outCfunction`, and `optical_depth_helpers`. Geometric units are built from `M = M_sun`, `L = G M / c^2`, `time = L / c`, `D = M / L^3`, `R = D / time`, and `Q = c^2 D / time` (converted to `MeV / (cm^3 s)`), while current numeric constants live exclusively in `GRHayL/include/ghl_nrpyleakage.h`.
- **`Tutorial-Neutrino_leakage-Numerical_implementation.ipynb`:** Annotated walkthrough of the older Fortran/C++ `ZelmaniLeak` thorn (drivers, interpolation, deleptonization profiles, neutrino pressure, scheduling); comparative evidence only, not the NRPy+ generator.
- **`Tutorial-Neutrino_leakage-Optically_thin_semi_analytic_calculation.ipynb`:** Independent semi-analytic optically thin evolution cross-check illuminating the intent of `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c`.

### Preserved Generation Flow And `tmp_*` Origin

In `Tutorial-Leakage_Scheme-Implementation.ipynb`, `generate_ET_thorn` selects output paths: `false` writes `standalone/NRPyLeakage/`, whereas `true` writes C, header, and `make.code.defn` into `standalone/NRPyLeakageET/src/` plus `param.ccl`, `schedule.ccl`, and `interface.ccl` one level up. Initialization recursively deletes and recreates `standalone/<module>`, so any regeneration must run in a disposable clone. Its driver generates in order: (1) Fermi-Dirac integrals, (2) low-level combined source/opacity kernels for NRPy and HARM constant sets, (3) combined selector driver, (4) grid-level optical depth, (5) low-level opacity kernels for both constant sets, (6) opacity driver, (7) header, and (8) ET CCL/build files when enabled. Checked-in upstream directories (`standalone/NRPyLeakage/` and `NRPyLeakageET/src/`) do not form a turnkey reproducible lockfile sequence; results vary with SymPy, Astropy, and NRPy+ versions.

Upstream `nrpy_core/outputC.py` enables SymPy common-subexpression elimination (`cse`) by default with prefix `tmp` starting at `tmp_0`, hoisting repeated expressions into `const` temporaries before final assignments. Individual `tmp_N` variables have no physical meaning and vary across symbolic ordering and SymPy versions; read named final assignments or notebook expressions instead. (`GRHayL/Flux_Source/generate_flux_source.py` uses NRPy 2 independently and does not establish NRPyLeakage's original settings.)

Evidence (generated evidence): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c::tmp_0`; corroboration: `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c::tmp_0`.

### GRHayL Adaptations And Hand-Applied C Corrections

GRHayL adapts the ancestral outputs into struct-based public calls with `ghl_error_codes_t` returns, tabulated EOS dispatch, `GHL_DISABLE_HDF5` returns, `NRPYLEAKAGE_FD_OR_RETURN`, and final sanitizers:

- `NRPyLeakage_Fermi_Dirac_integrals.c` returns `ghl_error_codes_t` via output pointer instead of aborting on unsupported keys.
- `NRPyLeakage_compute_neutrino_opacities.c` and `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c` retain recognizable `tmp_*` blocks from the standalone NRPy-constant kernels.
- `NRPyLeakage_compute_neutrino_luminosities.c` corresponds to upstream `NRPyLeakageET_compute_neutrino_luminosities_nrpy_constants.c` (no luminosity generator notebook exists upstream).
- `NRPyLeakage_optical_depths_PathOfLeastResistance.c` extracts the pointwise face-average minimum-path calculation from the ET routine.
- `NRPyLeakage_nucleon_blocking.h` and `NRPyLeakage_rate_helpers.h` are hand-maintained GRHayL headers tracked via `#! INCS`.

`GRHayL/include/ghl_nrpyleakage.h` exposes one active NRPy constant family and retains dormant `NRPyLeakage_ZL_*` macros unused by current C; GRHayL has no runtime NRPy/HARM constant-set selector.

Because no repo-local generator reproduces the complete file set, the following corrections were applied directly to checked-in C and must be preserved across any future regeneration: (1) shared cgs diffusion prefactor `6.0 / NRPyLeakage_c_light`, (2) `rho_cgs * rho_cgs` scaling in nucleon-nucleon bremsstrahlung, (3) single heavy-lepton `Q_free_nux` shared by numerator and inverse loss time, (4) stable `log1p` form for order-zero Fermi-Dirac integral, (5) density-derived free-nucleon blocking with `expm1` overlap algebra, and (6) paired shifted charged-current emission/absorption moments using $q = \widehat\mu - T(\eta_n - \eta_p)$ from a single Kirchhoff-balanced spectral parent. To regenerate a formula safely: edit symbolically in a disposable clone, diff against current C, reapply all GRHayL interfaces and hand-applied corrections, and re-verify tests.

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h::nrpyl_effective_emission_rate`; corroboration: `GRHayL/include/ghl_nrpyleakage.h::NRPyLeakage_ZL_Q_npmass`.

### Nucleon-Blocking Literature And Code Attribution

`NRPyLeakage_nucleon_blocking.h` follows Ardevol-Pulpillo et al. (ILEAS, *MNRAS* 485, 2019, `doi:10.1093/mnras/stz613`) Appendix B Eqs. (69)-(71) for density-derived free-nucleon degeneracies and Appendix C Eqs. (100)-(109) for thresholded shifted reaction pairing (implemented as an algebraic grey reduction without finite charged-lepton mass, EOS effective masses, or mean-field shifts). Scalar `fdm1h` and `ifd1h` rational fits are ported from Scott Maddox's `fdint/_fdint.pyx` (implementing T. Fukushima's 2015 minimax fits) with FDINT's full BSD-3 license retained in the header. Common density normalization, fraction ordering, `expm1` overlap, and equal-population limits are GRHayL adaptations relying on Sterbenz's exact-subtraction theorem (P. H. Sterbenz, *Floating-Point Computation*, 1974, Sec. 4.3).

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h::nrpyl_compute_population`; corroboration: `GRHayL/Neutrinos/NRPyLeakage/make.code.defn::NRPyLeakage_rate_helpers.h`.

## Sources

- [GRHayL/Neutrinos/NRPyLeakage/make.code.defn](../../../GRHayL/Neutrinos/NRPyLeakage/make.code.defn) - built C files and `#! INCS` private headers
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c) - `tmp_0` CSE sequence
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c) - `tmp_0` CSE sequence
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h) - `nrpyl_compute_population`, BSD-3 attribution
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h) - `nrpyl_effective_emission_rate`
- [GRHayL/include/ghl_nrpyleakage.h](../../../GRHayL/include/ghl_nrpyleakage.h) - active NRPy constants and dormant `NRPyLeakage_ZL_Q_npmass`
- [GRHayL/Flux_Source/generate_flux_source.py](../../../GRHayL/Flux_Source/generate_flux_source.py) - separate NRPy 2 generator
- [Tabulated_EOS_IllinoisGRMHD Repository](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD) - ancestral repository
- [Tutorial-Neutrino_leakage-Formulation.ipynb](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Formulation.ipynb) - physics derivation notebook
- [Tutorial-Leakage_Scheme-Implementation.ipynb](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Leakage_Scheme-Implementation.ipynb) - symbolic generator notebook
- [Tutorial-Neutrino_leakage-Numerical_implementation.ipynb](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Numerical_implementation.ipynb) - ZelmaniLeak walkthrough
- [Tutorial-Neutrino_leakage-Optically_thin_semi_analytic_calculation.ipynb](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Neutrino_leakage-Optically_thin_semi_analytic_calculation.ipynb) - optically thin semi-analytic check
- [standalone/NRPyLeakage Tree](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/tree/master/standalone/NRPyLeakage) - checked-in standalone C outputs
- [NRPyLeakageET/src Tree](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/tree/master/NRPyLeakageET/src) - checked-in ET thorn C outputs
- [nrpy_core Tree](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/tree/master/nrpy_core) - `outputC.py` and `optical_depth_helpers.py`
- [Ardevol-Pulpillo et al. (ILEAS, MNRAS 485)](https://doi.org/10.1093/mnras/stz613) - Appendices B and C
- [FDINT _fdint.pyx](https://github.com/scott-maddox/fdint/blob/master/fdint/_fdint.pyx) - ported `fdm1h` and `ifd1h` fits
- [Fukushima (2015) Half-Odd Fermi Integrals](https://doi.org/10.1016/j.amc.2015.03.009) - forward fit reference
- [Fukushima (2015) Inverse Fermi Integrals](https://doi.org/10.1016/j.amc.2015.03.015) - inverse fit reference

## See Also

- Implements: [Implementation Flow](implementation-flow.md)
- See also: [Generated Boundaries](../../architecture/generated-boundaries.md)
- See also: [Nucleon Blocking And EOS Conventions](nucleon-blocking-and-eos-conventions.md)
- See also: [Physics Model And Sources](physics-model-and-sources.md)
