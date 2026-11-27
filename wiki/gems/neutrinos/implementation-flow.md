# Neutrinos Implementation Flow

> Built `GRHayL/Neutrinos/NRPyLeakage/` source files, minimal link sets, shared error and sanitizer boundaries, and per-file execution flows. · Status: confirmed
> Up: [Neutrinos](index.md)

## Summary

`GRHayL/Neutrinos/make.code.defn` delegates to `GRHayL/Neutrinos/NRPyLeakage/make.code.defn`, which compiles every public `NRPyLeakage_*` translation unit and tracks the source-private headers (`NRPyLeakage_nucleon_blocking.h`, `NRPyLeakage_rate_helpers.h`) via `#! INCS`. Checked-in C files are authoritative over the external ancestral notebooks ([Generator Provenance](generator-provenance.md)); blocking equations and dense-matter qualification are defined in [Nucleon Blocking And EOS Conventions](nucleon-blocking-and-eos-conventions.md).

## Detail

### Build Manifest And Smallest Link Sets

Every `.c` file listed in `GRHayL/Neutrinos/NRPyLeakage/make.code.defn` matches a public declaration in `GRHayL/include/ghl_nrpyleakage.h`, with no unbuilt `.c` files in the directory. Private headers `NRPyLeakage_nucleon_blocking.h` and `NRPyLeakage_rate_helpers.h` are not installed publicly and serve every EOS-dependent entry point. Individual operations require narrower link sets than the full module:

| Operation | Required implementation files | Additional GRHayL dependencies | Not required to link |
| --- | --- | --- | --- |
| Fermi-Dirac approximation | `NRPyLeakage_Fermi_Dirac_integrals.c` | Radiation/leakage headers, error enum, `-lm` | EOS, HDF5, opacity, depth, source, luminosity files |
| Standalone opacities | `NRPyLeakage_compute_neutrino_opacities.c`, `NRPyLeakage_nucleon_blocking.h`, `NRPyLeakage_rate_helpers.h`, `NRPyLeakage_Fermi_Dirac_integrals.c` | Radiation structs, constants/macros, tabulated-EOS callback, HDF5 guard, `-lm` | Combined source/opacity, optical-depth, luminosity files |
| GRMHD sources + opacities | `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`, `NRPyLeakage_nucleon_blocking.h`, `NRPyLeakage_rate_helpers.h`, `NRPyLeakage_Fermi_Dirac_integrals.c` | Same EOS, type, constant, error/HDF5, `-lm` dependencies | Standalone-opacity, optical-depth, luminosity files |
| One optical-depth update | `NRPyLeakage_optical_depths_PathOfLeastResistance.c` | Opacity/depth typedefs, `-lm` | EOS, HDF5, Fermi helper, source, opacity, luminosity files |
| Pointwise luminosities | `NRPyLeakage_compute_neutrino_luminosities.c`, `NRPyLeakage_nucleon_blocking.h`, `NRPyLeakage_rate_helpers.h`, `NRPyLeakage_Fermi_Dirac_integrals.c` | Radiation structs, constants/macros, tabulated-EOS callback, metric inputs, `-lm` | Source, standalone-opacity, optical-depth files |

The optical-depth routine consumes opacity and neighbor-depth data without calling any opacity routine, and the combined source/opacity routine carries its own opacity write path rather than calling `NRPyLeakage_compute_neutrino_opacities`.

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/make.code.defn::SRCS`; corroboration: `GRHayL/Neutrinos/NRPyLeakage/make.code.defn::NRPyLeakage_nucleon_blocking.h`.

### Shared Failure And Finite-Value Boundaries

Every EOS-dependent routine executes the same error sequence before writing public outputs: (1) return `ghl_error_used_disabled_hdf5` under `GHL_DISABLE_HDF5`, (2) call `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T` and return any EOS error unchanged, (3) convert `rho` to cgs via `NRPyLeakage_units_geom_to_cgs_D` and run `NRPyLeakage_nucleon_blocking.h` validation and inversion (returning `ghl_error_nrpyleakage_blocking` on out-of-range inputs, nonfinite results, or shifted-beta-moment failures), and (4) evaluate Fermi-Dirac integrals via `NRPYLEAKAGE_FD_OR_RETURN`. All of these exits leave caller output structs untouched.

After roundoff normalization of `X_n` and `X_p`, both-zero composition takes a neutral path, while an exactly-one-zero state retains occupied-species scattering and non-beta channels and applies the analytic zero limit to charged-current products without constructing divergent $q$. Final-output sanitizers run after ordinary writeback: `nrpyl_sanitize_opacities` replaces nonfinite opacities with the converted cgs inverse-length floor, while `nrpyl_sanitize_sources` and `nrpyl_sanitize_luminosities` replace nonfinite sources or luminosities with `0.0`, returning `ghl_error_nrpyleakage_nonfinite_output` while keeping finite siblings intact. Local `EnsureFinite` wrappers replace selected nonfinite intermediate expressions with a small positive value without setting `ghl_error_nrpyleakage_nonfinite_output`.

All finite-value guards rely on `robust_isfinite` / `robust_isnan` (alias-safe IEEE binary64 `memcpy` bit classification with C99 fallback) paired with `configure`'s rejection of unsafe math / flush-to-zero flags (direct and Cactus builds must preserve `-fp-model=precise -no-ftz` on Intel LLVM).

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h::nrpyl_sanitize_opacities`; corroboration: `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h::nrpyl_sanitize_luminosities`.

### EOS-Dependent Routine Flows

- **`NRPyLeakage_compute_neutrino_luminosities.c` (`NRPyLeakage_compute_neutrino_luminosities`):** Runs the shared HDF5/EOS/blocking/shifted-beta-moment sequence, evaluates generated emissivity, opacity-like denominator, Fermi-Dirac, and optical-depth suppression blocks, applies `alpha`, spatial metric determinant, `W`, and `NRPyLeakage_units_cgs_to_geom_Q` scaling, writes `lum.nue`, `lum.anue`, and `lum.nux`, and runs `nrpyl_sanitize_luminosities`. Covered by `Unit_Tests/unit_test_nrpyleakage_luminosities.c` (row-matched unperturbed/perturbed binary replay aborting with row index on mismatch).
- **`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c` (`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`):** Runs the shared HDF5/EOS/blocking/shifted-beta-moment sequence, evaluates generated rate, opacity, and optical-depth-limited source blocks, writes `*R_source`, `*Q_source`, and all six `kappa` entries (`nue[0..1]`, `anue[0..1]`, `nux[0..1]`), and runs both `nrpyl_sanitize_sources` and `nrpyl_sanitize_opacities` before returning `ghl_error_nrpyleakage_nonfinite_output` if either group was sanitized. Covered by `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c` (RK4 evolution of `Y_e` and `eps`, checking every EOS lookup with `ghl_abort_if_error`; opacity outputs are executed there but not compared against reference binaries).
- **`NRPyLeakage_compute_neutrino_opacities.c` (`NRPyLeakage_compute_neutrino_opacities`):** Runs the shared HDF5/EOS/blocking sequence, evaluates $\nu_e$ and $\bar\nu_e$ scattering Fermi integrals first, and passes them into `nrpyl_compute_beta_absorption_moments_with_fermi` as $F_2$ and $F_3$ normalizers. Any channel whose neutrino-energy threshold $\max(0, -\mathrm{sign}\cdot q/T)$ is zero reuses all four integrals directly as its shifted moments (guaranteed for at least one channel by the sign of $q$, and both when $q = 0$); nonzero-threshold channels evaluate shifted integrals normally. It then runs the remaining opacity formula blocks, writes all six `kappa` entries, and applies `nrpyl_sanitize_opacities`. Covered by `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c` (serial mismatch check reporting grid coordinates; source-term tests do not cover this file).

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c::NRPyLeakage_compute_neutrino_opacities`; corroboration: `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`.

### Table-Free Routine Flows

- **`NRPyLeakage_Fermi_Dirac_integrals.c` (`NRPyLeakage_Fermi_Dirac_integrals`):** Zeroes `*Fermi_Dirac_integral`, dispatches strictly on `z > 1e-3` (high branch) vs. `z <= 1e-3` (low branch), evaluates keys `0` through `5`, returns `ghl_success` without checking null pointers or output finiteness, and returns `ghl_error_invalid_fermi_dirac_integral_key` (leaving output `0.0`) for other keys. Covered across all valid keys in `unit_test_nrpyleakage_physics.c`, selected keys in `unit_test_nrpyleakage_luminosities.c`, and invalid keys on both branches in `unit_test_code_error.c`.
- **`NRPyLeakage_optical_depths_PathOfLeastResistance.c` (`NRPyLeakage_optical_depths_PathOfLeastResistance`):** Accepts `dxx[0..2]` (`x, y, z`), three-point metric stencils `[minus-one, center, plus-one]` (`stencil_gxx`, `stencil_gyy`, `stencil_gzz`), six neighbor opacity and depth structs (`im1`, `ip1`, `jm1`, `jp1`, `km1`, `kp1`), local opacity, and output depth. It averages face diagonal metrics to compute six face path lengths, averages local and neighbor opacities on each face, forms six candidate optical depths per species and slot, selects the minimum via nested `fmin`, and writes all six `tau_i_j_k` entries.
  Because it returns `void` with no input validation, negative face-averaged metrics enter `sqrt` unguarded and nested `fmin` can mask a NaN candidate if another candidate is finite. Covered by `unit_test_nrpyleakage_constant_density_sphere.c` (including `check_sphere_caller_neighbor_order` with asymmetric x-face metrics) and the all-direction asymmetric stencil in `unit_test_nrpyleakage_physics.c`.

Evidence (descriptive behavior): `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c::NRPyLeakage_optical_depths_PathOfLeastResistance`; corroboration: `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c::NRPyLeakage_Fermi_Dirac_integrals`.

## Sources

- [GRHayL/Neutrinos/make.code.defn](../../../GRHayL/Neutrinos/make.code.defn) - `SUBDIRS` routing to `NRPyLeakage`
- [GRHayL/Neutrinos/NRPyLeakage/make.code.defn](../../../GRHayL/Neutrinos/NRPyLeakage/make.code.defn) - `SRCS`, `#! INCS` private header tracking
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c) - `NRPyLeakage_compute_neutrino_luminosities`
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c) - `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c) - `NRPyLeakage_compute_neutrino_opacities`
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c) - `NRPyLeakage_Fermi_Dirac_integrals`
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c) - `NRPyLeakage_optical_depths_PathOfLeastResistance`
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h) - `nrpyl_compute_beta_absorption_moments_with_fermi`
- [GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h](../../../GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h) - `nrpyl_sanitize_opacities`, `nrpyl_sanitize_sources`, `nrpyl_sanitize_luminosities`
- [GRHayL/include/ghl_nrpyleakage.h](../../../GRHayL/include/ghl_nrpyleakage.h) - `robust_isfinite`, `NRPYLEAKAGE_FD_OR_RETURN`
- [Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c](../../../Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c) - `check_sphere_caller_neighbor_order`
- [Unit_Tests/unit_test_nrpyleakage_physics.c](../../../Unit_Tests/unit_test_nrpyleakage_physics.c) - table-free routine checks
- [Tabulated_EOS_IllinoisGRMHD Repository](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD) - ancestral NRPy+ notebooks
- [Tutorial-Leakage_Scheme-Implementation.ipynb](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Leakage_Scheme-Implementation.ipynb) - symbolic C-generation provenance

## See Also

- Depends on: [Neutrinos API And Data](api-and-data.md)
- Depends on: [Nucleon Blocking And EOS Conventions](nucleon-blocking-and-eos-conventions.md)
- See also: [Generator Provenance](generator-provenance.md)
- Validated by: [Neutrinos Tests And Reference Data](tests-and-reference-data.md)
