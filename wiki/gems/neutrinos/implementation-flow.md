# Neutrinos Implementation Flow

This page maps the built NRPyLeakage implementation files listed by
`GRHayL/Neutrinos/NRPyLeakage/make.code.defn`. Source remains authority; the
large generated formula blocks are described by role instead of copied.

Use [Generator Provenance](generator-provenance.md) to trace notebook source
and regeneration limits. Use
[Physics And EOS Contract](physics-and-eos-contract.md) when replacing the
GRHayL tabulated EOS or interpreting leakage outputs.

## Build List

`GRHayL/Neutrinos/make.code.defn` routes Neutrinos to the `NRPyLeakage`
subdirectory. `GRHayL/Neutrinos/NRPyLeakage/make.code.defn` builds exactly
these source files:

- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_luminosities.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_compute_neutrino_opacities.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_Fermi_Dirac_integrals.c`
- `GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c`

All five names have matching public declarations in `ghl_nrpyleakage.h`; no
extra Neutrinos `.c` file sits outside the manifest.

## Smallest File Sets And Data Dependencies

The five manifest entries do not form one indivisible link unit. Current
source has these narrower boundaries:

| Requested operation | Required implementation files | Additional current-GRHayL dependencies | Not required by that entry point |
| --- | --- | --- | --- |
| Fermi-Dirac approximation | `NRPyLeakage_Fermi_Dirac_integrals.c` | Radiation/leakage headers, GRHayL error enum, math functions | EOS, HDF5, opacity, depth, source, luminosity files |
| Opacities only | `NRPyLeakage_compute_neutrino_opacities.c` and `NRPyLeakage_Fermi_Dirac_integrals.c` | Radiation structs, leakage constants/macros, tabulated-EOS callback, GRHayL errors/HDF5 guard, math functions | Combined source/opacity, optical-depth, and luminosity files |
| GRMHD sources plus opacities | `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c` and `NRPyLeakage_Fermi_Dirac_integrals.c` | Same EOS, type, constant, error/HDF5, and math adapter | Standalone-opacity, optical-depth, and luminosity files |
| One path-of-least-resistance depth update | `NRPyLeakage_optical_depths_PathOfLeastResistance.c` | Opacity/depth struct definitions and math functions | EOS, HDF5, Fermi helper, source, standalone-opacity, and luminosity files |
| Pointwise luminosities | `NRPyLeakage_compute_neutrino_luminosities.c` and `NRPyLeakage_Fermi_Dirac_integrals.c` | Radiation structs, leakage constants/macros, tabulated-EOS callback, GRHayL errors/HDF5 guard, metric inputs, math functions | Source, standalone-opacity, and optical-depth implementation files |

These are link and call boundaries, not a complete simulation workflow. The
optical-depth routine needs opacity and neighbor-depth *data*, but does not
call the standalone-opacity routine. Source/opacity and luminosity routines
consume optical-depth data, but do not call the optical-depth implementation.
The combined routine owns a separate generated opacity write path; it does not
delegate to `NRPyLeakage_compute_neutrino_opacities`.

For a host with its own EOS and hydrodynamics, retain the Fermi helper wherever
generated rate code calls `NRPYLEAKAGE_FD_OR_RETURN`; port the optical-depth
file independently if its six-neighbor update is wanted. The exact ABI,
callback, error, HDF5, and unit substitutions are mapped in
[API And Data](api-and-data.md).

## Generator Versus Current Source Authority

The separate
[Tabulated_EOS_IllinoisGRMHD repository](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD)
contains NRPy/Python development material absent from this tree.
[Generator Provenance](generator-provenance.md) preserves the notebook roles,
output paths, unit construction, `tmp_*` origin, and reviewed generator
coverage.
Use it for ancestry only. Current repo-local C, headers, manifests, and tests
remain authority; generate into a disposable directory, compare formulas and
constants, and deliberately reapply the GRHayL ABI, EOS, error, and unit
adapter before replacing any current file.

## Shared Failure Boundary

The three EOS-dependent routines return immediately with
`ghl_error_used_disabled_hdf5` in no-HDF5 builds. With HDF5, they return the
tabulated EOS error unchanged, and generated Fermi calls return an invalid-key
error through `NRPYLEAKAGE_FD_OR_RETURN`. Output writes occur only after those
calls, so these error exits leave caller outputs unchanged. None checks null
pointers, EOS initialization, table bounds independently of the EOS call, or
physical/finiteness preconditions.

`robust_isfinite` is used only by the combined source-term file;
`robust_isnan` has no active repo-local caller. Both public inline helpers
inspect a `double` by casting its address to `unsigned long *`, with no size or
representation check. Treat them as current platform-dependent implementation,
not a portable finiteness contract. Tests do not inject NaN/Inf into any
leakage routine or directly exercise these helpers.

## `NRPyLeakage_compute_neutrino_luminosities.c`

Public routine: `NRPyLeakage_compute_neutrino_luminosities`.

Flow:

1. Return `ghl_error_used_disabled_hdf5` when `GHL_DISABLE_HDF5` is set.
2. Query tabulated EOS composition and chemical potentials with
   `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`; return any EOS
   error directly.
3. Convert `rho` to cgs units through `NRPyLeakage_units_geom_to_cgs_D`, then
   derive proton/neutron fractions used by the leakage rates.
4. Run source-owned generated formula blocks for emissivity, opacity-like
   denominators, Fermi-Dirac factors, and optical-depth suppression. These
   blocks call `NRPYLEAKAGE_FD_OR_RETURN`, so Fermi-Dirac errors propagate.
5. Use metric/lapse/Lorentz input in the luminosity scaling:
   `alpha`, the six spatial metric components, and `W` enter the
   `NRPyLeakage_units_cgs_to_geom_Q` luminosity prefactor before writeback;
   `NRPyLeakage_units_geom_to_cgs_D` is the earlier density conversion.
6. Write `lum->nue`, `lum->anue`, and `lum->nux`.

Finite handling: local `EnsureFinite` wraps selected generated subexpressions
with `isfinite` fallback to a small positive value; there is no final output
array scrub in this file.

Nearest tests: `Unit_Tests/unit_test_nrpyleakage_luminosities.c` directly
checks selected Fermi-Dirac branches, generates luminosity fixtures, recomputes
`NRPyLeakage_compute_neutrino_luminosities`, and reads `nue`, `anue`, and `nux`
fixtures. Its three `ghl_pert_test_fail` return values are discarded, so
numerical luminosity mismatches do not currently fail the test.

## `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`

Public routine:
`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`.

Flow:

1. Return `ghl_error_used_disabled_hdf5` when `GHL_DISABLE_HDF5` is set.
2. Query tabulated EOS composition and chemical potentials with
   `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`; return any EOS
   error directly.
3. Convert `rho` to cgs units, then build proton/neutron fraction inputs.
4. Run source-owned generated formula blocks for rate terms, opacity terms,
   optical-depth limited source terms, and Fermi-Dirac factors. Calls through
   `NRPYLEAKAGE_FD_OR_RETURN` propagate invalid Fermi-Dirac keys if they ever
   occur.
5. Write `*R_source` and `*Q_source`.
6. Write all six opacity entries: `kappa->nue[0..1]`,
   `kappa->anue[0..1]`, and `kappa->nux[0..1]`.

Finite handling: this file's `EnsureFinite` uses `robust_isfinite` from
`GRHayL/include/ghl_nrpyleakage.h` and replaces non-finite intermediate terms
with a small positive value. There is no final output array scrub after
`R_source`, `Q_source`, or opacity writes.

Nearest tests: `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c` calls
this routine in its RHS, divides `R_source` and `Q_source` by `rho`, advances
`Y_e` and `eps` with RK4, and reads fixture replay. Its comparison helper also
discards every `ghl_pert_test_fail` result, so the test supplies execution and
file-shape evidence, not an effective numerical assertion. Opacity writes get
execution coverage there through the same call but are not compared.

## `NRPyLeakage_compute_neutrino_opacities.c`

Public routine: `NRPyLeakage_compute_neutrino_opacities`.

Flow:

1. Return `ghl_error_used_disabled_hdf5` when `GHL_DISABLE_HDF5` is set.
2. Query tabulated EOS composition and chemical potentials with
   `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`; return any EOS
   error directly.
3. Convert `rho` to cgs units through `NRPyLeakage_units_geom_to_cgs_D`.
4. Build the proton/neutron fraction inputs and run source-owned generated
   formula blocks for absorption/scattering opacity entries. These generated
   blocks call `NRPYLEAKAGE_FD_OR_RETURN`, so Fermi-Dirac errors propagate.
5. Write all six opacity entries: `kappa->nue[0..1]`,
   `kappa->anue[0..1]`, and `kappa->nux[0..1]`.
6. Scrub each written opacity entry with `isfinite`; any non-finite final
   value is reset to a small positive value.

Finite handling: local `EnsureFinite` handles selected generated
subexpressions, then the final loop handles non-finite output entries.

Nearest tests: `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
directly calls this routine for interior and exterior states, stores the six
opacity fields on the grid, and reads opacity/depth fixtures. Its comparison
return values are discarded, so those fixture values cannot currently fail the
test. Source-term tests do not cover this implementation: the combined
source-term routine has its own opacity write path.

## `NRPyLeakage_Fermi_Dirac_integrals.c`

Public routine: `NRPyLeakage_Fermi_Dirac_integrals`.

Flow:

1. Set `*Fermi_Dirac_integral = 0.0` before dispatch.
2. Use the high-`z` branch when `z > 1e-3`; otherwise use the low-`z` branch.
3. Support keys `0`, `1`, `2`, `3`, `4`, and `5` in both branches.
4. Overwrite the output pointer for supported keys and return `ghl_success`.
5. Return `ghl_error_invalid_fermi_dirac_integral_key` for unsupported keys.
   Because the output is zeroed before dispatch, invalid-key calls leave the
   caller-provided output at zero.

Generated formula role: each branch contains source-owned approximation
expressions for the selected key; do not duplicate those expressions into KB
pages.

Nearest tests: `Unit_Tests/unit_test_nrpyleakage_luminosities.c` checks valid
keys `0`, `1`, and `2` in the high-`z` branch and keys `0` and `1` in the
low-`z` branch. Keys `3` through `5` lack direct valid-result assertions.
`Unit_Tests/unit_test_code_error.c`
directly checks invalid-key behavior for both `z < 1e-3` and `z > 1e-3`, and
maps those cases to `ghl_error_invalid_fermi_dirac_integral_key`.

The helper performs no null or finiteness check. Supported keys return
`ghl_success` even if extreme `z` makes approximation arithmetic non-finite.
The branch boundary itself is exact: `z > 1e-3` uses the high branch;
`z <= 1e-3` uses the low branch.

## `NRPyLeakage_optical_depths_PathOfLeastResistance.c`

Public routine: `NRPyLeakage_optical_depths_PathOfLeastResistance`.

Flow:

1. Accept cell widths `dxx[0..2]`, three-point metric stencils
   `stencil_gxx`, `stencil_gyy`, and `stencil_gzz`, six neighbor opacity
   structs, six neighbor optical-depth structs, current-cell opacity, and the
   current-cell optical-depth output. Each metric stencil uses
   `[minus-one, center, plus-one]`; `dxx` uses coordinate order `x`, `y`, `z`.
2. Average same-direction metric stencil entries to get face-centered
   diagonal metric components.
3. Convert face metrics and cell widths into six face path lengths.
4. Average current-cell and neighbor opacities at each face for each species
   and each two-entry index.
5. Build six candidate optical depths per species/index by adding each
   neighbor depth to path length times face-averaged opacity.
6. Select the minimum candidate with nested `fmin` calls for each
   species/index.
7. Write `tau_i_j_k->nue[0..1]`, `tau_i_j_k->anue[0..1]`, and
   `tau_i_j_k->nux[0..1]`.

The routine returns `void` and has no failure channel. It does not validate
pointer lengths, metric sign, grid spacing, opacity sign, or finiteness.
Negative face-averaged diagonal metric values enter `sqrt`; any resulting
non-finite propagation is not scrubbed before output.

Input shape note: the six-neighbor shape is explicit in the signature:
`im1`, `ip1`, `jm1`, `jp1`, `km1`, and `kp1` opacity/depth pointers surround
the current cell.

Nearest tests: `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
directly computes opacities, iterates optical-depth updates with flat metric
stencils, calls `NRPyLeakage_optical_depths_PathOfLeastResistance`, writes all
six output depth fields back to grid storage, and reads fixture replay. That
caller passes each neighbor pair in plus-then-minus order, opposite the public
minus-then-plus signature. Its flat metric and minimum over symmetric
directions hide this reversal, so current test evidence does not verify
directional argument mapping for unequal plus/minus metrics. It also discards
all comparison results, as noted above.

## Ground Truth References

- Original NRPyLeakage development repository:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD
- NRPyLeakage symbolic implementation and C-generation notebook:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Leakage_Scheme-Implementation.ipynb
