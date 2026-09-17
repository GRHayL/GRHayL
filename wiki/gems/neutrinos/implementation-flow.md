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

Every listed name has a matching public declaration in `ghl_nrpyleakage.h`; no
extra Neutrinos `.c` file sits outside the manifest.

The same manifest records `NRPyLeakage_nucleon_blocking.h` and
`NRPyLeakage_rate_helpers.h` through `#! INCS`. These source-private headers
are not installed public headers and add no public functions. Both headers
serve all three EOS-dependent routines.

## Smallest File Sets And Data Dependencies

The manifest entries do not form one indivisible link unit. Current
source has these narrower boundaries:

| Requested operation | Required implementation files | Additional current-GRHayL dependencies | Not required by that entry point |
| --- | --- | --- | --- |
| Fermi-Dirac approximation | `NRPyLeakage_Fermi_Dirac_integrals.c` | Radiation/leakage headers, GRHayL error enum, math functions | EOS, HDF5, opacity, depth, source, luminosity files |
| Opacities only | `NRPyLeakage_compute_neutrino_opacities.c`, `NRPyLeakage_nucleon_blocking.h`, `NRPyLeakage_rate_helpers.h`, and `NRPyLeakage_Fermi_Dirac_integrals.c` | Radiation structs, leakage constants/macros, tabulated-EOS callback, GRHayL errors/HDF5 guard, math functions | Combined source/opacity, optical-depth, and luminosity files |
| GRMHD sources plus opacities | `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`, `NRPyLeakage_nucleon_blocking.h`, `NRPyLeakage_rate_helpers.h`, and `NRPyLeakage_Fermi_Dirac_integrals.c` | Same EOS, type, constant, error/HDF5, and math adapter | Standalone-opacity, optical-depth, and luminosity files |
| One path-of-least-resistance depth update | `NRPyLeakage_optical_depths_PathOfLeastResistance.c` | Opacity/depth struct definitions and math functions | EOS, HDF5, Fermi helper, source, standalone-opacity, and luminosity files |
| Pointwise luminosities | `NRPyLeakage_compute_neutrino_luminosities.c`, `NRPyLeakage_nucleon_blocking.h`, `NRPyLeakage_rate_helpers.h`, and `NRPyLeakage_Fermi_Dirac_integrals.c` | Radiation structs, leakage constants/macros, tabulated-EOS callback, GRHayL errors/HDF5 guard, metric inputs, math functions | Source, standalone-opacity, and optical-depth implementation files |

These are link and call boundaries, not a complete simulation workflow. The
optical-depth routine needs opacity and neighbor-depth *data*, but does not
call the standalone-opacity routine. Source/opacity and luminosity routines
consume optical-depth data, but do not call the optical-depth implementation.
The combined routine owns a separate generated opacity write path; it does not
delegate to `NRPyLeakage_compute_neutrino_opacities`.

For a host with its own EOS and hydrodynamics, retain the private blocking and
rate headers and the Fermi helper for each EOS-dependent entry point. Port the
optical-depth file independently if its six-neighbor update is wanted. The
exact ABI, callback, error, HDF5, and unit substitutions are mapped in
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
tabulated EOS error unchanged. The blocking helper then validates finite,
positive cgs density and temperature; finite free fractions in `[0,1]`, apart
from endpoint excursions within the interpolation forward-error bound; and
finite, bounded inversion and overlap results. Failure returns
`ghl_error_nrpyleakage_blocking`. Generated Fermi calls return an invalid-key
error through `NRPYLEAKAGE_FD_OR_RETURN`. Output writes occur only after these
calls, so all these error exits leave caller outputs unchanged. The routines
do not check null pointers, EOS initialization, or table bounds independently
of the EOS call.

## Shared Nucleon-Blocking Evaluator

`NRPyLeakage_nucleon_blocking.h` reconstructs kinetic degeneracies from cgs
density, temperature, and EOS free-neutron/free-proton fractions. This removes
dependence on the arbitrary common energy zero of `mu_n` and `mu_p`, avoids
counting bound nucleons as free targets, and removes the old equal-population
quotient pole. One private helper keeps opacity, source, and luminosity paths
on the same blocking model instead of maintaining three copies.

The evaluator uses piecewise rational fits for inverse $F_{1/2}$ and
$F_{-1/2}$, an `expm1` overlap identity, and an analytic equal-population
limit. It adds no quadrature, iterative root solve, or EOS/table lookup. This
choice keeps blocking cost compatible with a leakage approximation while
retaining finite and population-bound checks. The source header records ILEAS,
FDINT/Fukushima, Sterbenz, and BSD-3 provenance; the exact equations and model
limits are in [Physics And EOS Contract](physics-and-eos-contract.md).

The callers retain `muhat` in the grey equilibrium-neutrino degeneracy and use
the helper's kinetic degeneracy difference to form the reaction shift
`q = muhat - T*(eta_n-eta_p)`. Private algebraic moment evaluators apply that
same shift and its particle-energy threshold to the paired charged-current
emission and absorption kernels. This enforces their spectral Kirchhoff
relation without quadrature. It remains a grey leakage approximation: emitted
neutrino vacancy is sampled at the mean energy, and the independent physical
qualification bounds the resulting model error.

All leakage finite-value guards use `robust_isfinite` or `robust_isnan`. On
IEEE binary64 platforms the public inline helpers copy the representation with
`memcpy` and classify exponent/fraction bits through `uint64_t`; this is
alias-safe and remains effective when user flags enable finite-math
assumptions. Compile-time representation guards select the C99 predicates on
other platforms. The table-free physics executable directly checks finite
values, infinities, quiet and signaling NaNs, signed zeros, and signed minimum
subnormals.

## `NRPyLeakage_compute_neutrino_luminosities.c`

Public routine: `NRPyLeakage_compute_neutrino_luminosities`.

Flow:

1. Return `ghl_error_used_disabled_hdf5` when `GHL_DISABLE_HDF5` is set.
2. Query tabulated EOS composition and chemical potentials with
   `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`; return any EOS
   error directly.
3. Convert `rho` to cgs units through `NRPyLeakage_units_geom_to_cgs_D`, then
   derive scattering and transition populations and the kinetic degeneracy
   difference from `rho_cgs`, `T`, `X_n`, and `X_p`. Return helper errors before
   output writes and normalize only accepted endpoint roundoff for all later
   uses of the fractions.
4. Form `q = muhat - T*(eta_n-eta_p)` and evaluate the paired shifted beta
   emission and absorption moments algebraically when both free populations
   are positive. Zero populations retain zero charged-current moments.
5. Run the remaining source-owned generated formula blocks for emissivity,
   opacity-like denominators, Fermi-Dirac factors, and optical-depth
   suppression. Calls through `NRPYLEAKAGE_FD_OR_RETURN` propagate invalid
   Fermi-Dirac keys.
6. Use metric/lapse/Lorentz input in the luminosity scaling:
   `alpha`, the six spatial metric components, and `W` enter the
   `NRPyLeakage_units_cgs_to_geom_Q` luminosity prefactor before writeback;
   `NRPyLeakage_units_geom_to_cgs_D` is the earlier density conversion.
7. Write `lum->nue`, `lum->anue`, and `lum->nux`.

Finite handling: local `EnsureFinite` wraps selected generated subexpressions
with `robust_isfinite` fallback to a small positive value. After writeback,
`nrpyl_sanitize_luminosities` maps any non-finite luminosity output to the
neutral zero-emission value.

Nearest tests: `Unit_Tests/unit_test_nrpyleakage_luminosities.c` directly
checks selected Fermi-Dirac branches, generates luminosity fixtures, recomputes
`NRPyLeakage_compute_neutrino_luminosities`, and reads `nue`, `anue`, and `nux`
fixtures. It consumes all three local `luminosity_pert_test_fail` return values
and aborts with the row index on the first mismatch. Its generator draws each
base state once and evaluates it unperturbed and perturbed before drawing the
next row, so the two files are matched row by row.

## `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms.c`

Public routine:
`NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`.

Flow:

1. Return `ghl_error_used_disabled_hdf5` when `GHL_DISABLE_HDF5` is set.
2. Query tabulated EOS composition and chemical potentials with
   `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`; return any EOS
   error directly.
3. Convert `rho` to cgs units, then derive scattering and transition
   populations and the kinetic degeneracy difference from `rho_cgs`, `T`,
   `X_n`, and `X_p`. Return helper errors before output writes and normalize
   only accepted endpoint roundoff for all later fraction uses.
4. Form `q = muhat - T*(eta_n-eta_p)` and evaluate paired shifted beta
   emission and absorption moments algebraically when both free populations
   are positive. Zero populations retain zero charged-current moments.
5. Run the remaining source-owned generated formula blocks for rate terms,
   opacity terms, optical-depth limited source terms, and Fermi-Dirac factors.
   Calls through `NRPYLEAKAGE_FD_OR_RETURN` propagate invalid keys.
6. Write `*R_source` and `*Q_source`.
7. Write all six opacity entries: `kappa->nue[0..1]`,
   `kappa->anue[0..1]`, and `kappa->nux[0..1]`.

Finite handling: this file's `EnsureFinite` uses `robust_isfinite` from
`GRHayL/include/ghl_nrpyleakage.h` and replaces selected non-finite intermediate
terms with a small positive value. After writeback, `nrpyl_sanitize_sources`
maps either non-finite signed source to neutral zero, while
`nrpyl_sanitize_opacities` maps non-finite opacities to the established small
positive floor.

Nearest tests: `Unit_Tests/unit_test_nrpyleakage_optically_thin_gas.c` calls
this routine in its RHS, divides `R_source` and `Q_source` by `rho`, advances
`Y_e` and `eps` with RK4, and reads fixture replay. Its comparison helper also
consumes every `ghl_pert_test_fail` result and aborts with the evolution time
on the first mismatch. Every EOS temperature inversion and the generator's
initial energy lookup are checked with `ghl_abort_if_error`, so a failed lookup
cannot feed a later right-hand side. Opacity writes get execution coverage
there through the same call but are not compared.

## `NRPyLeakage_compute_neutrino_opacities.c`

Public routine: `NRPyLeakage_compute_neutrino_opacities`.

Flow:

1. Return `ghl_error_used_disabled_hdf5` when `GHL_DISABLE_HDF5` is set.
2. Query tabulated EOS composition and chemical potentials with
   `ghl_tabulated_compute_muhat_mue_mup_mun_Xn_Xp_from_T`; return any EOS
   error directly.
3. Convert `rho` to cgs units through `NRPyLeakage_units_geom_to_cgs_D`, then
   derive scattering and transition populations and the kinetic degeneracy
   difference from `rho_cgs`, `T`, `X_n`, and `X_p`. Return helper errors before
   output writes and normalize only accepted endpoint roundoff for all later
   fraction uses.
4. Form `q = muhat - T*(eta_n-eta_p)` and evaluate paired shifted beta
   absorption moments algebraically when both free populations are positive.
   Zero populations retain zero charged-current moments.
5. Run the remaining source-owned generated formula blocks for
   absorption/scattering opacity entries. Calls through
   `NRPYLEAKAGE_FD_OR_RETURN` propagate invalid keys.
6. Write all six opacity entries: `kappa->nue[0..1]`,
   `kappa->anue[0..1]`, and `kappa->nux[0..1]`.
7. Pass all written opacity entries through `nrpyl_sanitize_opacities`; any
   non-finite final value is reset to a small positive value.

Finite handling: local `EnsureFinite` handles selected generated
subexpressions, then `nrpyl_sanitize_opacities` handles non-finite output
entries.

Nearest tests: `Unit_Tests/unit_test_nrpyleakage_constant_density_sphere.c`
directly calls this routine for interior and exterior states, stores the six
opacity fields on the grid, and reads opacity/depth fixtures. Its comparison
return values are consumed, and the serial validation loop aborts with the grid
coordinates on the first mismatch. Source-term tests do not cover this
implementation: the combined
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

Nearest tests: `Unit_Tests/unit_test_nrpyleakage_luminosities.c` checks selected
valid keys. `Unit_Tests/unit_test_nrpyleakage_physics.c` checks every valid key
from `0` through `5` in both approximation branches.
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
caller passes each neighbor pair in the public minus-then-plus order. Its flat
metric makes the minimum invariant under a paired swap, so this test still does
not verify directional argument mapping for unequal plus/minus face metrics.

## Ground Truth References

- Original NRPyLeakage development repository:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD
- NRPyLeakage symbolic implementation and C-generation notebook:
  https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD/blob/master/Tutorial-Leakage_Scheme-Implementation.ipynb
