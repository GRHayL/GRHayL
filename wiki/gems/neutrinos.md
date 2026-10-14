# Neutrinos Gem

## Purpose

This page routes Neutrinos KB work. It summarizes where to read first for
NRPyLeakage physics, generated-code provenance, public data, implementation
flow, portability, and tests without replacing source, headers, tests, CI, or
Doxygen source.

## Read Order

1. [Physics And EOS Contract](neutrinos/physics-and-eos-contract.md) for
   species, chemical potentials, density-derived free-nucleon blocking,
   producer energy conventions, physical qualification limits, units, and the
   number and energy source-term conventions.
2. [Nucleon Blocking And EOS Conventions](neutrinos/nucleon-blocking-and-eos-conventions.md)
   for the installed approximation, detailed-balance pairing, performance
   motivation, qualification evidence, and dense-matter limitation.
3. [Generator Provenance](neutrinos/generator-provenance.md) for the original
   Python notebooks, symbolic common-subexpression temporaries, and the
   boundary between ancestral generators and current checked-in C.
4. [API And Data](neutrinos/api-and-data.md) for radiation structs, public
   `NRPyLeakage_*` declarations, constants ownership, errors, and HDF5/EOS
   dependency.
5. [CompOSE EOS Adapter How-To](neutrinos/compose-eos-adapter-how-to.md) for
   adapting CompOSE state, composition, and chemical-potential outputs to the
   NRPyLeakage EOS callback and validating that boundary.
6. [Implementation Flow](neutrinos/implementation-flow.md) for the
   [GRHayL/Neutrinos/NRPyLeakage/](../../GRHayL/Neutrinos/NRPyLeakage) source set, its writeback paths, and
   the minimal direct-compilation boundary.
7. [Tests And Fixtures](neutrinos/tests-and-fixtures.md) for unit tests,
   fixture pairs, EOS table setup, and CI downloads.

## Ground Truth

- Source: `GRHayL/Neutrinos/NRPyLeakage/`
- Installed shared blocking implementation header:
  [GRHayL/include/ghl_nrpyleakage_nucleon_blocking.h](../../GRHayL/include/ghl_nrpyleakage_nucleon_blocking.h),
  with shared rate algebra in
  [GRHayL/include/ghl_nrpyleakage_rate_helpers.h](../../GRHayL/include/ghl_nrpyleakage_rate_helpers.h)
- Installed blocking model and qualification:
  [Nucleon Blocking And EOS Conventions](neutrinos/nucleon-blocking-and-eos-conventions.md)
- Public radiation structs: [GRHayL/include/ghl_radiation.h](../../GRHayL/include/ghl_radiation.h)
- Public leakage declarations and constants: [GRHayL/include/ghl_nrpyleakage.h](../../GRHayL/include/ghl_nrpyleakage.h)
- Error codes and EOS parameter types: [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h)
- Direct tabulated EOS dependencies: [GRHayL/include/ghl_nrpyeos_tabulated.h](../../GRHayL/include/ghl_nrpyeos_tabulated.h)
  and [GRHayL/include/ghl_eos_functions.h](../../GRHayL/include/ghl_eos_functions.h)
- Tests: [Unit_Tests/nrpyleakage_main.h](../../Unit_Tests/nrpyleakage_main.h) and
  `Unit_Tests/unit_test_nrpyleakage_*.c`
- Fermi-Dirac error coverage: [Unit_Tests/unit_test_code_error.c](../../Unit_Tests/unit_test_code_error.c)
- Build and CI gates: `configure`, [.github/run_tests.sh](../../.github/run_tests.sh), and
  [.github/workflows/](../../.github/workflows)
- Ancestral derivation and generator evidence: the
  [Tabulated_EOS_IllinoisGRMHD repository](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD),
  with durable roles and conventions summarized in
  [Generator Provenance](neutrinos/generator-provenance.md) and
  [Physics And EOS Contract](neutrinos/physics-and-eos-contract.md)

If this page, a child page, or external evidence conflicts with current
repo-local source, headers, manifests, or tests, trust the repo-local ground
truth and update the KB.

## Public Surface

- `ghl_neutrino_luminosities`
- `ghl_neutrino_opacities`
- `ghl_neutrino_optical_depths`
- `NRPyLeakage_Fermi_Dirac_integrals`
- `NRPyLeakage_compute_neutrino_opacities`
- `NRPyLeakage_compute_neutrino_luminosities`
- `NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms`
- `NRPyLeakage_optical_depths_PathOfLeastResistance`

The `NRPyLeakage_*` spelling is the exported family; no `ghl_*` wrapper family
exists for these calls. Radiation container types retain the `ghl_` prefix.

The frozen M1 grey neutrino provider consumes the same NRPyLeakage rate algebra
through the installed `ghl_neutrino_rate_provider.h` boundary and
[ghl_m1_nrpyleakage_kernel.c](../../GRHayL/Radiation/Neutrinos/ghl_m1_nrpyleakage_kernel.c);
route its ownership and caching questions to the
[radiation-m1 rate-provider contract](radiation-m1/rate-provider-contract.md).

## Contract Summary

- Radiation struct order is `nue`, `anue`, `nux`; opacity/depth has two slots
  per species. Slot meanings are used differently by number and energy source
  formulas but remain unnamed in the public header, so keep `[0]`/`[1]`
  wording in caller contracts.
- Endpoint-roundoff normalization of EOS free-nucleon fractions goes through
  the installed `ghl_nrpyleakage_normalize_nucleon_fractions` in
  [GRHayL/include/ghl_nrpyleakage.h](../../GRHayL/include/ghl_nrpyleakage.h),
  shared by the EOS-dependent leakage callers and the M1 provider kernel.
- EOS-dependent routines require an HDF5-enabled build and an initialized,
  compatible chemistry callback with its required state; the built-in
  StellarCollapse backend requires a loaded table. Inputs include geometric
  density, EOS-compatible positive temperature, and valid output pointers. They return before writeback on HDF5, EOS, or generated Fermi
  errors. They also return `ghl_error_nrpyleakage_blocking` before writeback
  when density-derived nucleon blocking receives invalid state or cannot
  produce finite, bounded factors, or when shifted beta-moment evaluation
  rejects invalid or nonfinite inputs/results.
  If final output sanitization replaces a nonfinite value, it retains finite
  fallback outputs and returns `ghl_error_nrpyleakage_nonfinite_output`.
- All EOS-dependent routines obtain blocking from `rho`, `T`, `Xn`, and
  `Xp`. This avoids using producer-dependent absolute `mu_n` and `mu_p` as
  kinetic occupations. They combine the density-derived kinetic degeneracy
  difference with `muhat` to form one reaction-energy shift, then apply it to
  paired algebraic charged-current emission and absorption moments. The
  remaining qualification boundary is the accuracy of that grey approximation,
  recorded in the blocking-conventions page.
- Optical-depth update is a `void` six-neighbor stencil call with no validation
  or failure channel. Metric stencil order is minus/center/plus.
- Every implementation file matches its manifest and header declarations.
- The aggregate runner invokes each HDF5 test binary, and every compiler
  workflow family configures those commands. Command presence is not a
  historical execution result. When executed, their main fixture comparisons
  consume `ghl_pert_test_fail` results and abort on the first mismatch, so a
  completed run establishes numerical agreement with the fixtures actually
  supplied.

## Scope Notes

- Neutrinos KB pages live under [wiki/](..) only; Doxygen source under `docs/**`
  is a separate authority (currently no dedicated Neutrinos page exists there).
- Keep generated formula blocks in `GRHayL/Neutrinos/NRPyLeakage/*.c`.
- The shared nucleon-blocking algebra lives in the installed
  `static inline` implementation headers
  [GRHayL/include/ghl_nrpyleakage_nucleon_blocking.h](../../GRHayL/include/ghl_nrpyleakage_nucleon_blocking.h)
  and
  [GRHayL/include/ghl_nrpyleakage_rate_helpers.h](../../GRHayL/include/ghl_nrpyleakage_rate_helpers.h).
  One algebraic density inversion serves every EOS-dependent routine, avoiding
  duplicated corrections and the runtime cost of per-call quadrature or root
  solves. These headers are caller-compiled by `#include`, not tracked through the
  NRPyLeakage `make.code.defn` (which lists only `SRCS`), and they do not add
  separately linked library symbols beyond the installed `NRPyLeakage_*` family.

- The external notebooks are provenance, not the authority for current GRHayL
  signatures or behavior. Preserve the durable derivation and interface facts
  in child pages so routine KB use does not depend on that repository remaining
  available.
- Keep constants and unit conversions in `GRHayL/include/ghl_nrpyleakage.h`;
  do not copy them into KB tables.
- Treat HDF5-enabled tabulated EOS access as part of the direct dependency
  surface for the opacity, combined source/opacity, and luminosity routines;
  the Fermi helper and optical-depth update do not call the EOS.

## Ground Truth References

- [Original Tabulated_EOS_IllinoisGRMHD repository](https://github.com/leowerneck/Tabulated_EOS_IllinoisGRMHD)
