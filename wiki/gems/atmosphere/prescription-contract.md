# Atmosphere Prescription Contract

## Purpose

This page records the active Atmosphere prescription contract and the known
source drift around the unused radial-falloff routine. Source, Doxygen, tests,
and build lists remain the authority; this page routes readers to those files
without copying full APIs or implementation blocks.

Read with [wiki/gems/atmosphere.md](../atmosphere.md) for the module router.

## Surface Status

| Path | Current evidence status |
| --- | --- |
| Constant | Public installed-header declaration; manifest-built definition; Doxygen entry; test calls. |
| Radial | Public installed-header declaration and incomplete source; absent from Atmosphere manifest, Doxygen, calls, and tests. |
| Internal-only | No evidence supports this label for radial: its declaration is installed, but implementation intent remains unknown. |
| Broken runtime path | Standard manifest does not build radial, so no standard radial runtime path exists; direct/manual compilation exposes the incomplete behavior below. |

## Active Prescription

`ghl_set_prims_to_constant_atm` is the only Atmosphere prescription that is
declared, built, documented, and called by repo tests.

- Public declaration:
  [GRHayL/include/ghl_atmosphere.h](../../../GRHayL/include/ghl_atmosphere.h).
- Built source:
  [GRHayL/Atmosphere/set_prims_to_constant_atm.c](../../../GRHayL/Atmosphere/set_prims_to_constant_atm.c).
- Build list:
  [GRHayL/Atmosphere/make.code.defn](../../../GRHayL/Atmosphere/make.code.defn).
- Installed-header list:
  [GRHayL/include/make.code.defn](../../../GRHayL/include/make.code.defn).
- Doxygen source:
  [docs/raw/Atmosphere.dox](../../../docs/raw/Atmosphere.dox).

The constant prescription resets thermodynamic and velocity primitive fields
from `ghl_eos_parameters`:

- `rho = eos->rho_atm`
- `press = eos->press_atm`
- `eps = eos->eps_atm`
- `entropy = eos->entropy_atm`
- `Y_e = eos->Y_e_atm`
- `temperature = eos->T_atm`
- `vU[0..2] = 0`

The primitive struct also contains `u0` and `BU[0..2]` in
[GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h). The constant Atmosphere
source does not assign either field. Callers that need a complete primitive
state must preserve magnetic fields separately and compute or refresh `u0`
through the primitive-limit path.

The routine returns `void`, dereferences both pointers, and performs no null,
EOS-family, bounds, or finiteness validation. It always reads all six
atmosphere thermodynamic/composition fields listed above, regardless of
`eos->eos_type`.

For simple and hybrid EOS setup,
[GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c)
populates `rho_atm`, `press_atm`, `eps_atm`, and `entropy_atm` but does not
populate `Y_e_atm` or `T_atm`. Because constant reset still reads those two
fields, callers must initialize them separately (or zero-initialize the full
EOS struct before setup). Zero-initialization makes the values defined as zero;
it does not establish that zero has valid physical meaning for the caller.
This is a current contract gap, not evidence that the fields are intentionally
irrelevant to simple/hybrid callers.

## Unsupported Radial Falloff

`ghl_set_prims_to_radial_falloff_atm` is declared and source-present, but
current repository evidence does not establish a built or supported
prescription.

- It is declared in
  [GRHayL/include/ghl_atmosphere.h](../../../GRHayL/include/ghl_atmosphere.h).
- Its source file exists at
  [GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c](../../../GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c).
- It is absent from
  [GRHayL/Atmosphere/make.code.defn](../../../GRHayL/Atmosphere/make.code.defn),
  so the Atmosphere gem build list does not compile it.
- It is absent from
  [docs/raw/Atmosphere.dox](../../../docs/raw/Atmosphere.dox), which says only
  one prescription is provided and links only the constant-density routine.
- Repository declaration/definition/call search finds no call or unit test for
  the radial name.

This status was reproduced by comparing the two `*.c` basenames in
`GRHayL/Atmosphere/` with `SRCS` in its manifest: only
`set_prims_to_radial_falloff_atm.c` remains in the source-minus-manifest set.
That proves standard manifest absence only; it does not decide whether the
routine should be completed, internalized, or removed.

If compiled and called directly, current source does not implement a radial
thermodynamic falloff: `r` is unused, and `rho`, `press`, and `eps` assignments
are commented out. It copies `entropy`, `Y_e`, and `temperature`, zeros
`vU[0..2]`, and leaves `BU[0..2]`, `u0`, `rho`, `press`, and `eps` at their
entry values. This makes the source incomplete for its stated reset purpose;
it is not evidence for a usable dormant algorithm. Maintainer intent remains
unknown.

## EOS Field Coupling

Atmosphere values live in `ghl_eos_parameters`, declared in
[GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h). EOS initialization fills
the atmosphere and floor fields that the active Atmosphere prescription later
copies into primitives.

- Simple EOS: `ghl_initialize_simple_eos` accepts `rho_atm` and `press_atm`,
  stores them on `ghl_eos_parameters`, then derives `eps_atm`,
  `entropy_atm`, and `tau_atm` from the one-piece ideal-fluid setup in
  [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).
  It does not initialize `Y_e_atm` or `T_atm`.
- Hybrid EOS: `ghl_initialize_hybrid_eos` accepts `rho_atm`, stores it, then
  computes `press_atm`, `eps_atm`, `entropy_atm`, and `tau_atm` through the
  hybrid cold-pressure/entropy helpers in
  [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).
  It does not initialize `Y_e_atm` or `T_atm`.
- Tabulated EOS: `ghl_initialize_tabulated_eos` accepts `rho_atm`, `Y_e_atm`,
  and `T_atm`, clamps min/max bounds against table bounds, stores the
  atmosphere values on `ghl_eos_parameters`, and computes `press_atm`,
  `eps_atm`, and `entropy_atm` from the table-backed `(rho,Y_e,T)` route in
  [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).

For wrapper-level EOS setup, the public
`ghl_initialize_*_eos_functions_and_params` functions in
[GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) route into the same
initialization source.

## Tests And Coverage

Direct constant-atmosphere coverage lives in
[Unit_Tests/unit_test_grhayl_core_test_suite.c](../../../Unit_Tests/unit_test_grhayl_core_test_suite.c).
That test poisons primitive fields, calls `ghl_set_prims_to_constant_atm`, and
checks density, pressure, epsilon, entropy, and zero velocity for simple and
hybrid EOS cases. Its `Y_e` and temperature assertions are inside an
`eos_type == 2` branch, but the loop stops before that value; those assignments
are not directly checked. The simple/hybrid EOS locals are not zero-initialized,
so their setup does not establish the two values before constant reset reads
them.

Do not claim direct tabulated Atmosphere test coverage from that file. Its
tabulated setup branch is commented out, includes a TODO about adding
table-based default checks, and the Atmosphere loop currently stops before the
tabulated case.

[Unit_Tests/unit_test_tabulated_eos.c](../../../Unit_Tests/unit_test_tabulated_eos.c)
covers tabulated EOS initialization and table behavior, not a direct Atmosphere
reset call.

[Unit_Tests/unit_test_ET_Legacy_primitives.c](../../../Unit_Tests/unit_test_ET_Legacy_primitives.c)
calls constant reset in its negative-conservative-density branch. The source
contains a TODO to validate that reset explicitly, so treat this as indirect
legacy-path coverage rather than a focused Atmosphere contract test.

No direct radial test or call is present in repo-local test evidence.

## Impact Links Only

Primitive-limit and recovery pages are downstream impact routes, not extra
Atmosphere prescription contracts:

- [wiki/gems/con2prim/limits-and-conversions.md](../con2prim/limits-and-conversions.md)
  documents `ghl_enforce_primitive_limits_and_compute_u0`, including the final
  `u0` computation path and EOS-bound enforcement after primitives already
  exist.
- [wiki/gems/con2prim/recovery-flow.md](../con2prim/recovery-flow.md) routes
  recovery order and shows where primitive limiting happens after recovery.

## Source Of Truth

- [AGENTS.md](../../../AGENTS.md)
- [wiki/index.md](../../index.md)
- [wiki/catalog.md](../../catalog.md)
- [wiki/gems/atmosphere.md](../atmosphere.md)
- [GRHayL/include/ghl_atmosphere.h](../../../GRHayL/include/ghl_atmosphere.h)
- [GRHayL/Atmosphere/set_prims_to_constant_atm.c](../../../GRHayL/Atmosphere/set_prims_to_constant_atm.c)
- [GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c](../../../GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c)
- [GRHayL/Atmosphere/make.code.defn](../../../GRHayL/Atmosphere/make.code.defn)
- [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h)
- [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c)
- [docs/raw/Atmosphere.dox](../../../docs/raw/Atmosphere.dox)
- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
- [Unit_Tests/unit_test_tabulated_eos.c](../../../Unit_Tests/unit_test_tabulated_eos.c)
- [Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c](../../../Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c)
