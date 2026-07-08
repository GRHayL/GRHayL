# Induction HLL Flux Contract

## Purpose

This page routes callers packing `ghl_HLL_vars` for the Induction magnetic HLL
flux helpers. It is a caller contract for `A_x`, `A_y`, and `A_z` flux
contributions, not a replacement for the source, tests, or Doxygen evidence.

## Public API

`ghl_HLL_vars` lives in
[GRHayL/include/ghl_induction.h](../../../GRHayL/include/ghl_induction.h). It
stores two magnetic-field components, two velocity components reconstructed in
two directions, and characteristic speeds for those two directions. The labels
`1` and `2` are cross-product component labels for the selected `A_i`; they
are not fixed coordinate indices.

Component packing:

| RHS component | `A_dir` in helpers | label `1` packs | label `2` packs |
| --- | ---: | --- | --- |
| `A_x` | `1` | `y` component | `z` component |
| `A_y` | `2` | `z` component | `x` component |
| `A_z` | `3` | `x` component | `y` component |

`ghl_HLL_flux_with_B` is implemented in
[GRHayL/Induction/HLL_flux_with_B.c](../../../GRHayL/Induction/HLL_flux_with_B.c).
It takes undensitized staggered `B` input through `ghl_HLL_vars` and a
caller-provided `psi6`. The function multiplies the HLL result by `psi6`;
callers must supply `psi6` at the selected `A_i` grid point.

`ghl_HLL_flux_with_Btilde` is implemented in
[GRHayL/Induction/HLL_flux_with_Btilde.c](../../../GRHayL/Induction/HLL_flux_with_Btilde.c).
It expects densitized magnetic input, `Btilde`, already packed into the same
`ghl_HLL_vars` fields. Do not pass undensitized `B` to this path unless the
caller has already applied consistent densitization.

Both source files are compiled by
[GRHayL/Induction/make.code.defn](../../../GRHayL/Induction/make.code.defn).

## Packing Audit

The test helpers show the concrete array-index contract:
[Unit_Tests/compute_A_flux_with_B.c](../../../Unit_Tests/compute_A_flux_with_B.c)
and
[Unit_Tests/compute_A_flux_with_Btilde.c](../../../Unit_Tests/compute_A_flux_with_Btilde.c).

For each `A_dir`, the helpers compute:

| RHS component | `v1*`, `v2*` sample offset | `B1*` sample offset | `B2*` sample offset | speed fields |
| --- | --- | --- | --- | --- |
| `A_x` | `(i, j+1, k+1)` | `(i, j, k+1)` | `(i, j+1, k)` | `c1_*` from label `1` at `B2` offset; `c2_*` from label `2` at `B1` offset |
| `A_y` | `(i+1, j, k+1)` | `(i+1, j, k)` | `(i, j, k+1)` | `c1_*` from label `1` at `B2` offset; `c2_*` from label `2` at `B1` offset |
| `A_z` | `(i+1, j+1, k)` | `(i, j+1, k)` | `(i+1, j, k)` | `c1_*` from label `1` at `B2` offset; `c2_*` from label `2` at `B1` offset |

The direction-dependent offset logic is in the helper-local `v_offset`,
`B1_offset`, and `B2_offset` arrays. The public flux functions receive only
the packed `ghl_HLL_vars`; they do not know `A_dir`, grid layout, or whether
the caller's arrays came from reconstruction, interpolation, or stored staggered
variables.

`compute_A_flux_with_B.c` also computes caller-side `psi6` from `phi_bssn`
before calling `ghl_HLL_flux_with_B`. The `Btilde` helper has no `psi6` path.

## Caller Hazards

- Swapped `B1/B2` changes the cross-product orientation.
- Wrong left/right reconstruction suffixes change the four velocity states:
  `rr`, `rl`, `lr`, and `ll`.
- Wrong characteristic-speed direction breaks the HLL weighting. Induction
  depends on Flux_Source characteristic-speed conventions from
  [GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h)
  and `ghl_calculate_characteristic_speed_dirn*` sources under
  [GRHayL/Flux_Source/](../../../GRHayL/Flux_Source/); keep detailed
  Flux_Source behavior routed there.
- Inconsistent densitization mixes the `B` and `Btilde` contracts. Use
  `ghl_HLL_flux_with_B` only when `psi6` correctly densitizes undensitized `B`;
  use `ghl_HLL_flux_with_Btilde` only when magnetic inputs already match the
  densitized convention.

## Test Evidence

Direct HLL coverage is in
[Unit_Tests/unit_test_HLL_flux.c](../../../Unit_Tests/unit_test_HLL_flux.c).
It runs both `ghl_HLL_flux_with_B` and `ghl_HLL_flux_with_Btilde` through the
shared helper paths and checks `Ax_rhs`, `Ay_rhs`, and `Az_rhs`.

ET Legacy HLL coverage is in
[Unit_Tests/unit_test_ET_Legacy_HLL_flux.c](../../../Unit_Tests/unit_test_ET_Legacy_HLL_flux.c).
It exercises the `B` path against ET Legacy reference data.

Generators:

- [Unit_Tests/data_gen/unit_test_data_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_HLL_flux.c)
  writes `HLL_flux_input.bin`,
  `HLL_flux_with_B_output.bin`, `HLL_flux_with_B_output_pert.bin`,
  `HLL_flux_with_Btilde_output.bin`, and
  `HLL_flux_with_Btilde_output_pert.bin`.
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c)
  writes `ET_Legacy_HLL_flux_input.bin` and
  `ET_Legacy_HLL_flux_input_pert.bin`; the test reads
  `ET_Legacy_HLL_flux_output.bin` and
  `ET_Legacy_HLL_flux_output_pert.bin`.

## Evidence Links

- [AGENTS.md](../../../AGENTS.md)
- [wiki/gems/induction.md](../induction.md)
- [GRHayL/include/ghl_induction.h](../../../GRHayL/include/ghl_induction.h)
- [GRHayL/Induction/HLL_flux_with_B.c](../../../GRHayL/Induction/HLL_flux_with_B.c)
- [GRHayL/Induction/HLL_flux_with_Btilde.c](../../../GRHayL/Induction/HLL_flux_with_Btilde.c)
- [GRHayL/Induction/make.code.defn](../../../GRHayL/Induction/make.code.defn)
- [Unit_Tests/unit_test_HLL_flux.c](../../../Unit_Tests/unit_test_HLL_flux.c)
- [Unit_Tests/unit_test_ET_Legacy_HLL_flux.c](../../../Unit_Tests/unit_test_ET_Legacy_HLL_flux.c)
- [Unit_Tests/compute_A_flux_with_B.c](../../../Unit_Tests/compute_A_flux_with_B.c)
- [Unit_Tests/compute_A_flux_with_Btilde.c](../../../Unit_Tests/compute_A_flux_with_Btilde.c)
- [Unit_Tests/data_gen/unit_test_data_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_HLL_flux.c)
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_HLL_flux.c)
- [docs/raw/Induction.dox](../../../docs/raw/Induction.dox)
