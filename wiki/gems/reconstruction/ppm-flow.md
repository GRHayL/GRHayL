# PPM Flow

## Purpose

Route PPM reconstruction questions here when the work touches wrapper/helper
call flow, flattening, steepening, stencil alignment, or PPM parameters. Source
remains authoritative: start with the public declarations in
[ghl_reconstruction.h](../../../GRHayL/include/ghl_reconstruction.h), the PPM
sources in [GRHayL/Reconstruction/PPM](../../../GRHayL/Reconstruction/PPM/),
and the PPM fields in [ghl.h](../../../GRHayL/include/ghl.h).

## Call Flow

The public wrappers use 6-point stencils and return the left-face values used
by Reconstruction callers:

- `ghl_ppm_reconstruction` in
  [ppm_reconstruction.c](../../../GRHayL/Reconstruction/PPM/ppm_reconstruction.c)
  takes `ftilde[2]` and `var_data[6]`. It calls
  `ghl_ppm_compute_for_cell(ftilde[0], var_data, ...)`, stores the first
  helper right value in `var_datal`, then calls
  `ghl_ppm_compute_for_cell(ftilde[1], &var_data[1], ...)` and stores the
  shifted helper left value in `var_datar`.
- `ghl_ppm_reconstruction_with_steepening` in
  [ppm_reconstruction_with_steepening.c](../../../GRHayL/Reconstruction/PPM/ppm_reconstruction_with_steepening.c)
  follows the same two-call wrapper pattern, but passes `params`,
  `pressure[6]`, `Gamma_eff`, and each aligned 5-point pressure slice into
  `ghl_ppm_compute_for_cell_with_steepening`.
- `ghl_ppm_compute_for_cell` in
  [ppm_compute_for_cell.c](../../../GRHayL/Reconstruction/PPM/ppm_compute_for_cell.c)
  is the 5-point helper. It slope-limits neighboring differences, computes
  parabolic values at both faces of the helper cell, flattens them with one
  `ftilde` value, then applies monotonicity limiting. For helper center `c`,
  `Ur_ptr` receives the left state at `c+1/2` and `Ul_ptr` receives the right
  state at `c-1/2`; wrapper output names regain the one-face convention.
- `ghl_ppm_compute_for_cell_with_steepening` in
  [ppm_compute_for_cell_with_steepening.c](../../../GRHayL/Reconstruction/PPM/ppm_compute_for_cell_with_steepening.c)
  is the steepened 5-point helper. It computes the same provisional face
  values, calls `ghl_steepen_var`, then applies flattening and monotonicity
  limiting.
- `ghl_slope_limit` in
  [slope_limit.c](../../../GRHayL/Reconstruction/PPM/slope_limit.c)
  supplies the helper slopes. It returns a limited centered slope when adjacent
  differences have the same sign and returns zero otherwise.

The wrapper stencil is 6 points, documented in source comments as `(i-3)` to
`(i+2)` for the left face of cell `i`. Each helper sees a 5-point slice using
the `MINUS2`, `MINUS1`, `PLUS_0`, `PLUS_1`, and `PLUS_2` stencil indices from
[ghl_reconstruction.h](../../../GRHayL/include/ghl_reconstruction.h). Keep
wrapper and helper indexing aligned when changing caller loops.

`ppm_compute_for_cell.c` has an older comment that labels its helper `Ur` as a
left-face right state. That conflicts with its `U[c]`/`U[c+1]` arithmetic, the
steepened helper docs, and wrapper filtering. Use the implementation mapping
above.

## Flattening

`ghl_compute_ftilde` in
[ppm_compute_ftilde.c](../../../GRHayL/Reconstruction/PPM/ppm_compute_ftilde.c)
builds `ftilde[2]` from matching `pressure[6]` and `v_flux_dirn[6]` stencils.
It calls `ghl_shock_detection_ftilde` once on the first 5-point stencil and
once on the shifted 5-point stencil. Those two `ftilde[2]` entries must match
the pressure and velocity stencils used for the two wrapper helper calls:
`ftilde[0]` belongs with the first 5-point slice, and `ftilde[1]` belongs with
the shifted slice.

`ghl_shock_detection_ftilde` in
[shock_detection_ftilde.c](../../../GRHayL/Reconstruction/PPM/shock_detection_ftilde.c)
computes one flattening value from a 5-point pressure stencil and a matching
5-point velocity-in-flux-direction stencil. Larger returned values flatten the
parabolic face values more strongly toward the cell-centered value. The
flattening controls live in `ghl_parameters` as `ppm_flattening_epsilon`,
`ppm_flattening_omega1`, and `ppm_flattening_omega2`.

## Steepening

`ghl_ppm_reconstruction_with_steepening` and
`ghl_ppm_compute_for_cell_with_steepening` carry the density-specific
steepening convention. The source documents this path as normally used for
density reconstruction, while pressure and velocity in the ET Legacy test use
the non-steepened `ghl_ppm_reconstruction` path.

`ghl_steepen_var` in
[steepen_var.c](../../../GRHayL/Reconstruction/PPM/steepen_var.c) receives a
5-point variable stencil, matching 5-point pressure stencil, `Gamma_eff`, and
the provisional right/left face values. It only changes those face values when
its contact-discontinuity, second-derivative, and relative-change checks pass.
The steepening checks couple density-like input, pressure, `Gamma_eff`, and PPM
fields in `ghl_parameters`: `ppm_shock_k0`, `ppm_shock_eta1`,
`ppm_shock_eta2`, and `ppm_shock_epsilon`.

`Gamma_eff` is a caller-provided steepening input. In
[unit_test_ET_Legacy_reconstruction.c](../../../Unit_Tests/unit_test_ET_Legacy_reconstruction.c),
the test computes it from hybrid EOS data before calling the steepened density
path.

## Parameters

PPM parameter storage is in
[ghl.h](../../../GRHayL/include/ghl.h). Defaults are set by
[initialize_params.c](../../../GRHayL/GRHayL_Core/initialize_params.c):

- `ppm_flattening_epsilon`
- `ppm_flattening_omega1`
- `ppm_flattening_omega2`
- `ppm_shock_k0`
- `ppm_shock_eta1`
- `ppm_shock_eta2`
- `ppm_shock_epsilon`

`ghl_initialize_params` sets these defaults: flattening epsilon `0.33`, omega1
`0.75`, omega2 `10.0`; shock `k0` `0.1`, eta1 `20.0`, eta2 `0.05`, and epsilon
`0.01`. GRHayLib exposes matching defaults and then overwrites the initialized
fields from Cactus parameters. PPM routines do not validate these values;
caller-owned `ghl_parameters` is a precondition.

Changing these fields affects both flattening and steepening behavior. Downstream
parameter plumbing is visible in
[GRHayLib initialize_and_shutdown.c](../../../implementations/GRHayLib/src/initialize_and_shutdown.c)
and [GRHayLib param.ccl](../../../implementations/GRHayLib/param.ccl).

## Build And Tests

The PPM build list is
[GRHayL/Reconstruction/PPM/make.code.defn](../../../GRHayL/Reconstruction/PPM/make.code.defn),
included through the parent Reconstruction build list.

ET Legacy exercises PPM paths; no dedicated `unit_test_PPM_reconstruction.c` or
PPM data generator is visible. The relevant test is
[unit_test_ET_Legacy_reconstruction.c](../../../Unit_Tests/unit_test_ET_Legacy_reconstruction.c),
and its generator is
[unit_test_data_ET_Legacy_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c).
The broader test routing is in [test-map.md](../../test-map.md).

Execution status is split:

- `.github/run_tests.sh` executes ET Legacy reconstruction, which directly
  exercises both PPM wrappers and `ghl_compute_ftilde`.
- All five compiler/OS workflow files include `reconstruction` in their
  ET-Legacy matrix. Standalone Reconstruction matrices contain only PLM and
  WENOZ, not PPM.
- ET Legacy data generation writes input and perturbed-input files only. Replay
  consumes downloaded input, trusted output, and perturbed-output fixtures;
  it never consumes `ET_Legacy_reconstruction_input_pert.bin`.
- No direct test isolates flattening, slope limiting, steepening, or either
  5-point helper. Current evidence is transitive through ET Legacy fixture
  replay.

## Repo-Local References

- [docs/raw/Reconstruction.dox](../../../docs/raw/Reconstruction.dox) as
  read-only evidence for Reconstruction method intent and face convention.
- [GRHayL/Reconstruction/PPM/ppm_reconstruction.c](../../../GRHayL/Reconstruction/PPM/ppm_reconstruction.c)
- [GRHayL/Reconstruction/PPM/ppm_reconstruction_with_steepening.c](../../../GRHayL/Reconstruction/PPM/ppm_reconstruction_with_steepening.c)
- [GRHayL/Reconstruction/PPM/ppm_compute_for_cell.c](../../../GRHayL/Reconstruction/PPM/ppm_compute_for_cell.c)
- [GRHayL/Reconstruction/PPM/ppm_compute_for_cell_with_steepening.c](../../../GRHayL/Reconstruction/PPM/ppm_compute_for_cell_with_steepening.c)
- [GRHayL/Reconstruction/PPM/ppm_compute_ftilde.c](../../../GRHayL/Reconstruction/PPM/ppm_compute_ftilde.c)
- [GRHayL/Reconstruction/PPM/shock_detection_ftilde.c](../../../GRHayL/Reconstruction/PPM/shock_detection_ftilde.c)
- [GRHayL/Reconstruction/PPM/slope_limit.c](../../../GRHayL/Reconstruction/PPM/slope_limit.c)
- [GRHayL/Reconstruction/PPM/steepen_var.c](../../../GRHayL/Reconstruction/PPM/steepen_var.c)
- [GRHayL/Reconstruction/PPM/make.code.defn](../../../GRHayL/Reconstruction/PPM/make.code.defn)
- [GRHayL/include/ghl_reconstruction.h](../../../GRHayL/include/ghl_reconstruction.h)
- [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h)
- [GRHayL/GRHayL_Core/initialize_params.c](../../../GRHayL/GRHayL_Core/initialize_params.c)
- [Unit_Tests/unit_test_ET_Legacy_reconstruction.c](../../../Unit_Tests/unit_test_ET_Legacy_reconstruction.c)
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_reconstruction.c)
- [wiki/gems/reconstruction.md](../reconstruction.md)
- [wiki/test-map.md](../../test-map.md)
