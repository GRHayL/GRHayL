# Reconstruction Face And Stencil Contract

## Routing Purpose

Use this page when checking Reconstruction face orientation, stencil width,
wrapper/helper boundaries, or caller impact. Start from the hub
[`wiki/gems/reconstruction.md`](../reconstruction.md), then use the public
declarations in
[`GRHayL/include/ghl_reconstruction.h`](../../../GRHayL/include/ghl_reconstruction.h)
and the method sources under
[`GRHayL/Reconstruction/`](../../../GRHayL/Reconstruction/) as the ground truth.

## Public Left-Face Contract

The public reconstruction wrappers return the right and left states on the left
face of the current cell. For cell `i`, that face is `i-1/2`:

- `Ur` is the value on the right side of that left face, from the current-cell
  side: `U(i-1/2 + epsilon)`.
- `Ul` is the value on the left side of that left face, from the neighbor-cell
  side: `U(i-1/2 - epsilon)`.

This contract is stated in the Reconstruction Doxygen source
[`docs/raw/Reconstruction.dox`](../../../docs/raw/Reconstruction.dox) and in
the method comments for PLM, PPM, and WENOZ wrappers. It matters because
Flux_Source and Induction consume reconstructed right/left face states rather
than interpolated cell-center values.

## Stencil Shapes

### PLM 4-Point Inputs

The PLM public wrappers
`ghl_minmod_reconstruction`, `ghl_mc_reconstruction`, and
`ghl_superbee_reconstruction` take `const double U[4]`. For the public left-face
contract at `i-1/2`, the stencil is `i-2` through `i+1`:

| Array slot | Cell |
| --- | --- |
| `U[0]` | `i-2` |
| `U[1]` | `i-1` |
| `U[2]` | `i` |
| `U[3]` | `i+1` |

The PLM test and data generator pass `&var[index-2]`, matching this 4-point
shape in
[`Unit_Tests/unit_test_PLM_reconstruction.c`](../../../Unit_Tests/unit_test_PLM_reconstruction.c)
and
[`Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c`](../../../Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c).

### PPM And WENOZ 6-Point Wrappers

The public PPM wrappers take a 6-point variable stencil. For the face `i-1/2`,
comments in
[`GRHayL/Reconstruction/PPM/ppm_reconstruction.c`](../../../GRHayL/Reconstruction/PPM/ppm_reconstruction.c)
and
[`GRHayL/Reconstruction/PPM/ppm_reconstruction_with_steepening.c`](../../../GRHayL/Reconstruction/PPM/ppm_reconstruction_with_steepening.c)
state that `var_data[6]` spans `i-3` through `i+2`. The steepened wrapper also
takes `pressure[6]`, `Gamma_eff`, `ftilde[2]`, and `ghl_parameters`.

`ghl_compute_ftilde` uses the same 6-point pressure and flux-direction velocity
stencils and fills `ftilde[2]` through two 5-point shock-detection calls. Keep
`ftilde[2]`, pressure, velocity, and reconstructed variable stencils aligned.

The WENOZ public wrapper
[`GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c`](../../../GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c)
takes `const double U[6]`. Its tests pass `&var[index-3]`, so the wrapper input
has the same 6-point `i-3` through `i+2` shape used by PPM:
[`Unit_Tests/unit_test_WENOZ_reconstruction.c`](../../../Unit_Tests/unit_test_WENOZ_reconstruction.c)
and
[`Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c`](../../../Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c).

### 5-Point Helper Routines

PPM and WENOZ both adapt 6-point public inputs through 5-point helper calls:

- `ghl_ppm_compute_for_cell(ftilde, U[5], Ur_ptr, Ul_ptr)` computes one
  cell-centered PPM right/left face pair.
- `ghl_ppm_compute_for_cell_with_steepening(params, pressure[5], Gamma_eff,
  ftilde, U[5], Ur_ptr, Ul_ptr)` adds steepening before flattening and
  monotonization.
- `ghl_shock_detection_ftilde(params, P[5], v_flux_dirn[5])`,
  `ghl_slope_limit`, and `ghl_steepen_var` support the PPM helpers.
- `ghl_wenoz_reconstruction_right_left_faces(U[5], Ur, Ul)` computes a
  centered WENOZ right/left face pair; the 6-point wrapper calls it on adjacent
  windows and filters the outputs into the public left-face `Ur`/`Ul` contract.

For 5-point helper stencils, use `enum reconstruction_stencil` from
[`GRHayL/include/ghl_reconstruction.h`](../../../GRHayL/include/ghl_reconstruction.h):

| Enum | Value | Relative cell |
| --- | ---: | --- |
| `MINUS2` | `0` | center `-2` |
| `MINUS1` | `1` | center `-1` |
| `PLUS_0` | `2` | center |
| `PLUS_1` | `3` | center `+1` |
| `PLUS_2` | `4` | center `+2` |

## Public Declarations And Internal Helpers

[`GRHayL/include/ghl_reconstruction.h`](../../../GRHayL/include/ghl_reconstruction.h)
declares both user-facing wrappers and lower-level support routines. Treat the
PLM wrappers, PPM wrappers, WENOZ wrapper, and `ghl_compute_ftilde` as the normal
public call surface. The Doxygen/source comments mark `ghl_minmod`,
`ghl_maxmod`, `ghl_shock_detection_ftilde`, `ghl_slope_limit`, and
`ghl_steepen_var` as reconstruction-internal support. The 5-point PPM and WENOZ
helpers are declared, but they expose centered helper behavior rather than the
standard public left-face wrapper contract.

## Caller Impact

- Flux_Source: HLLE flux routines operate on right/left states at a left face,
  so any `Ur`/`Ul` orientation change must be coordinated with
  [`GRHayL/Flux_Source/`](../../../GRHayL/Flux_Source/) and flux tests. The
  source map records reconstructed face primitives as a Flux_Source dependency:
  [`wiki/source-map.md`](../../source-map.md).
- Induction: HLL magnetic flux support also uses reconstructed right/left edge
  states and reconstructs velocities in direction-dependent ways. Orientation
  or stencil changes can break the assumptions documented in
  [`GRHayL/Induction/HLL_flux_with_Btilde.c`](../../../GRHayL/Induction/HLL_flux_with_Btilde.c).
- Tests: PLM, WENOZ, and ET Legacy reconstruction tests encode stencil offsets
  and orientation. PPM paths are exercised through
  [`Unit_Tests/unit_test_ET_Legacy_reconstruction.c`](../../../Unit_Tests/unit_test_ET_Legacy_reconstruction.c),
  which builds 6-point PPM stencils from `i-3` through `i+2`.
- Public header users: signature, array-width, or output-orientation changes in
  `ghl_reconstruction.h` are public API changes for downstream users and for
  `GRHayL/include/ghl_unit_tests.h`.
- Stencil width: changing PLM 4-point, PPM 6-point, WENOZ 6-point, or helper
  5-point expectations changes ghost-zone and caller-packing requirements.
- Build lists: new or renamed source files must stay in the method-local
  `make.code.defn` files and in the parent
  [`GRHayL/Reconstruction/make.code.defn`](../../../GRHayL/Reconstruction/make.code.defn)
  subdirectory list.
- PPM parameter coupling: flattening and steepening depend on
  `ppm_flattening_epsilon`, `ppm_flattening_omega1`,
  `ppm_flattening_omega2`, `ppm_shock_k0`, `ppm_shock_eta1`,
  `ppm_shock_eta2`, and `ppm_shock_epsilon` in
  [`GRHayL/include/ghl.h`](../../../GRHayL/include/ghl.h). Parameter changes
  affect `ghl_compute_ftilde`, `ghl_shock_detection_ftilde`,
  `ghl_steepen_var`, ET Legacy reconstruction coverage, and downstream
  initialization paths noted in [`wiki/change-impact.md`](../../change-impact.md).

## Repo-Local References

- [`docs/raw/Reconstruction.dox`](../../../docs/raw/Reconstruction.dox)
- [`GRHayL/include/ghl_reconstruction.h`](../../../GRHayL/include/ghl_reconstruction.h)
- [`GRHayL/include/ghl.h`](../../../GRHayL/include/ghl.h)
- [`GRHayL/Reconstruction/PLM/`](../../../GRHayL/Reconstruction/PLM/)
- [`GRHayL/Reconstruction/PPM/`](../../../GRHayL/Reconstruction/PPM/)
- [`GRHayL/Reconstruction/WENOZ/`](../../../GRHayL/Reconstruction/WENOZ/)
- [`wiki/gems/reconstruction.md`](../reconstruction.md)
- [`wiki/source-map.md`](../../source-map.md)
- [`wiki/change-impact.md`](../../change-impact.md)
