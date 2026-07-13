# PLM Limiters

## Routing Purpose

Use this page to route Reconstruction PLM questions about minmod, monotonized-central (MC), superbee, shared limiter helpers, PLM build membership, and direct PLM test coverage. Keep detailed derivations in source comments and Doxygen source; repo-local C source and public headers remain ground truth.

## Public PLM Entry Points

Declarations live in [`GRHayL/include/ghl_reconstruction.h`](../../../GRHayL/include/ghl_reconstruction.h). The public PLM reconstruction routines are:

- `ghl_minmod_reconstruction`: implemented in [`GRHayL/Reconstruction/PLM/minmod_reconstruction.c`](../../../GRHayL/Reconstruction/PLM/minmod_reconstruction.c). Uses `ghl_minmod` to choose the limited slope on each side of the reconstructed face.
- `ghl_mc_reconstruction`: implemented in [`GRHayL/Reconstruction/PLM/mc_reconstruction.c`](../../../GRHayL/Reconstruction/PLM/mc_reconstruction.c). Uses repeated `ghl_minmod` calls for the monotonized-central slope choice.
- `ghl_superbee_reconstruction`: implemented in [`GRHayL/Reconstruction/PLM/superbee_reconstruction.c`](../../../GRHayL/Reconstruction/PLM/superbee_reconstruction.c). Combines `ghl_minmod` and `ghl_maxmod` to select the superbee slope.

Each routine takes `const double U[4]` plus `Ur` and `Ul` output pointers. For
the left face of cell `i`, the 4-point stencil is centered around that face and
corresponds to values from `i-2` through `i+1`. Use
[`wiki/gems/reconstruction/face-and-stencil-contract.md`](face-and-stencil-contract.md)
for the shared stencil-width and left-face `Ur`/`Ul` orientation contract.

## Limiter Helpers

The helper routines are shared slope-selection utilities, not reconstruction entry points:

- `ghl_minmod`: implemented in [`GRHayL/Reconstruction/PLM/minmod.c`](../../../GRHayL/Reconstruction/PLM/minmod.c). Used by minmod, MC, and superbee reconstruction.
- `ghl_maxmod`: implemented in [`GRHayL/Reconstruction/PLM/maxmod.c`](../../../GRHayL/Reconstruction/PLM/maxmod.c). Used by superbee reconstruction.

Changing either helper can affect multiple PLM methods. Treat helper edits as shared-behavior changes and verify every direct PLM method, not only the method that motivated the edit.

PLM has no runtime parameter struct. Limiter coefficients are fixed by each
source: minmod compares adjacent one-sided slopes, MC includes centered and
twice-one-sided candidates, and superbee combines minmod candidates with
maxmod. Exact arithmetic and tie behavior stay in the five PLM sources.

## Build Membership

PLM source membership is listed in [`GRHayL/Reconstruction/PLM/make.code.defn`](../../../GRHayL/Reconstruction/PLM/make.code.defn). The build list includes:

- [`GRHayL/Reconstruction/PLM/minmod.c`](../../../GRHayL/Reconstruction/PLM/minmod.c)
- [`GRHayL/Reconstruction/PLM/maxmod.c`](../../../GRHayL/Reconstruction/PLM/maxmod.c)
- [`GRHayL/Reconstruction/PLM/minmod_reconstruction.c`](../../../GRHayL/Reconstruction/PLM/minmod_reconstruction.c)
- [`GRHayL/Reconstruction/PLM/mc_reconstruction.c`](../../../GRHayL/Reconstruction/PLM/mc_reconstruction.c)
- [`GRHayL/Reconstruction/PLM/superbee_reconstruction.c`](../../../GRHayL/Reconstruction/PLM/superbee_reconstruction.c)

## Test Evidence

Direct PLM evidence lives in
[`Unit_Tests/unit_test_PLM_reconstruction.c`](../../../Unit_Tests/unit_test_PLM_reconstruction.c).
The test dispatches over `ghl_minmod_reconstruction`,
`ghl_mc_reconstruction`, and `ghl_superbee_reconstruction`, then compares
right/left face outputs against trusted and perturbed PLM fixture data.

Fixture generation evidence lives in
[`Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c`](../../../Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c).
The generator writes `PLM_reconstruction_input.bin`,
`PLM_reconstruction_output.bin`, and `PLM_reconstruction_output_pert.bin` for
the three PLM methods.

The broader test routing table is [`wiki/test-map.md`](../../test-map.md),
which records PLM reconstruction as direct Reconstruction/PLM coverage and
notes that this PLM test does not cover PPM.

Coverage limits:

- Replay checks only `index=2` through `arraylength-2`, the valid four-point
  stencil range, and branches to `ghl_error` when a comparison fails.
- The data generator instead calls `&var[index-2]` for every index from `0`
  through `arraylength-1`. Its first two and final calls access outside the
  allocated input. Boundary fixture entries are ignored by replay, but the
  generator itself has an out-of-bounds coverage/tooling gap.
- No direct test isolates `ghl_minmod` or `ghl_maxmod`; wrapper tests exercise
  them transitively.

## Repo-Local References

- [`docs/raw/Reconstruction.dox`](../../../docs/raw/Reconstruction.dox)
- [`GRHayL/include/ghl_reconstruction.h`](../../../GRHayL/include/ghl_reconstruction.h)
- [`GRHayL/Reconstruction/PLM/minmod.c`](../../../GRHayL/Reconstruction/PLM/minmod.c)
- [`GRHayL/Reconstruction/PLM/maxmod.c`](../../../GRHayL/Reconstruction/PLM/maxmod.c)
- [`GRHayL/Reconstruction/PLM/minmod_reconstruction.c`](../../../GRHayL/Reconstruction/PLM/minmod_reconstruction.c)
- [`GRHayL/Reconstruction/PLM/mc_reconstruction.c`](../../../GRHayL/Reconstruction/PLM/mc_reconstruction.c)
- [`GRHayL/Reconstruction/PLM/superbee_reconstruction.c`](../../../GRHayL/Reconstruction/PLM/superbee_reconstruction.c)
- [`GRHayL/Reconstruction/PLM/make.code.defn`](../../../GRHayL/Reconstruction/PLM/make.code.defn)
- [`Unit_Tests/unit_test_PLM_reconstruction.c`](../../../Unit_Tests/unit_test_PLM_reconstruction.c)
- [`Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c`](../../../Unit_Tests/data_gen/unit_test_data_PLM_reconstruction.c)
- [`wiki/gems/reconstruction.md`](../reconstruction.md)
- [`wiki/gems/reconstruction/face-and-stencil-contract.md`](face-and-stencil-contract.md)
- [`wiki/test-map.md`](../../test-map.md)
