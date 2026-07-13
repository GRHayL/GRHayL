# WENOZ Reconstruction Contract

## Routing Purpose

Use this page when working on WENO-z reconstruction call orientation, stencil
width, helper/wrapper behavior, WENOZ fixtures, or build/test routing. For
derivation and implementation details, use source as authority.

## Public Calls

- `ghl_wenoz_reconstruction`: public 6-point wrapper declared in
  [`GRHayL/include/ghl_reconstruction.h`](../../../GRHayL/include/ghl_reconstruction.h)
  and implemented in
  [`GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c`](../../../GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c).
- `ghl_wenoz_reconstruction_right_left_faces`: 5-point helper declared in
  [`GRHayL/include/ghl_reconstruction.h`](../../../GRHayL/include/ghl_reconstruction.h)
  and implemented in
  [`GRHayL/Reconstruction/WENOZ/wenoz_reconstruction_right_left_faces.c`](../../../GRHayL/Reconstruction/WENOZ/wenoz_reconstruction_right_left_faces.c).

`ghl_wenoz_reconstruction` follows the common Reconstruction contract: return
right/left values for the left face of the current cell. The WENOZ helper is
different: `ghl_wenoz_reconstruction_right_left_faces` returns centered
right/left face values from one 5-point stencil. The wrapper adapts helper
output to the common left-face contract.

## Stencil And Face Contract

`ghl_wenoz_reconstruction` accepts `U[6]`. It calls the 5-point helper twice:

- first on `&U[0]`, then assigns helper right-face output into wrapper `Ul`
- then on `&U[1]`, then assigns helper left-face output into wrapper `Ur`

This 6-point wrapper built from two 5-point helper calls lets callers use one
public left-face interface while the helper computes both faces around the
center of its local stencil.

## Algorithm Notes

The helper computes WENO-z candidate blends, then blends them with a
monotonized-central fallback near nonsmooth data. Keep KB math high level; do
not copy source formulas here. Source comments and code own detailed
provenance, license, constants, and arithmetic.

Source provenance and license details are in
[`GRHayL/Reconstruction/WENOZ/wenoz_reconstruction_right_left_faces.c`](../../../GRHayL/Reconstruction/WENOZ/wenoz_reconstruction_right_left_faces.c).
Route readers there instead of duplicating the license block in KB text.

There is one built WENOZ variant and no caller-set WENO parameter object. The
helper fixes linear weights to `{0.1, 0.6, 0.3}`, regularization epsilon to
`1e-100`, and MC fallback coefficient to `2.0` in source. Any change to these
is an implementation change, not runtime configuration.

`docs/raw/Reconstruction.dox` defines a WENO group and says only WENO-z is
supported, but an earlier sentence says current categories are only PLM and
PPM. Header, manifest, two source files, direct test, and workflows resolve
that wording conflict in favor of current WENOZ build/test membership.

## Build And Coverage

WENOZ source membership is listed in
[`GRHayL/Reconstruction/WENOZ/make.code.defn`](../../../GRHayL/Reconstruction/WENOZ/make.code.defn).

Direct WENOZ unit coverage is in
[`Unit_Tests/unit_test_WENOZ_reconstruction.c`](../../../Unit_Tests/unit_test_WENOZ_reconstruction.c).
It reads `WENOZ_reconstruction_input.bin`,
`WENOZ_reconstruction_output.bin`, and
`WENOZ_reconstruction_output_pert.bin`, then checks wrapper `Ur` and `Ul`.

WENOZ fixture generation is in
[`Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c`](../../../Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c).
The generator writes `WENOZ_reconstruction` input, output, and perturbed output
fixtures using `ghl_wenoz_reconstruction`.

The generator initializes wrapper outputs only for indices `3` through
`arraylength-3`, then writes full allocated arrays. Replay reads full arrays but
compares only that valid interior range. Boundary fixture values are therefore
indeterminate and untested. No test calls the 5-point helper directly.

CI nuance: standalone reconstruction workflow jobs run `WENOZ_reconstruction`
through the matrix in
[`github-actions-Ubuntu-gcc.yml`](../../../.github/workflows/github-actions-Ubuntu-gcc.yml)
and peer workflow files. The aggregate
[`run_tests.sh`](../../../.github/run_tests.sh) script runs PLM reconstruction
but does not run WENOZ.

## Repo-Local References

- [`docs/raw/Reconstruction.dox`](../../../docs/raw/Reconstruction.dox)
- [`GRHayL/include/ghl_reconstruction.h`](../../../GRHayL/include/ghl_reconstruction.h)
- [`GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c`](../../../GRHayL/Reconstruction/WENOZ/wenoz_reconstruction.c)
- [`GRHayL/Reconstruction/WENOZ/wenoz_reconstruction_right_left_faces.c`](../../../GRHayL/Reconstruction/WENOZ/wenoz_reconstruction_right_left_faces.c)
- [`GRHayL/Reconstruction/WENOZ/make.code.defn`](../../../GRHayL/Reconstruction/WENOZ/make.code.defn)
- [`Unit_Tests/unit_test_WENOZ_reconstruction.c`](../../../Unit_Tests/unit_test_WENOZ_reconstruction.c)
- [`Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c`](../../../Unit_Tests/data_gen/unit_test_data_WENOZ_reconstruction.c)
- [`wiki/gems/reconstruction.md`](../reconstruction.md)
- [`wiki/test-map.md`](../../test-map.md)
