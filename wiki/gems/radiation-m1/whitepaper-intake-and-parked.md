# Radiation M1 Design Notes And Parked Work

Current source, headers, manifests, tests, and the integration contract are
the authorities for this checkout. Design notes or proposed file inventories
are not part of the build unless the active manifests and installed headers
say otherwise.

## Current repository boundary

- Shared and neutrino M1 kernels remain in `GRHayL/Radiation/`.
- The public declarations are installed through `GRHayL/include/`.
- The neutrino M1 route has one primary path: metric light-cone speeds, the
  `E^2(1-epsilon)/F^2` realizability rescale, full four-dimensional closure, and
  four-point blended Rusanov transport without a separate diffusion correction.
  Finite non-PSD closure candidates use a flagged Eulerian Minerbo
  admissibility fallback; this is not a user-selectable alternative.
- The existing finite-difference/Newton source solver remains in use; any
  source-regime handling is internal rather than a user-selected alternative
  solver.
- Host loops, stage sequencing, matter recovery, and the coupled limiter stay
  outside the library.

## Parked work

Files present in the working tree but absent from the active manifests are
parked until their contracts, dependencies, and tests are complete. They must
not be described as production APIs or verification evidence merely because
they are present on disk.

The same rule applies to generated outputs, downstream integrations, and
external comparison material: none is assumed to exist unless it is present in
this repository and exercised by a checked-in command.

## Local authority

- [`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h)
- [`make.code.defn`](../../../GRHayL/Radiation/make.code.defn)
- [`Neutrinos/make.code.defn`](../../../GRHayL/Radiation/Neutrinos/make.code.defn)
- [`M1_INTEGRATION_CONTRACT.md`](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md)
- [`TRACEABILITY.md`](../../../GRHayL/Radiation/TRACEABILITY.md)
