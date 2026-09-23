# Radiation M1 Maintenance and Consumer Boundary

This page records library-local follow-up for the canonical M1 implementation.
GRHayL is host agnostic; production mesh integration and complete evolution
campaigns are downstream consumer responsibilities, not pending library work.

## Library-local follow-up

- Keep the canonical metric-light-cone, four-point blended Rusanov transport
  and existing finite-difference/Newton source route unchanged.
- Exercise the four-point transport operation across the full intended stencil
  and its configured limiter/opacity-suppression ranges.
- Exercise every source branch, including failure and terminal no-update
  behavior.
- Keep the three-species limiter reference in the host boundary and verify one
  shared scalar across all coupled increments.

## Intentionally outside GRHayL

Mesh traversal, reconstruction, boundaries, AMR, schedules, RK/IMEX stage
sequencing, matter updates, EOS/Con2Prim calls, and the coupled limiter remain
host responsibilities. They must not be added to `GRHayL/Radiation`.

## Current status

The core route and public API are present. Ten focused M1 unit-test sources
are shipped under `Unit_Tests/`, selected by
[`Unit_Tests/run_m1_tests.sh`](../../../Unit_Tests/run_m1_tests.sh), and
invoked by the dedicated
[Radiation CI action](../../../.github/actions/run_m1/action.yml). The
repository's normal `.github/run_tests.sh` path does not select the scoped M1
runner.
These sources and runner provide library-level test routes; source and runner
selection alone do not establish a remote pass, downstream host integration,
or physical validation. No production evolution host is required or selected in
this checkout, and no external comparison or full-evolution result is claimed.
A downstream project may perform its own composed or multidimensional
validation without becoming a GRHayL build dependency.

See [`M1_INTEGRATION_CONTRACT.md`](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md)
and [`host-integration-and-downstream.md`](host-integration-and-downstream.md).
