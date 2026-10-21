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

The core route and public API are present. The focused M1 unit-test sources
are shipped under `Unit_Tests/`, selected by
[`.github/run_tests.sh`](../../../.github/run_tests.sh), and
executed by the ordinary test runner used by compiler/OS workflows.
These sources and runner provide library-level test routes; source and runner
selection alone do not establish a remote pass, downstream host integration,
or physical validation. The scoped tests replay retained THC_M1
discrete-operation outputs in repository-local comparisons; see
[M1 tests and fixtures](tests-and-fixtures.md). No live cross-code run,
production evolution host, or full-evolution result is established in this
checkout.
A downstream project may perform its own composed or multidimensional
validation without becoming a GRHayL build dependency.

See [`Radiation_integration_contract.md`](../../../docs/raw/Radiation_integration_contract.md)
and [`host-integration-and-downstream.md`](host-integration-and-downstream.md).
