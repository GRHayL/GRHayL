# Flux Source Source-Term Contract

## Routing Purpose

Use this page when checking `ghl_calculate_source_terms`, caller-owned metric
derivatives, or generated source-term provenance. It records inputs, outputs,
and evidence routes; equations remain in source, Doxygen, and derivation notes.

## Public API

`ghl_calculate_source_terms` and `ghl_calculate_source_terms_checked` are declared in
[GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h)
and implemented in
[GRHayL/Flux_Source/ghl_calculate_source_terms.c](../../../GRHayL/Flux_Source/ghl_calculate_source_terms.c).
The parent build list
[GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn)
compiles the source-term kernel.

Shared structs come from [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h):
`ghl_eos_parameters`, `ghl_primitive_quantities`,
`ghl_metric_quantities`, `ghl_extrinsic_curvature`, and
`ghl_conservative_quantities`.

## Caller Contract

The caller supplies:

- `eos`: EOS parameters compatible with `ghl_compute_h_and_cs2`.
- `prims`: primitive state at the point where source terms are evaluated,
  including density, pressure, velocity, magnetic field, and `u0`.
- `metric`: ADM metric at the same point.
- `metric_derivs_x`, `metric_derivs_y`, `metric_derivs_z`: caller-computed
  spatial derivatives of lapse, shift, and spatial metric in each coordinate
  direction, packed into `ghl_metric_quantities` containers.
- `curv`: extrinsic curvature at the same point.
- `cons`: caller-owned conservative output. The kernel writes source
  contributions to `cons->SD[0..2]` and `cons->tau`.

`cons->rho`, `cons->Y_e`, and `cons->entropy` are untouched. The routine does
not zero or initialize the complete conservative struct.

The checked kernel calls `ghl_compute_h_and_cs2(eos, prims, &h, &cs2)` before
building the source terms; `cs2` is unused. Its `prims` argument is mutable:
production tabulated dispatch
limits `rho`, `Y_e`, and `temperature` to table bounds, then recomputes `press` and `eps` in
place. The checked routine returns the callback's exact error and leaves `cons`
unchanged on failure; the legacy `void` wrapper aborts on that error. Callback
mutation of `prims` is not rolled back. Callers
needing unchanged primitives must pass a copy. EOS dispatch details otherwise
belong to EOS/Core routes.

Metric derivatives are not computed inside Flux_Source. Callers own their
discretization, staggering, and consistency with `metric`, `curv`, and
`prims`. Passing cell-centered metric data with face-centered derivatives, or
mixing derivative directions, changes the physical source term even though the
function signature still matches.

## Equation And Generated-Source Routes

Source-term equations are described in read-only evidence:
[docs/raw/derivation.md](../../../docs/raw/derivation.md) and
[docs/raw/Flux_Source.dox](../../../docs/raw/Flux_Source.dox). Do not copy the
long expressions into KB pages.

Local generator provenance:

- [GRHayL/Flux_Source/GRHayL_rhs.py](../../../GRHayL/Flux_Source/GRHayL_rhs.py)
  generates source-term C into an explicit empty staging directory. See the
  [generated NRPy boundary](generated-nrpy-boundary.md) for the command.
- [GRHayL/Flux_Source/IGM_All_Source_Terms.py](../../../GRHayL/Flux_Source/IGM_All_Source_Terms.py)
  contains the `ghl_calculate_source_terms` generation path.
- [GRHayL/Flux_Source/nrpy/](../../../GRHayL/Flux_Source/nrpy/) supplies local
  NRPy support modules used by the generator scripts.

Treat checked-in Python and generated C as coupled source evidence. If source
terms change, review both the generated kernel and generator path instead of
editing only one expression copy.

## Test Evidence

- [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../../Unit_Tests/unit_test_ET_Legacy_flux_source.c)
  covers ET Legacy flux/source behavior with reference data.
- [Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_flux_source.c)
  shows fixture packing for primitives, metric, metric_derivs, extrinsic curv,
  and conservative outputs.

Evidence status: declaration, definition, and normal/no-HDF5 manifest membership
are present. `unit_test_ET_Legacy_flux_source` directly calls the routine inside
a combined flux-divergence/source replay. It installs a test-local
`compute_h_and_cs2` callback, so this is not
a production EOS-dispatch test. The ordinary runner and all compiler workflows
configure that replay. `unit_test_tabulated_eos_compose` also covers a
production tabulated source call in the Ubuntu-GCC `compose-regularized-eos`
job. `unit_test_hybrid_flux` isolates the checked source callback failure and
verifies untouched output.

## Evidence Links

- [GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h)
- [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h)
- [GRHayL/Flux_Source/ghl_calculate_source_terms.c](../../../GRHayL/Flux_Source/ghl_calculate_source_terms.c)
- [GRHayL/Flux_Source/IGM_All_Source_Terms.py](../../../GRHayL/Flux_Source/IGM_All_Source_Terms.py)
- [GRHayL/Flux_Source/GRHayL_rhs.py](../../../GRHayL/Flux_Source/GRHayL_rhs.py)
- [GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn)
- [docs/raw/Flux_Source.dox](../../../docs/raw/Flux_Source.dox)
- [docs/raw/derivation.md](../../../docs/raw/derivation.md)
