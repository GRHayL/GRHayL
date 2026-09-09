# Radiation M1 Integration Contract

This document defines the host-neutral boundary for the Radiation M1 code. The
library provides pointwise operations on one local state or one prepared face;
the caller owns mesh traversal, storage, reconstruction, staging, boundaries,
and matter evolution.

## Transport

`ghl_m1_compute_neutrino_four_point_transport_flux` is the canonical prepared
face operation. It receives:

- four neighboring states in the public component order `{N,E,Fx,Fy,Fz}`;
- undensitized physical fluxes on the two cells adjacent to the face;
- the two adjacent, uncapped wave speeds;
- face opacity, coordinate spacing, and the face metric;
- the configured limiter parameters.

The host assembles the four-point stencil and supplies the physical operands.
The operation computes the low Rusanov flux, the componentwise limiter and
sawtooth flag, the opacity suppression, and the high/low blend for all five
components. It multiplies the final face result by `sqrt_detgamma` exactly
once. `delta_x` is coordinate spacing; the separate
`ghl_m1_compute_face_normal_delta_l` helper returns proper normal spacing when
that quantity is required by the host.

The four-point operation always uses metric light-cone speeds. It does not
apply an optical-depth wavespeed cap or a separate diffusion correction. The
generic two-state Rusanov operation remains available to other shared callers,
but is not an alternative neutrino M1 transport path. The host performs the
final flux divergence.

When electron pair channels are active, the host includes their
partner-dependent absorption in the prepared opacity; scalar provider
`kappa_tr` contains only independent absorption and scattering. See
[the pair transport-opacity contract](PAIR_SOURCE_MODEL.md#opacity-supplied-to-transport).

## Fixed M1 method

The neutrino M1 implementation has one primary numerical path. Realizability
repair uses the linear `E^2(1-epsilon)/F^2` rescale. Closure is fully
four-dimensional. Every published pressure tensor satisfies
`gamma_ij P^ij = E`, symmetry, and positive semidefiniteness within the existing
floating-point validation tolerances. Finite candidates that fail the physical
PSD check, or exact-zero-flux candidates that fail tensor validation, use the
built-in Eulerian Minerbo admissibility fallback; that result is marked
`four_point_compatibility=false` and `solve_status=endpoint_fallback`. The
fallback is not user-selectable, and nonfinite, invalid-root, and failed
fallback states remain errors.
The local implicit source update uses the existing finite-difference Jacobian
and Newton solver, including its line search and fallback substepping.
At an exact-zero-flux start, an admissible source predictor may seed Newton
away from the closure's directional discontinuity. It changes only the initial
iterate; residuals and convergence normalization retain the original source
base. There is no small-flux cutoff or relaxed solve tolerance.

The realizability margin, energy floor, closure root and residual tolerances,
closure iteration limits, finite-difference step sizes, Newton iteration limits
and tolerances, and four-point limiter/opacity-suppression parameters remain
configurable numerical controls within their documented ranges. Legacy fields
and arguments retained for source or ABI compatibility do not select another
wave-speed, repair, closure, or diffusion method.

## Electron-flavor pair source update

Electron-flavor pair reactions require both species in one pointwise call.
The conserving operation composes the independent charged-current/scattering
update with a joint pair update, then publishes both species and their exchange
packets together. Its grey collision equations, first-order splitting,
baryon-number normalization, and rejection behavior are defined in
[PAIR_SOURCE_MODEL.md](PAIR_SOURCE_MODEL.md). Single-species source operations
reject electron-flavor pair fields and unsupported aggregate non-charged-current
number rates. Hosts using the default production provider must select the
paired operation for the electron flavors; the lumped heavy species retains
its single-species source path.

## Single-species source update

`ghl_m1_solve_neutrino_source_update` is a transactional, pointwise source
dispatcher. It receives frozen metric and fluid primitives, a frozen rate
bundle, the pre-transport state, the transport-predicted state, the stage
timestep, and the conserved baryon-density normalization.

The transport-predicted state is the source base. Returned radiation
increments are therefore measured from that state. The source update retains
its internal thin, thick-equilibrium, scattering-dominated, general implicit,
and recovery handling, but these branches do not select a different M1
transport or closure method. It does not own stage counters, matter updates,
`Con2Prim`, or rate refreshes.

`dt` is the coordinate-time stage timestep. Every source-update branch uses
the corresponding `dt_alpha = alpha * dt`, including stiffness thresholds,
proper-time predictors, explicit updates, and the implicit residual/Jacobian.

The optional `enforce_mean_energy_bounds` control is a final-endpoint check. It
is applied after state repair to the endpoint produced by the ordinary implicit
route, each thin/thick/scattering compatibility branch, and both thin-update
wrappers. Each positive bound is checked independently against the comoving
endpoint ratio `J*Gamma_N/N`; a nonpositive lower or upper bound disables that
bound. The bounds reject an invalid endpoint and never clamp it; when `N == 0`,
the existing ratio-check skip is retained.

For the branched source policy, `thermalized_number_threshold < 0` disables the
optional equilibrium mean-energy number projection. A threshold of zero selects
that projection whenever `dt_alpha*kappa_a_N >= 0`, including zero opacity and
`dt == 0`, so it can change `N` even without a time-integrated number source.
The projection is separate from both backward-Euler number integration and the
subsequent `N_floor` repair. The default policy keeps the threshold negative.

On failure, the caller's output state remains the transport state and the
exchange packet is zero. A successful packet contains the undensitized
radiation increments, already-densitized equal-and-opposite matter energy and
momentum increments, the total radiation number increment, and the separate
charged-current lepton increment. The canonical source update uses only the
charged-current increment for the `Y_e` recommendation; total radiation-number
exchange remains a separate diagnostic.

## Coupled limiter

The host solves all three species into temporary states and packets before
publishing anything. It aggregates the packets, intersects radiation,
matter, `Y_e`, EOS, and `Con2Prim` admissibility intervals, and applies one
scalar `theta` to every radiation, matter, and `Y_e` increment. Matter recovery
is performed on a temporary candidate. A provider, solve, limiter, or recovery
failure leaves all caller outputs unchanged. Rates are refreshed only at the
next stage declared by the host.

No host driver is included in this library. A downstream host must implement
the coupled limiter and publication sequence described above; that host work
is outside the GRHayL build and dependency boundary.

## Ownership rules

`GRHayL/Radiation` does not include host-framework headers or own grid arrays,
ghost zones, schedules, AMR operations, boundary conditions, matter variables,
or stage sequencing. The host owns those responsibilities and supplies the
canonical transport operands and stage inputs. Downstream consumers may
implement that boundary in their own projects; GRHayL remains host agnostic.

## Current evidence boundary

The focused invariant tests are built and run with
`python3 scripts/test_radiation.py --hdf5 disabled --mode debug` (or
`--hdf5 enabled --mode opt`). The runner follows the Radiation source manifests
and compiles only the required Core, Flux_Source, and EOS dependencies. These
regressions cover the zero-flux closure and coupled pair-source contracts. The
production library boundary remains defined by the installed header
`GRHayL/include/ghl_m1.h`, the active manifests
`GRHayL/Radiation/make.code.defn` and
`GRHayL/Radiation/Neutrinos/make.code.defn`, and the source files named by
those manifests. This library-level boundary does not establish a complete
mesh evolution, framework integration, or physical validation.
