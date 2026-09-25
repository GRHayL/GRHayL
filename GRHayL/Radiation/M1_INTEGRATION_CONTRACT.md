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

On the fallback path `chi` is the Eulerian Minerbo Eddington factor built from
the Eulerian flux factor, `xi` remains the comoving reduced-flux magnitude of
the published tensor, and `root_residual` is zero because no scalar
four-dimensional root was solved. Those three fields are therefore not related
by the Minerbo `chi(xi)` relation that holds for the primary construction.

The two states that reach the fallback are counted separately by
`ghl_m1_closure_counters::admissibility_fallback_zero_flux` and
`ghl_m1_closure_counters::admissibility_fallback_psd`. The zero-flux case is an
expected consequence of the covariant thin dyad vanishing identically at exact
zero Eulerian flux and carries no admissibility concern. The PSD case reflects
the \f$O(v)\f$-accurate relativistic thick tensor losing positive semidefiniteness:
its smallest eigenvalue can reach a sizable negative fraction of the tensor norm,
so the check is a genuine rejection rather than a tolerance artifact. Repository
measurement places that regime at an Eulerian flux transverse to the fluid
velocity with fluid Eulerian speed above roughly 0.5c; flux parallel or
antiparallel to the velocity does not reach it at any speed tested up to 0.9c,
and no state in the repository validation campaign reaches it. A host running
fast transverse flows should monitor `admissibility_fallback_psd` rather than
assume the primary construction is always published.
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

The dispatcher validates the shared `terminal_fallback_policy` on every route,
not only where the ordinary implicit solver inspects it, so a value other than
`no_update_all` is rejected identically by the default implicit policy and by
every opt-in branched compatibility branch.

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

For the opt-in branched source policy, `thermalized_number_threshold` applies to
the endpoint-number step on every branch: thin, thick-equilibrium,
scattering-dominated, and the general implicit fallback. A negative value
disables the optional equilibrium mean-energy projection. A nonnegative value
selects it when `dt_alpha*kappa_a_N >= thermalized_number_threshold`; zero
therefore selects it even for zero opacity or `dt == 0`, so it can change `N`
without a time-integrated number source. The configured thick/scattering
shortcut and thermalized-number policy comparisons use scaled products,
avoiding overflow or underflow solely during route selection. This protects
policy selection, not endpoint publication: a genuinely nonrepresentable final
endpoint remains an error and is rolled back.
The projection is separate from backward-Euler number integration and the
subsequent `N_floor` repair. The default policy keeps the threshold negative,
and the public homogeneous implicit API remains unchanged and continues its
ordinary backward-Euler number update.

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

The checkout ships these focused M1 unit-test sources under `Unit_Tests/`:

- `unit_test_m1_closure_fallback.c`
- `unit_test_m1_diffusion_flux.c`
- `unit_test_m1_error_handling.c`
- `unit_test_m1_fd_jacobian.c`
- `unit_test_m1_neutrino_rusanov_flux.c`
- `unit_test_m1_neutrino_seeded_invariants.c`
- `unit_test_m1_neutrino_source_update.c`
- `unit_test_m1_rate_provider.c`
- `unit_test_m1_thcm1_blended_rusanov.c`
- `unit_test_rusanov_flux.c`

`configure` discovers `Unit_Tests/unit_test_*.c` for its generated `tests`
target, subject to its HDF5 filtering, and maps discovered sources to
`test/unit_test_*` targets. Compilation alone is target-selection evidence.
The dedicated `.github/actions/run_m1/action.yml` invokes
`Unit_Tests/run_m1_tests.sh` to build and execute the ten listed tests in the
Radiation jobs of the compiler/OS workflows. The broad `.github/run_tests.sh`
invokes `make tests datagen` to compile discovered tests, then executes its
separate test list; it does not invoke the scoped M1 runner. The M1 runner
selects the rate-provider test's generated-table mode in HDF5 builds and its
available table-free checks without HDF5. See
[the M1 test guide](../../Unit_Tests/README.m1.md)
for scoped build commands and the stored-reference boundary. CI selection alone
does not establish a remote pass, measured coverage, complete mesh evolution, or
framework integration.

The production library boundary remains defined by the installed header
`GRHayL/include/ghl_m1.h`, the active manifests
`GRHayL/Radiation/make.code.defn` and
`GRHayL/Radiation/Neutrinos/make.code.defn`, and the source files named by
those manifests. These configured library-level checks do not establish a
complete mesh evolution, framework integration, physical validation, or
line/branch coverage.

The test sources provide local checks for closure and realizability, transport,
source updates, error handling, finite-difference Jacobians, and the rate
provider boundary. Those available source checks do not establish downstream
integration or external physical validation.
