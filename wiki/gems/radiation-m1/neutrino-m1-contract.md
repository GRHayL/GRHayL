# Neutrino M1 Contract

This page defines the host-neutral grey, one-group, three-species local
contract. The source of truth is `GRHayL/include/ghl_m1.h` and the
implementations named by the two Radiation make manifests.

## State and ordering

Species are ordered `{nue, anue, nux}`. Each local state is
`{N, E, Fx, Fy, Fz}` and is undensitized. `N` is the transported number
density, `E` is Eulerian radiation energy density, and `F_i` are covariant
Eulerian flux components. A host may store densitized values but must
undensitize before calling local kernels and apply the declared metric factor
once when storing returned face or matter quantities.

The number floor is separate from the energy/flux realizability repair. The
default electron-fraction convention uses only the charged-current lepton
exchange; total radiation-number exchange remains a separate diagnostic.

## Frozen inputs

Each local source call receives a frozen metric, frozen fluid primitives,
frozen rate bundle, validated input state, validated transport-predicted state,
the coordinate-time stage `dt`, and positive baryon-density normalization. The
library forms `dt_alpha = alpha * dt` and uses that lapse-weighted interval for
all local radiation-source evolution, stiffness decisions, and implicit
residuals. The library performs no EOS lookup, host-matter update, grid
traversal, boundary operation, or stage tracking.

The optional `enforce_mean_energy_bounds` setting checks the comoving endpoint
ratio `J*Gamma_N/N` after repair. This applies to the ordinary implicit update,
the thin/thick/scattering compatibility paths, and both thin-update wrappers.
Each positive bound is enforced independently; a nonpositive lower or upper
bound disables that bound. An out-of-bounds endpoint is rejected without
clamping, while `N == 0` retains the existing skip of the ratio check.

In the branched source policy, `thermalized_number_threshold < 0` disables the
equilibrium mean-energy projection for the number update. Zero selects it even
when opacity and `dt` are zero, so the projection can change `N` without a
time-integrated number source. It remains distinct from backward-Euler number
integration and the separate `N_floor` repair.

## Fixed numerical path

Every neutrino M1 call uses the same primary numerical method:

- metric light-cone wave speeds, with no optical-depth cap;
- realizability repair by the `E^2(1-epsilon)/F^2` flux rescale;
- a full four-dimensional primary closure with the radiation trace invariant
  enforced on every published tensor. A finite tensor that fails the PSD
  check, or an exact-zero-flux candidate that fails tensor validation, uses
  the built-in Eulerian Minerbo admissibility fallback and is reported with
  `four_point_compatibility=false`; and
- four-point blended Rusanov transport, with no separate diffusion correction.

The existing finite-difference Jacobian and Newton source solver remains in
use for the local implicit source update, including its line search and
fallback substepping.

The realizability margin, energy floor, closure root and residual tolerances,
closure iteration limits, finite-difference step sizes, Newton iteration limits
and tolerances, and four-point limiter/opacity-suppression parameters remain
ordinary configurable numerical controls within their documented ranges. They
tune this method; they do not select an alternative method.

## Four-point transport

The canonical face operation consumes four neighboring values for each
component, physical cell-centered fluxes on the two sides, adjacent uncapped
metric light-cone speeds, face opacity, coordinate spacing, and the face metric
factor. Inputs are undensitized and use `{N,E,Fx,Fy,Fz}`. The result is one
final face flux, densitized exactly once. The host assembles the stencil and
computes the final divergence.

The operation returns limiter/sawtooth and opacity-suppression diagnostics. A
wavespeed cap or separately enabled diffusion correction is not part of the
canonical route. The generic two-state Rusanov operation remains available to
other shared callers but is not an alternative neutrino M1 path.

## Transactional source update

Electron pair, plasmon, and bremsstrahlung channels use
`ghl_m1_solve_neutrino_pair_source_update`, with inputs ordered `{nue, anue}`.
It advances independent charged-current/scattering sources and then a shared
pair reaction. Both number increments from the pair stage are equal; its
electron-fraction increment is zero. Either both final states and exchange
packets are published or both source bases are retained with zero packets.
Number-floor injections are rejected. The
[pair collision model](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md) defines
the grey occupancy approximation and the frozen-normalization split.

Single-species source calls reject separated pair coefficients and legacy
aggregate electron-number coefficients containing non-charged-current
reactions. The following dispatcher options describe the single-species path;
they do not select a different paired collision model.

The source update is transactional and uses the transport-predicted state as its
source base. Its existing finite-difference/Newton machinery handles the local
implicit solve, with line search and fallback substepping. Any thin, thick, or
scattering-regime reductions are internal source handling; they do not expose a
second M1 transport or closure method.

`interaction_sources_already_applied` is rejected with a dedicated
double-application error so an explicit source is not applied again by the
local implicit stage. Existing closure-fallback and terminal-recovery behavior
continues to follow the source-update contract.

Every successful update is measured relative to `state_transport`, not
`state_input`. The dispatcher returns the updated state, radiation increments,
the total radiation-number increment, the charged-current lepton increment,
equal-and-opposite densitized matter energy/momentum increments, and path
diagnostics. It never owns an RK counter or publishes a Newton guess.

## Failure and fallback behavior

Outputs are initialized transactionally to `state_transport` and a zero
exchange packet. Invalid inputs, rate/provider failure, closure rejection,
Newton failure, disallowed fallback, and double application leave those values
unchanged/zero. The existing terminal no-update return is preserved and is
identified separately from a hard failure. A successful branch reports its
source-regime handling and whether a closure fallback was used.

## Exchange and units

Radiation increments in the exchange packet are undensitized. Matter energy
and momentum increments are densitized with the supplied metric determinant;
the electron-fraction recommendation is derived from charged-current lepton
exchange and baryon normalization. The total number increment is not folded
into the electron-fraction field.

The baryon normalization is the undensitized Eulerian baryon number density
`n_b_cons = W*rho/m_b`. Convert a densitized host variable before passing it;
the charged-current radiation-number increment used in `dYe` is undensitized.

The host aggregates three species before limiting. It computes one admissible
scalar interval from radiation positivity/realizability, matter energy and
momentum, electron-fraction bounds, and EOS/Con2Prim restrictions. That same
scalar multiplies every species' radiation increment and every paired matter
and electron-fraction increment before matter recovery and publication.

## Provider and stage boundary

The provider owns rate construction, caching, units, and provider diagnostics.
The host chooses when to refresh the provider and passes one frozen bundle to
each local solve. Rates are refreshed only at the declared next stage. The
host initializes caller-owned `ghl_m1_neutrino_diagnostics` and
`ghl_neutrino_rate_provider_cache` records before first use and does not share
either mutable record across concurrent calls without synchronization. The
default provider is table-free and deterministic; the tabulated NRPyLeakage
backend is selected explicitly.

## Ground truth

- `GRHayL/include/ghl_m1.h`
- `GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md`
- `GRHayL/Radiation/make.code.defn`
- `GRHayL/Radiation/Neutrinos/make.code.defn`
- `GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_source_update.c`
- `GRHayL/include/make.code.defn`
