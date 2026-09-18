# Source-update branches, pair atomicity, and rollback

The neutrino source dispatcher is a local policy boundary around one physical
source model. Its branch labels describe how a source update is evaluated;
they do not select a second transport, closure, or wave-speed method. The
public options, path values, and diagnostics are declared in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h), and the dispatcher is
[`ghl_m1_neutrino_source_update.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_source_update.c).

## Default and opt-in policies

The zero-valued default is
`ghl_m1_neutrino_source_grhayl_implicit`. It invokes the established local
implicit E/F solve, followed by ordinary backward-Euler endpoint number
evolution, repair, endpoint checks, and exchange assembly. The default is not a
thin/thick branch selector. The public
`ghl_m1_solve_neutrino_implicit_homogeneous_update` entry point retains this
same signature and behavior; the thermalized-number projection is introduced
only by the opt-in branched dispatcher.

The opt-in
`ghl_m1_neutrino_source_branched_compatibility` policy exposes the following
internal choices:

| path | selection and operation |
| --- | --- |
| thin explicit | selected when `dt_alpha*kappa_a_E < 1` and `dt_alpha*kappa_s < 1`; E/F receives an explicit interaction update, then N uses the shared endpoint-number policy |
| thick equilibrium | enabled only by a positive host threshold and selected when `dt_alpha*sqrt(kappa_a_E*kappa_tr)` exceeds it; a stiff comoving predictor is used before the shared endpoint-number policy |
| scattering dominated | enabled only by a positive host threshold and selected when `dt_alpha*kappa_s` exceeds it; it uses the same stiff predictor machinery and endpoint-number policy |
| general implicit | the fallback branch when no compatibility shortcut is selected; it uses the policy-aware implicit E/F solve and the same endpoint-number policy |

All branches use `dt_alpha = alpha*dt`, frozen matter/rates, endpoint repair,
endpoint mean-energy checks when enabled, and the same exchange assembly. A
nonpositive thick or scattering threshold disables that shortcut. The default
options use zero for both and therefore do not opt into those compatibility
branches.

The configured thick/scattering shortcut and thermalized-number policy products
are compared in a scaled representation, so finite inputs are not
misclassified solely because a direct multiplication overflows or underflows.
This selection protection does not make a final state valid: a genuinely
nonrepresentable endpoint remains an error and is rolled back.

The branched policy is fail-closed on a closure fallback unless
`allow_closure_fallback` is explicitly enabled. That choice affects whether a
candidate is accepted; it does not change the closure algorithm.

## Number projection is separate

`thermalized_number_threshold` controls an optional equilibrium mean-energy
projection in every opt-in branched endpoint: thin, thick, scattering, and the
general implicit fallback. A negative value disables it. A nonnegative value
selects the projection when `dt_alpha*kappa_a_N` reaches the threshold; zero
therefore selects it even for zero opacity or zero timestep. The policy
comparison uses scaled products to avoid overflow or underflow during
selection. It is distinct from backward-Euler number evolution and from the
separate `N_floor` repair. A nonrepresentable final endpoint remains an error;
the default is negative.

The optional endpoint mean-energy bound checks use the final repaired ratio

```text
mean_energy_endpoint = J * Gamma_N / N.
```

Positive lower and upper bounds reject an out-of-range endpoint; they never
clamp it. The endpoint checker rejects nonfinite `N` and `N < N_floor`; after
that floor check, `N == 0` skips the ratio calculation, while any other
accepted nonzero `N` (including `N == N_floor`) is checked against the enabled
bounds. The separate mean-energy diagnostic uses its own `N <= N_floor`
invalidity rule.

## Pair source operation

Electron-flavor pair channels cannot be evaluated correctly by a single-species
call. The paired API
[`ghl_m1_solve_neutrino_pair_source_update`](../../../GRHayL/include/ghl_m1.h)
receives `{nue, anue}` together and uses
[`PAIR_SOURCE_MODEL.md`](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md).
It first advances independent charged-current/scattering terms into private
temporary states, then applies the shared pair number reaction and paired grey
energy/flux update.

The pair number increments are equal, preserving the electron-neutrino minus
antineutrino number difference. Pair exchange contributes no `Y_e` change;
charged-current lepton exchange is retained separately from the independent
source stage. The partner-dependent inverse-pair energy opacity must be added
by the host when it prepares an electron-flavor transport face opacity.

The pair operation publishes both species and both exchange packets together.
If either species fails validation, source solve, repair, endpoint bounds, or
exchange assembly, neither candidate is published. Single-species calls reject
separated pair coefficients and unsupported aggregate electron-number rates.

## Rollback and source base

Before validation completes, the dispatcher sets

```text
state_out = state_transport
exchange  = zero.
```

`state_input` is checked as the pre-transport state, but all source candidates
are based on `state_transport`. On a hard failure, closure rejection, rate
failure, double-application request, or disallowed endpoint, the output
remains that source base and the exchange packet remains zero. A terminal
implicit fallback is reported distinctly from a hard failure, but has the same
no-update publication behavior.

Diagnostics are assembled in candidate records and published with a
successful result or the documented failure accounting. They are observability
state, not a license to publish a partial state. Hosts should therefore stage
state, radiation exchange, matter exchange, and pair outputs in temporary
storage until the operation returns success.

## Host composition rule

For three species, solve `nue` and `anue` through the paired operation when
pair channels are active, solve `nux` through the single-species operation, and
aggregate only after all local results are valid. The host then applies its
common admissibility limiter to every species and all matter/lepton increments.
The library does not own that multi-species publication or matter recovery.

The [M1 integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md)
defines the frozen-input and coupled-limiter boundary. Focused local checks
are [`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c),
[`unit_test_m1_error_handling.c`](../../../Unit_Tests/unit_test_m1_error_handling.c),
and the source fixtures documented in
[`README.m1.md`](../../../Unit_Tests/README.m1.md). These tests do not prove a
full three-species host evolution or a continuum pair-transport result.
