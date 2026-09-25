# Electron-flavor pair source model

Electron-flavor pair, plasmon, and bremsstrahlung emission is a coupled
`nu_e`/`anti-nu_e` operation. The current equations live in the
[pair-source contract](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md) and
the implementation is
[`ghl_m1_neutrino_pair_source.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_pair_source.c).
This leaf explains the physical ownership and approximation; the public API
signatures remain in [`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h).

## Which rates belong to the pair operator

For each electron flavor, the provider returns process-indexed fields

```text
eta_N_pair[c], eta_E_pair[c],
c in {pair, plasmon, bremsstrahlung}.
```

The number emissivity for a process is shared by `nu_e` and `anti-nu_e`; the
energy emissivity is allowed to differ. The provider preserves both raw
moments and does not impose

```text
eta_E_pair[c] / eta_N_pair[c] = J_eq / n_eq.
```

That difference matters because number and energy are separate grey moments.
Absent channels must be represented by zero process fields in the complete
rate struct. The provider-side mapping is in
[`ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c),
and the field contract is in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h).

The scalar electron-flavor `kappa_tr` contains the independent absorption and
scattering contribution. It is not overwritten with a partner-dependent pair
opacity. For `nu_x`, the process arrays are zero in the public bundle and its
already-summed pair/plasmon/bremsstrahlung content is represented by aggregate
scalar rates; the heavy-flavor multiplicity is already applied. See the
[provider boundary](../../../GRHayL/include/ghl_neutrino_rate_provider.h)
and [rate validation](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rates.c).

## Shared number reaction

The pair operation receives both electron-flavor states in `{nue, anue}`
order. During one pair substep, let

```text
n_s = N_s / Gamma_s
```

be the comoving number density represented by the Eulerian number state and its
M1-derived number-current normalization. For process `c`, with common number
emissivity `q_c`, the grey occupancy-product reaction is

```text
R_c = q_c * (1 - n_e*n_a/(n_eq_e*n_eq_a)).
```

With `h = alpha*dt_sub`, backward Euler gives one shared Eulerian number
increment `d`:

```text
d = h * sum_c(q_c) *
    (1 - (N_e+d)*(N_a+d)
         /(Gamma_e*Gamma_a*n_eq_e*n_eq_a)).
```

The admissible root is continuous with `d = 0` as `h` tends to zero. Both
species receive exactly that same increment, so the pair operator preserves
`N_e - N_a` up to arithmetic roundoff. The implementation solves the
equivalent quadratic using the smaller population when necessary to avoid
catastrophic cancellation during near-complete annihilation. Number-floor
violations reject the candidate; the operation does not inject a floor.

These equations and the numerical stabilization are specified in
[`PAIR_SOURCE_MODEL.md`](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md#shared-number-reaction)
and implemented at
[`ghl_m1_neutrino_pair_source.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_pair_source.c).

## Energy and momentum after number

After the number reaction, the partner's updated comoving number is frozen for
the E/F solve. For species `s` with partner `p`, define

```text
K_s = sum_c(eta_E_pair_s[c]/J_eq_s) * (n_p_new/n_eq_p)
Q_s = sum_c(eta_E_pair_s[c]) - K_s*J_s.
```

The existing M1 source projection then gives

```text
S_E_s = W*Q_s + K_s*H_n_s
S_i_s = W*Q_s*V_i - K_s*H_i_s.
```

At paired comoving equilibrium (`n_s = n_eq_s`, `J_s = J_eq_s`, and `H_i = 0`)
these number and E/F sources vanish. The E/F solve uses the same implicit
substep schedule as the ordinary local source route. Provider rates, the metric,
and fluid primitives remain frozen across the pointwise call. The partner
occupancy is frozen within each E/F subsolve; when fallback substepping advances
to another substep, the current states and partner comoving number are
recomputed.

## Pair opacity in transport

The inverse pair energy contribution depends on the partner state. For an
electron species `s`, with partner `p`, the host prepares the transport
coefficient

```text
kappa_transport_s = rates_s.kappa_tr
  + sum_c(eta_E_pair_s[c]/J_eq_s) * (n_p/n_eq_p).
```

The partner density uses the same current normalization as the number update,
`n_p = N_p/(W - H_n_p/J_p)`. The host adds this contribution while preparing
the face opacity; it must not mutate the validated rate bundle. Omitting it
would suppress the pair inverse reaction from the opacity-aware four-point
transport path. The current transport rule is in
[`PAIR_SOURCE_MODEL.md`](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md#opacity-supplied-to-transport)
and the [M1 integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md#electron-flavor-pair-source-update).

## Composition and publication

The paired call first applies the independent charged-current/scattering source
to temporary states, then applies the pair reaction. Pair increments are equal
and have weights `+1` and `-1`, so pair reactions contribute no `Y_e` change.
Charged-current lepton exchange from the first stage is retained. Both species
must pass validation, endpoint bounds, source solves, and exchange assembly
before either output is published. A failure restores both source bases and
zeroes both exchange packets.

The public implementation is therefore not a pair of independent one-species
calls. Single-species source operations reject separated electron pair fields;
use [`ghl_m1_solve_neutrino_pair_source_update`](../../../GRHayL/include/ghl_m1.h)
for the coupled case.

## Approximation boundary

This is a first-order, grey, isotropic collision approximation with frozen
number-current normalizations and partner occupancies during subsolves. It is
not an energy-resolved or angle-resolved annihilation kernel and does not
provide exact spectral redistribution.

The table-free reference provider's pair model is a deterministic test model,
not a calibrated weak-interaction prescription. The table-backed path retains
the production provider's raw pair/plasmon/bremsstrahlung emissivities, but the
M1 pair operator still applies the grey model above.

## Evidence

The coupled pair tests in
[`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c)
check equal number increments, an independent quadratic oracle, inactive-pair
delegation, endpoint bounds, and transactional failures. The current source
and pair ownership boundary is also summarized in
[`TRACEABILITY.md`](../../../GRHayL/Radiation/TRACEABILITY.md#production-boundary).

## Historical status

- **Current:** process-indexed electron pair fields, shared number extent,
  separate energy weights, partner-dependent transport opacity, and zero pair
  `Y_e` exchange.
- **Adaptable context:** the neutrino interaction whitepaper's motivation for
  pair annihilation, plasmon decay, and bremsstrahlung as number/energy source
  channels.
- **Superseded:** aggregate electron pair rates passed to a one-species source
  call, or any design that clips the two number populations independently.
- **Future/non-claim:** exact inverse pair kinetics, spectral pair transport,
  and fully coupled matter/radiation pair equilibration.
