# Neutrino lepton and matter exchange

Radiation and matter exchange different conserved quantities. The current
implementation keeps total radiation-number change, charged-current
electron-lepton change, energy/momentum exchange, and `Y_e` bookkeeping as
separate fields in the exchange packet. The public definitions are in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h).

## Signed species weights

The validated species weights are

```text
nu_e       : lepton_weight = +1
anti-nu_e  : lepton_weight = -1
nu_x       : lepton_weight =  0.
```

For a charged-current number increment `dN_cc,s`, the signed radiation
electron-lepton increment is

```text
dL_rad_cc,s = lepton_weight_s * dN_cc,s.
```

The sign belongs to the species, not to the absolute number field. Validation
rejects a bundle whose weight does not match its species; see
[`ghl_m1_neutrino_rates.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rates.c).

For the ordinary backward-Euler source path, the charged-current increment is
the difference between the physical, un-repaired number endpoint and its source
base:

```text
dN_cc = N_physical_endpoint - N_initial.
```

For validated single-species charged-current rates, this equals the endpoint
source expression in exact arithmetic. Subtracting the number endpoints avoids
cancellation in the stiff limit. The optional mean-energy number projection
instead uses `alpha*dt * (eta_N_cc -
kappa_a_N_cc*N_physical_endpoint/Gamma_N_physical_endpoint)`. The heavy-flavor
charged-current increment is zero. The source then applies the validated signed
weight. This increment excludes the subsequent pair, plasmon, or
bremsstrahlung number change. The branch is implemented in
[`ghl_m1_neutrino_lepton_increment.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_lepton_increment.c),
called by
[`ghl_m1_neutrino_implicit_solve.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit_solve.c)
and the compatibility dispatcher in
[`ghl_m1_neutrino_source_update.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_source_update.c).

## Electron-fraction update

The canonical matter-composition recommendation is

```text
Delta Y_e,matter = - Delta L_rad_cc / n_b,cons.
```

The required `n_b,cons` is the undensitized Eulerian conserved baryon number
density, `W*rho/m_b`, in compatible units. A host that stores a densitized
baryon variable converts it before calling the helper. The sign consequences
are:

| radiation event | `Delta L_rad_cc` | matter `Delta Y_e` |
| --- | ---: | ---: |
| emit `nu_e` | positive | negative |
| absorb `nu_e` | negative | positive |
| emit `anti-nu_e` | negative | positive |
| absorb `anti-nu_e` | positive | negative |
| `nu_x` | zero | zero |

The implementation validates a positive baryon normalization and publishes the
composition increment only after the signed lepton input is valid in
[`ghl_m1_neutrino_lepton_increment.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_lepton_increment.c).

The default source policy uses charged-current exchange. An opt-in policy can
derive a signed total-number composition increment, but that is a caller
selection and does not redefine the default contract. The policy is declared
in [`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h) and applied in
[`ghl_m1_neutrino_source_update.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_source_update.c).

## Energy and momentum conservation

For a local source update, the radiation packet records undensitized changes

```text
dE_rad = E_out - E_base,
dF_rad = F_out - F_base.
```

The matter packet is the equal-and-opposite densitized increment:

```text
dTau_matter = -sqrt(det(gamma)) * dE_rad,
dS_matter_i = -sqrt(det(gamma)) * dF_rad_i.
```

The source projection that supplies the radiation E/F change carries the
coordinate-time factor and metric volume in the host/source stage. The exchange
assembly itself applies the matter negation exactly once; see
[`ghl_m1_neutrino_exchange.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_exchange.c)
and the [matter-coupling helper](../../../GRHayL/Radiation/ghl_m1_matter_coupling_sources.c).

`dN_rad_total` is the total radiation-number difference between the output and
source-base states. It is separate from `dL_rad_cc`, which is the signed
charged-current subset. Pair reactions can change total radiation number while
leaving `Y_e` unchanged because their electron and antielectron number
increments are equal.

## Three-species aggregation

The host accumulates the two electron-flavor signed contributions and the
zero-weight heavy-flavor contribution across a coupled matter update. It must
not apply `nu_x`'s multiplicity a second time and must not count pair number as
charged-current lepton number. The paired source contract explicitly retains
the independent charged-current packet and reports no pair `Y_e` increment;
see [`PAIR_SOURCE_MODEL.md`](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md#composition-exchange-and-failure).

The host also owns final matter recovery and publication. For the normal
three-species stage, all species are solved into temporary states and exchange
packets before the common admissible limiter is applied. This prevents a
species-by-species limiter from breaking the energy, momentum, or lepton
bookkeeping. The ownership and common-limiter rule are in
[`M1_INTEGRATION_CONTRACT.md`](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md#coupled-limiter).

## Failure and conservation boundary

Exchange packets are transactional. A failed local source update leaves its
output at the transport/source base and returns a zero exchange packet; a
paired failure leaves both species unpublished. No repair-floor injection is
silently converted into a physical matter source. The host may track input
repair separately and owns any `Y_e` bounds or final coupled limiter.

## Historical status

- **Current:** signed weights, charged-current-only default `Y_e` exchange,
  separate total-number diagnostics, equal-and-opposite energy/momentum
  packets, and pair cancellation in electron lepton number.
- **Adaptable context:** the whitepapers' conservation argument and sign
  motivation for electron-neutrino versus antineutrino exchange.
- **Superseded:** deriving `Y_e` from an unsigned total number increment, using
  a densitized baryon normalization without conversion, or publishing each
  species independently before the common limiter.
- **Future/non-claim:** a fully coupled matter EOS/Con2Prim solve inside the
  local Radiation kernel; that remains host-owned.

## Evidence

The source-update tests check the lepton identity and exchange signs in
[`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c),
while seeded invariant tests exercise source, matter-coupling, number, and
lepton paths in
[`unit_test_m1_neutrino_seeded_invariants.c`](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c).
