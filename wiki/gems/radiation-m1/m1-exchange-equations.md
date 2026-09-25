# M1 radiation and matter exchange equations

Radiation source updates return a local exchange packet so the host can apply
energy, momentum, and composition changes consistently. This leaf distinguishes
instantaneous source right-hand sides from integrated state increments. The
current assembly is
[ghl_m1_neutrino_exchange.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_exchange.c),
with public fields in [ghl_m1.h](../../../GRHayL/include/ghl_m1.h).

## Radiation increments

All local source states are undensitized. Let `state_base` be the state against
which the source update is measured and `state_out` the accepted endpoint:

$$
\Delta N_{\rm rad}=N_{\rm out}-N_{\rm base},
$$

$$
\Delta E_{\rm rad}=E_{\rm out}-E_{\rm base},
\qquad
\Delta F_{i,\rm rad}=F_{i,\rm out}-F_{i,\rm base}.
$$

For the ordinary implicit solver, `state_base` is its `state_in`; for the
source-policy dispatcher, the transport-predicted `state_transport` is the
source base and the pre-transport `state_input` is retained for validation.
The dispatcher therefore reports increments from transport, not from the
pre-transport stage.

Because the host’s conservative radiation variables contain `sqrt_detgamma`,
the corresponding integrated conservative radiation increments are

$$
\Delta\widetilde E_{\rm rad}=\sqrt{\gamma}\,\Delta E_{\rm rad},
\qquad
\Delta\widetilde F_{i,\rm rad}=\sqrt{\gamma}\,\Delta F_{i,\rm rad}.
$$

The number increment is also available as a diagnostic. It is not necessarily
the charged-current lepton increment when pair, plasmon, bremsstrahlung, or
other non-CC number channels are present.

## Energy and momentum conservation

The local source equations use the radiation-side projections `S_E` and `S_i`.
Their instantaneous conservative RHS contributions are

$$
\left(\partial_t\widetilde E\right)_{\rm int}
=\alpha\sqrt{\gamma}\,S_E,
\qquad
\left(\partial_t\widetilde F_i\right)_{\rm int}
=\alpha\sqrt{\gamma}\,S_i.
$$

The accepted integrated exchange packet instead uses the final state
differences above. The equal-and-opposite matter conservative increments are

$$
\Delta\widetilde\tau_{\rm matter}
=-\sqrt{\gamma}\,\Delta E_{\rm rad},
$$

$$
\Delta\widetilde S_{i,\rm matter}
=-\sqrt{\gamma}\,\Delta F_{i,\rm rad}.
$$

There is no additional lapse or timestep in these increment formulas: those
factors have already entered the accepted radiation endpoint. If a host is
assembling instantaneous RHS terms instead, it uses the explicit
`alpha*sqrt_detgamma` factors shown above. The distinction is implemented by
[ghl_m1_neutrino_assemble_exchange](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_exchange.c)
and [the matter coupling helper](../../../GRHayL/Radiation/ghl_m1_matter_coupling_sources.c).

## Charged-current lepton exchange

For an ordinary backward-Euler electron-flavor update, the charged-current
number increment uses the physical, un-repaired number endpoint and the source
base:

$$
\Delta N_{\rm cc}=N_{\rm phys,out}-N_{\rm base}.
$$

For validated single-species charged-current rates, this is algebraically
equivalent to the endpoint source expression below. The implementation uses
the endpoint difference to avoid cancellation when the source is stiff. When
the optional mean-energy number projection is selected, it instead evaluates
that endpoint source expression with the physical projected number and its
current normalization:

$$
\Delta N_{\rm cc}=\Delta t_\alpha
\left(\eta_{N,\rm cc}-\kappa_{a,N,\rm cc}
\frac{N_{\rm phys,out}}{\Gamma_{N,\rm phys,out}}\right),
\qquad
\Delta t_\alpha=\alpha\,\Delta t.
$$

The heavy-flavor charged-current increment is zero. The branch selection is
implemented in
[ghl_m1_neutrino_lepton_increment.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_lepton_increment.c).

The species-signed radiation lepton increment is

$$
\Delta L_{\rm rad,cc}=w_s\,\Delta N_{\rm cc},
$$

with exact current species weights

$$
w_{\nu_e}=+1,
\qquad
w_{\bar\nu_e}=-1,
\qquad
w_{\nu_x}=0.
$$

The provider/rates validator owns the species-weight consistency. The default
source policy uses this CC subset for composition, so pair processes do not
produce an electron-fraction increment. The total radiation-number increment
`dN_rad_total` remains separately observable.

## Matter composition increment

Given the signed CC radiation lepton increment and the positive undensitized
Eulerian conserved baryon-number density,

$$
\Delta Y_{e,\rm matter}
=-\frac{\Delta L_{\rm rad,cc}}{n_{b,\rm cons}}.
$$

For the current host convention,

$$
n_{b,\rm cons}=\frac{W\rho}{m_b},
$$

before any host-specific unit normalization. A host that stores a densitized
baryon variable must convert it before calling the local helper. The minus sign
means that emitting \(\nu_e\) lowers matter `Y_e`, while absorbing \(\nu_e\)
raises it; antineutrino exchange has the opposite weight. The implementation
is [ghl_m1_neutrino_lepton_increment.c](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_lepton_increment.c).

An opt-in dispatcher policy can instead use the signed total-number increment,

$$
\Delta L_{\rm rad,total}=w_s\,\Delta N_{\rm rad},
$$

but this is a bookkeeping choice and must not be mixed with the default CC
packet for the same update. The source dispatcher applies the selected policy
after a successful endpoint and exchange assembly.

## Coupled limiting boundary

The local solver reports radiation increments, equal-and-opposite densitized
matter energy/momentum increments, and the composition recommendation. The
host owns the coupled admissibility limiter and matter recovery. When a common
scalar \(\theta\) is applied, it must scale the radiation and corresponding
matter/lepton increments together so the equal-and-opposite and signed-lepton
relationships remain exact at the limited endpoint. The local Radiation code
does not own grid loops, EOS/Con2Prim recovery, or publication of host matter
variables.

## Transactional publication

Exchange assembly computes a complete local candidate and writes the packet only
after all differences, normalization, and composition arithmetic are valid. A
source-update failure initializes its output state to the source base and its
exchange packet to zero. This includes invalid inputs, rejected closure
fallbacks, failed Newton trials, endpoint bound failures, and disallowed double
application of interaction sources. An exhausted retry schedule returns the
distinct terminal no-update status with the same base/zero publication.

The direct matter-coupling helper has the same atomic publication boundary for
its output arguments. After required pointer validation, every non-success
return leaves `source_tildetau` and all three `source_tildeS` components
unchanged. It computes and validates the energy and momentum candidates in
locals first, then publishes all four values only when they are finite. A
successful call retains the projections
`-alpha*sqrt_detgamma*S_E` and `-alpha*sqrt_detgamma*S[i]`.

The paired `{nue, anue}` source operation extends the transaction across both
species: either both final states and packets publish or both source bases and
zero packets remain. Pair number increments are equal and therefore preserve
the pair number difference; pair reactions contribute no `dYe` directly. See
[PAIR_SOURCE_MODEL.md](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md#composition-exchange-and-failure)
for the paired collision details.

## Photon and obsolete-method boundary

The equal-and-opposite projection is shared radiation/matter mathematics. The
photon whitepaper’s LTE energy target and photon opacity prescription are not
needed to form a neutrino exchange packet. Likewise, an HLL or reduced-number
current formula from an earlier design does not change the current increment,
lepton-weight, or transaction rules.

## Focused evidence

- [Neutrino source-update tests](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c)
  covers endpoint exchange, lepton signs, pair publication, terminal no-update,
  and transactional failures.
- [Seeded M1 invariants](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
  covers source projections, exchange, and species invariants.
- [M1 error handling](../../../Unit_Tests/unit_test_m1_error_handling.c)
  covers unchanged outputs on invalid source/exchange inputs.
- [Neutrino M1 contract](../../../wiki/gems/radiation-m1/neutrino-m1-contract.md#transactional-source-update)
  records the host/library ownership and publication boundary.
