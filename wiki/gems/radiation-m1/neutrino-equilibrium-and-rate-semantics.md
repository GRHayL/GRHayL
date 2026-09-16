# Neutrino equilibrium and rate semantics

The M1 kernels consume a frozen, per-species rate bundle. The provider's
equilibrium targets and the Radiation source equations are related, but they
are different responsibilities: the provider computes microphysics and
targets, while Radiation validates and applies the grey operator.

## Provider-side weak equilibrium

For the deterministic reference provider, the weak-equilibrium chemical
potential convention is

```text
mu_nu_e      = mu_e + mu_p - mu_n
mu_anti_nu_e = -mu_nu_e
mu_nu_x      = 0.
```

The reference backend converts these degeneracies into its grey equilibrium
targets. The table-backed NRPyLeakage path receives the EOS chemical-potential
state and uses the corresponding electron-flavor degeneracy convention in its
raw-rate adapter. This is provider-side behavior in
[`ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c#L599-L627)
and [`ghl_m1_nrpyleakage_kernel.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_nrpyleakage_kernel.c#L100-L203),
not a Radiation callback that recomputes chemical equilibrium during a source
solve.

`n_eq` and `J_eq` are the provider's comoving grey equilibrium number and
energy targets. They are not the photon-LTE expression `a_R T^4`. The public
bundle also carries a positive `mean_energy` satisfying

```text
J_eq = n_eq * mean_energy.
```

The exact public field meanings are documented in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h#L930-L965), while provider
ownership and backend selection are in
[`ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h#L14-L27).

The chemical-potential relation is useful physical context from the neutrino
whitepaper. The current interface does not expose a neutrino chemical
potential field; downstream code should use the validated targets rather than
re-deriving them from an assumed EOS convention.

## Separate number and energy data

The bundle carries separate number and energy coefficients:

```text
eta_N,    kappa_a_N       number emissivity and absorption opacity
eta_E,    kappa_a_E       energy emissivity and absorption opacity
eta_N_cc, kappa_a_N_cc    charged-current number subset
kappa_s                  isoenergetic scattering opacity
kappa_tr                 transport opacity.
```

The aggregate one-body identities enforced by the current validator are

```text
eta_N    = kappa_a_N    n_eq,
eta_E    = kappa_a_E    J_eq,
eta_N_cc = kappa_a_N_cc n_eq,
kappa_tr = kappa_a_E + kappa_s.
```

The comparisons use the repository's binary64 validation rules; a product that
underflows to exactly zero is reported diagnostically rather than replaced by
an arbitrary floor. See
[`ghl_m1_neutrino_rates.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_rates.c#L49-L181).

Number and energy therefore do not share a hidden fixed mean-energy closure.
The number source uses `kappa_a_N` and the endpoint number-current
normalization, while the E/F source uses `kappa_a_E`, `kappa_s`, `J`, and `H`.
For a single species, the current grey source is

```text
N_source = eta_N - kappa_a_N * N / Gamma_N,
Q        = eta_E - kappa_a_E * J,
S_E      = Q W + kappa_tr H_n,
S_i      = Q W V_i - kappa_tr H_i.
```

The number update uses the endpoint normalization in its backward-Euler form,

```text
N_new = (N_old + alpha*dt*eta_N)
        / (1 + alpha*dt*kappa_a_N/Gamma_N_new).
```

`Gamma_N` is the radiation number-current normalization, not the fluid Lorentz
factor. The equations and the endpoint update are implemented in
[`ghl_m1_neutrino_sources.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_sources.c#L6-L22)
and [`ghl_m1_neutrino_number_flux.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_number_flux.c#L15-L78).

## Pair fields are not ordinary one-species rates

`eta_N_pair[c]` and `eta_E_pair[c]` are process-resolved emissivities for the
electron-flavor pair operation. They are not to be added independently to a
single-species source call. The current validator rejects an electron bundle
with nonzero pair fields when the partner is unavailable. The joint operation
uses the two electron states and the partner occupancy; see the
[pair-source leaf](neutrino-pair-source-model.md) and the
[pair source contract](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md).

The separate fields are physically important: the emitted number and energy
averages may differ. The implementation deliberately does not force

```text
eta_E_pair / eta_N_pair = J_eq / n_eq.
```

The number emissivity is shared by the electron pair, while each species keeps
its raw energy emissivity. The provider mapping preserves this distinction in
[`ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c#L561-L595).

## What equilibrium means here

At a grey equilibrium fixed point for an independent one-body source,

```text
J = J_eq,  N/Gamma_N = n_eq,  H_i = 0,
```

and the corresponding collision source vanishes. This is a relaxation
statement about the supplied grey targets. It is not a claim of a resolved
Fermi-Dirac spectrum, exact angular equilibrium, or exact detailed balance for
every unresolved reaction channel.

For paired electron flavors, equilibrium is joint: the number reaction
vanishes when both partner occupancies match their targets, and the E/F source
vanishes when the pair energy moments and comoving fluxes satisfy the paired
grey conditions. The paired equations are specified in
[`PAIR_SOURCE_MODEL.md`](../../../GRHayL/Radiation/PAIR_SOURCE_MODEL.md#shared-number-reaction).

## Historical whitepaper status

- **Current:** provider-supplied weak-equilibrium targets, separate number and
  energy identities, endpoint-`Gamma_N` number update, and pair fields kept
  separate for electron flavors.
- **Adaptable context:** the whitepaper explanation of beta equilibrium and
  Kirchhoff consistency. It is explanatory; current validation and source code
  decide the exact accepted fields.
- **Superseded:** photon `J_eq = a_R T^4`, callback-owned photon opacity, and
  Phase 1's reduced number source/current conventions.
- **Future/non-claim:** group-resolved spectral equilibrium, inelastic
  redistribution, and an exact coupled matter/radiation microphysics solve.

## Evidence

Use [`unit_test_m1_rate_provider.c`](../../../Unit_Tests/unit_test_m1_rate_provider.c#L880-L1035)
for equilibrium/rate construction and
[`unit_test_m1_neutrino_source_update.c`](../../../Unit_Tests/unit_test_m1_neutrino_source_update.c#L1011-L1288)
for independent pair source oracles and equilibrium behavior. The scoped test
entry point is [`run_m1_tests.sh`](../../../Unit_Tests/run_m1_tests.sh).
